configfile: "config.yaml"

import glob
import re
from collections import defaultdict

# set up lists with all files and their folders
raw_data=config["raw_data"]
samples = []
merged_samples = []
folders = {}
replicates = defaultdict(list)
print(raw_data)
for file in glob.glob(raw_data + "/**", recursive=True):
    if '_R1.fq.gz' in file:
        split = re.split('/|_R1', file)
        filename, folder = split[-2], split[-3]
        split = re.split('_', filename)
        merged_filename, replicate = split[-2], split[-1]
        if folder in config["sets"]: # only run files that are defined in config.yaml
            folders[filename] = folder  # we will need this one later
            samples.append(filename)
            replicates[merged_filename].append(replicate)
            merged_samples.append(merged_filename)

onstart:
    shell(expand("mkdir {output}/logs", output=config["output"]))

rule all:
    input:
        expand("{output}/{sample}/s6_{resolution}/{sample}.validPairs.binPair.uniMulti",
        output=config["output"], resolution=config["resolution"], sample=samples)

rule generate_genome_file_without_unitigs:
    input:
        config["genome"]
    output:
        expand("{output}/genome_files/{genome}_without_unitigs.fa", output=config["output"], genome=config["genomeName"])
    shell:
        "/work/project/ladsie_002/seqkit seq -w 0 {input} | grep -A 1 -E \"^>BES|^>Chr\" > {output}"

rule generate_genome_size_file:
    input:
        expand("{output}/genome_files/{genome}_without_unitigs.fa", output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/genome_files/{genome}.sizes", output=config["output"], genome=config["genomeName"])
    shell:
        """
        python bin/generate_genome_size_file.py {input} > {output}
        """

rule generate_digestion_file:
    input:
        expand("{output}/genome_files/{genome}_without_unitigs.fa", output=config["output"], genome=config["genomeName"])
    params:
        restriction_site="^GATC"
    output:
        expand("{output}/genome_files/{genome}_{enzyme}.bed", output=config["output"], genome=config["genomeName"], enzyme=config["enzyme"])
    shell:
        """
        set +u
        source /work/project/ladsie_002/myPython-2/bin/activate
        set -u
        python /work/project/ladsie_002/HiC-Pro_2.11.1/bin/utils/digest_genome.py \
        -r {params.restriction_site} -o {output} {input}
        """

rule bwa_index:
    input:
        expand("{output}/genome_files/{genome}_without_unitigs.fa",
            output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/genome_files/{genome}_without_unitigs.fa.{indices}",
            output=config["output"], genome=config["genomeName"], indices=["amb","ann","bwt","pac","sa"])
    shell:
        "bwa index {input}"


## ******************
## step 1: Alignment
## ******************
rule alignment:
    input:
        read=lambda wildcards: expand("{data}/{set}/{sample}_R{i}.fq.gz",
            data=config["raw_data"], set=folders[wildcards.sample], sample=wildcards.sample, i={1,2}),
        indices=expand("{output}/genome_files/{genome}_without_unitigs.fa.{indices}",
            output=config["output"], genome=config["genomeName"], indices=["amb","ann","bwt","pac","sa"]),
        genome=expand("{output}/genome_files/{genome}_without_unitigs.fa",
            output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/{{sample}}/s1/{{sample}}_{i}.bam",
            output=config["output"], i={1,2})
    params:
        read_folder=lambda wildcards: expand("{data}/{set}",
            data=config["raw_data"], set=folders[wildcards.sample]),
        out_folder=lambda wildcards: expand("{output}/{sample}",
            output=config["output"], sample=wildcards.sample),
        cutsite=config["enzyme_cutsite"]
    threads: 16
    shell:
        """
        name={wildcards.sample}
        ref={input.genome}
        bwaDir=$(dirname $(which bwa))
        samtoolsDir=$(dirname $(which samtools))
        fastqDir={params.read_folder}
        resultsDir={params.out_folder}/s1
        bin=$PWD/bin
        nCores={threads}
        saveFiles=0
        seqLength=25
        cutsite={params.cutsite}

        ## alignment
        echo "Start step 1 - alignment!"
        bash s1_bwaAlignment.sh \
            "$name" \
            "$ref" \
            "$bwaDir" \
            "$samtoolsDir" \
            "$fastqDir" \
            "$resultsDir" \
            "$bin" \
            "$nCores" \
            "$resultsDir/mHiC.summary_s1" \
            "$saveFiles"  \
            "$seqLength" \
            "${{cutsite[@]}}"
        """

## **************************
## step 2: Read ends pairing
## **************************
rule read_ends_pairing:
    input:
        read1=expand("{output}/{{sample}}/s1/{{sample}}_1.bam",
            output=config["output"]),
        read2=expand("{output}/{{sample}}/s1/{{sample}}_2.bam",
            output=config["output"])
    output:
        expand("{output}/{{sample}}/s2/{{sample}}.bam",
            output=config["output"])
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}",
            output=config["output"], sample=wildcards.sample)
    shell:
        """
        name={wildcards.sample}
        resultsDir={params.out_folder}

        echo "Start step 2 - joining read ends!"
        python s2_joinEnd.py \
            -r1 {input.read1} \
            -r2 {input.read2} \
            -o $resultsDir/s2/$name.bam \
            -sf $resultsDir/mHiC.summary_s2
        """

## *********************************
## step 3: Valid fragment filtering
## *********************************
rule valid_fragment_filtering:
    input:
        sam=expand("{output}/{{sample}}/s2/{{sample}}.bam",
            output=config["output"]),
        refrag=expand("{output}/genome_files/{genome}_{enzyme}.bed",
            output=config["output"], genome=config["genomeName"], enzyme=config["enzyme"])
    output:
        expand("{output}/{{sample}}/s3_{{resolution}}/{{sample}}.uniMulti",
            output=config["output"]),
        expand("{output}/{{sample}}/s3_{{resolution}}/{{sample}}.validPairs",
            output=config["output"])
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}",
            output=config["output"], sample=wildcards.sample),
        resolution=lambda wildcards: wildcards.resolution
    shell:
        """
        name={wildcards.sample}
        resultsDir={params.out_folder}
        resolution={params.resolution}
        bin=$PWD/bin
        refrag={input.refrag}
        lowerBound=$((resolution * 2))
        refragL=50 #$((seqLength * 2))
        refragU=500

        echo "Start step 3 - categorize read pairs!"
        python s3_categorizePairs.py \
            -f ${{refrag}} \
            -r ${{resultsDir}}/s2/${{name}}.sam \
            -o ${{resultsDir}}/s3_${{resolution}} \
            -l $refragL \
            -u $refragU \
            -d $lowerBound \
            -m "window" \
            -b $resolution \
            -sf $resultsDir/mHiC.summary_w${{resolution}}_s3
        """

## ***************************************
## step 4 - Remove duplicates and binning.
## ***************************************
rule duplicates_removal:
    input:
        pairs=expand("{output}/{{sample}}/s3_{{resolution}}/{{sample}}.validPairs",
            output=config["output"]),
        sizes=expand("{output}/genome_files/{genome}.sizes",
            output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/{{sample}}/s4_{{resolution}}/{{sample}}.validPairs.binPair.Marginal",
            output=config["output"]),
        expand("{output}/{{sample}}/s4_{{resolution}}/{{sample}}.validPairs.binPairCount.uni",
            output=config["output"]),
        expand("{output}/{{sample}}/s4_{{resolution}}/{{sample}}.validPairs.MULTI.binPair.multi",
            output=config["output"])
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}",
            output=config["output"], sample=wildcards.sample),
        resolution=lambda wildcards: wildcards.resolution,
        chr_list=config["chr_list"]
    shell:
        """
        name={wildcards.sample}
        resultsDir={params.out_folder}
        resolution={params.resolution}
        bin=$PWD/bin
        validP={input.pairs}
        validI={output}
        minCount=1 #min contact counts allowed

        normMethod="KR" #1. "ICE" 2. "KR" 3."None"
        ICEmappFile="" ## mappability file for ICE method
        ICEminMap="" ## min mappability threshold for ICE method
        ICEmaxIter="" ## maximum number of iteration for ICE method
        KRchromSizeFile="{input.sizes}" ## chromosome size file for KR method
        KRsparsePerc=5 ## remove *% of sparse regions for KR method
        splitByChrom=0
        saveSplitContact=0
        chrList=({params.chr_list})


        echo "Start step 4 - duplicates removal and binning!"
        bash s4_bin.sh \
            "$validP" \
            "$validI" \
            "$bin" \
            "$resolution" \
            "$minCount" \
            "$normMethod" \
            "$ICEmappFile" \
            "$ICEminMap" \
            "$ICEmaxIter" \
            "whole" \
            "$KRchromSizeFile" \
            "$KRsparsePerc" \
            "$resultsDir/mHiC.summary_w${{resolution}}_s4" \
            "$splitByChrom" \
            "$saveSplitContact" \
            "${{chrList[@]}}"
        """

## **********************
## step 5 - Build prior.
## **********************
rule generative_model_prior:
    input:
        pairs=expand("{output}/{{sample}}/s4_{{resolution}}/{{sample}}.validPairs.binPair.Marginal",
            output=config["output"]),
        sizes=expand("{output}/genome_files/{genome}.sizes",
            output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/{{sample}}/s5_{{resolution}}/s5_prior.mhic",
            output=config["output"])
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}",
            output=config["output"], sample=wildcards.sample),
        resolution=lambda wildcards: wildcards.resolution
    shell:
        """
        name={wildcards.sample}
        resultsDir={params.out_folder}
        resolution={params.resolution}
        bin=$PWD/bin
        validI=$(echo {input.pairs} | sed 's/.binPair.Marginal//')
        splineBin=150
        lower=$((resolution * 2))
        priorName="uniPrior"
        normMethod="KR" #"ICE" ##"KR"
        chromSizeFile={input.sizes}
        contactFile=$validI.binPairCount.uni

        echo "Starts step 5 - prior construction based on uni-reads only!"
        python $bin/createFitHiCFragments-fixedsize.py \
            --chrLens "$chromSizeFile" \
            --resolution "$resolution" \
            --outFile "$resultsDir/s5_${{resolution}}/$name.$resolution.uni.fragments.mHiC"

        python s5_prior.py \
            -f $resultsDir/s5_${{resolution}}/$name.$resolution.uni.fragments.mHiC \
            -i $contactFile \
            -o $resultsDir/s5_${{resolution}} \
            -t $validI.binPairCount.uni.KRnorm.bias \
            -b $splineBin \
            -L $lower \
            -r $resolution \
            -p 2
        """

## ************************************************************************************
## step 6 - Generative model to assign probability to multi-reads potential alignments.
## ************************************************************************************
rule assign_multi_reads:
    input:
        multi=expand("{output}/{{sample}}/s4_{{resolution}}/{{sample}}.validPairs.MULTI.binPair.multi",
            output=config["output"]),
        uni=expand("{output}/{{sample}}/s4_{{resolution}}/{{sample}}.validPairs.binPairCount.uni",
            output=config["output"]),
        prior=expand("{output}/{{sample}}/s5_{{resolution}}/s5_prior.mhic",
            output=config["output"]),
        sizes=expand("{output}/genome_files/{genome}.sizes",
            output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/{{sample}}/s6_{{resolution}}/{{sample}}.validPairs.binPair.multi.mHiC",
            output=config["output"])
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}",
            output=config["output"], sample=wildcards.sample),
        resolution=lambda wildcards: wildcards.resolution
    shell:
        """
        name={wildcards.sample}
        resultsDir={params.out_folder}
        resolution={params.resolution}
        prior={input.prior}
        multi={input.multi}
        multiKeys={input.multi}Keys
        uni={input.uni}
        filename="$name.validPairs.binPair.multi"
        threshold=0.5

        echo "Starts step 6 - assign probability to multi-reads potential alignment positions !"
        # awk -v OFS="=" '{{print $2, $3, $4, $5}}' $multi | sort -u >$multiKeys
        python s6_em.py -p $prior -u $uni -m $multi -mk $multiKeys -t $threshold -o "$resultsDir/s6_${{resolution}}" -f $filename
        """

## ************************************************************************************
## step 7 - Post-mHiC processing.
## ************************************************************************************
rule post_processing:
    input:
        multi=expand("{output}/{{sample}}/s4_{{resolution}}/{{sample}}.validPairs.MULTI.binPair.multi",
            output=config["output"]),
        uni=expand("{output}/{{sample}}/s4_{{resolution}}/{{sample}}.validPairs.binPairCount.uni",
            output=config["output"]),
        s6=expand("{output}/{{sample}}/s6_{{resolution}}/{{sample}}.validPairs.binPair.multi.mHiC",
            output=config["output"])
    output:
        expand("{output}/{{sample}}/s6_{{resolution}}/{{sample}}.validPairs.binPair.uniMulti",
            output=config["output"])
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}",
            output=config["output"], sample=wildcards.sample),
        resolution=lambda wildcards: wildcards.resolution
    shell:
        """
        name={wildcards.sample}
        resultsDir={params.out_folder}
        resolution={params.resolution}
        filterT=0.5
        multiOut={input.s6}
        multi={input.multi}
        multiKeys={input.multi}Keys
        uni={input.uni}

        awk -v OFS="\t" -v fT=$filterT '$6>fT {{print $2, $3, $4, $5}}' $multiOut | \
        sort | \
        uniq -c | \
        awk -v OFS="\t" '{{print $2, $3, $4, $5, $1}}' > \
        $multiOut.binPairCount.multi ## get binPair Count for multi-reads

        cat $uni $multiOut.binPairCount.multi | \
        sort -k1,1V -k2,2n -k3,3V -k4,4n | \
        awk -v OFS="\t" '{{a[$1" "$2" "$3" "$4]+=$5}}END{{for (i in a) print i,a[i]}}' | \
        sort -k1,1V -k2,2n -k3,3V -k4,4n > \
        {output} ## merged with uni-reads binpair count
        """
