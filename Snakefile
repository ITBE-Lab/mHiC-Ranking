configfile: "config.yaml"

import os
import re
from collections import defaultdict

# set up lists with all files and their folders
samples = []
merged_samples = set()
folders = {}
replicates = defaultdict(list)
data_sets = config["reads"]
for data_set in data_sets:
    folder = data_sets[data_set]
    files = [f for f in os.listdir(folder) if re.search(r'^[^.].*_R.\.f*', f)]
    for file in files:
        if '_R1' in file:
            my_split = re.split('_R1', file)
            filename = my_split[-2]
            my_split = re.split('_', filename)
            merged_filename, replicate = my_split[-2], my_split[-1]
            folders[filename] = folder  # we will need this one later
            samples.append(filename)
            replicates[merged_filename].append(replicate)
            merged_samples.add(merged_filename)
merged_samples=list(merged_samples)

gff_files = []
folder = config["gff_folder"]
files = [f for f in os.listdir(folder) if re.search(r'^[^.].*\.gff?', f)]
for file in files:
    my_split = re.split('.gff', file)
    filename = my_split[-2]
    gff_files.append(filename)


rule all:
    input:
        # expand("{output}/{sample}_merged/s7_{resolution}_4C/{sample}.interactions.homer",
        #         output=config["output"], resolution=config["resolution"], sample=merged_samples)
        expand("{output}/{sample}_merged/s7_{resolution}_{threshold}_4C/Virtual_4C/{gff_file}.wig",
           output=config["output"], resolution=config["resolution"], sample=merged_samples,
           threshold=config["probability"], gff_file=gff_files)
        # expand("{output}/{sample}_merged/s4_{resolution}/{genome}_without_unitigs_mappability.ice",
        #     resolution=config["resolution"], sample=merged_samples, output=config["output"], genome=config["genomeName"])

rule generate_genome_file_without_unitigs:
    input:
        config["genome"]
    output:
        expand("{output}/genome_files/{genome}_without_unitigs.fa", output=config["output"], genome=config["genomeName"])
    conda:
        "conda_environments/mhic_env.yaml"
    shell:
        "seqkit seq -w 0 {input} | grep -A 1 -E \"^>BES|^>Chr\" > {output}"

rule generate_genome_size_file:
    input:
        expand("{output}/genome_files/{genome}_without_unitigs.fa", output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/genome_files/{genome}.sizes", output=config["output"], genome=config["genomeName"])
    conda:
        "conda_environments/mhic_env.yaml"
    shell:
        "python bin/generate_genome_size_file.py {input} > {output}"

rule faidx:
    input:
        expand("{output}/genome_files/{genome}_without_unitigs.fa", output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/genome_files/{genome}_without_unitigs.fa.fai", output=config["output"], genome=config["genomeName"])
    conda:
        "conda_environments/mhic_env.yaml"
    shell:
        "samtools faidx {input}"

rule generate_digestion_file:
    input:
        expand("{output}/genome_files/{genome}_without_unitigs.fa", output=config["output"], genome=config["genomeName"])
    params:
        restriction_site="^GATC"
    output:
        expand("{output}/genome_files/{genome}_{enzyme}.bed", output=config["output"], genome=config["genomeName"], enzyme=config["enzyme"])
    conda:
        "conda_environments/hic_pro.yaml"
    shell:
        """
        python bin/digest_genome.py -r {params.restriction_site} -o {output} {input}
        """

rule bwa_index:
    input:
        expand("{output}/genome_files/{genome}_without_unitigs.fa",
            output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/genome_files/{genome}_without_unitigs.fa.{indices}",
            output=config["output"], genome=config["genomeName"], indices=["amb","ann","bwt","pac","sa"])
    conda:
        "conda_environments/mhic_env.yaml"
    shell:
        "bwa index {input}"

rule compile_cutsite:
    input:
        expand("bin/cutsite_trimming_mHiC.cpp")
    output:
        expand("bin/cutsite_trimming_mHiC")
    conda:
        "conda_environments/hic_pro.yaml"
    shell:
        "g++ -std=c++0x -o {output} {input}"

## ******************
## step 1: Alignment
## ******************
rule alignment:
    input:
        read=lambda wildcards: expand("{set}/{sample}_R{i}.fq.gz",
            set=folders[wildcards.sample], sample=wildcards.sample, i={1,2}),
        indices=expand("{output}/genome_files/{genome}_without_unitigs.fa.{indices}",
            output=config["output"], genome=config["genomeName"], indices=["amb","ann","bwt","pac","sa"]),
        genome=expand("{output}/genome_files/{genome}_without_unitigs.fa",
            output=config["output"], genome=config["genomeName"]),
        cutsite=expand("bin/cutsite_trimming_mHiC")
    output:
        temp(expand("{output}/{{sample}}/s1/{{sample}}_{i}.bam",
            output=config["output"], i={1,2}))
    params:
        read_folder=lambda wildcards: expand("{set}",
             set=folders[wildcards.sample]),
        out_folder=lambda wildcards: expand("{output}/{sample}",
            output=config["output"], sample=wildcards.sample),
        cutsite=config["enzyme_cutsite"]
    conda:
        "conda_environments/mhic_env.yaml"
    threads: 64
    shell:
        """
        name={wildcards.sample}
        ref={input.genome}
        bwaDir=$(dirname $(which bwa))
        samtoolsDir=$(dirname $(which samtools))
        fastqDir={params.read_folder}
        resultsDir={params.out_folder}/s1
        bin=$(pwd)/bin
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
        temp(expand("{output}/{{sample}}/s2/{{sample}}.bam",
            output=config["output"]))
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}",
            output=config["output"], sample=wildcards.sample)
    conda:
        "conda_environments/mhic_env.yaml"
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
        bam=expand("{output}/{{sample}}/s2/{{sample}}.bam",
            output=config["output"]),
        refrag=expand("{output}/genome_files/{genome}_{enzyme}.bed",
            output=config["output"], genome=config["genomeName"], enzyme=config["enzyme"])
    output:
        temp(expand("{output}/{{sample}}/s3_{{resolution}}/{{sample}}.uniMulti",
            output=config["output"])),
        temp(expand("{output}/{{sample}}/s3_{{resolution}}/{{sample}}.validPairs",
            output=config["output"]))
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}",
            output=config["output"], sample=wildcards.sample),
        resolution=lambda wildcards: wildcards.resolution
    conda:
        "conda_environments/mhic_env.yaml"
    shell:
        """
        name={wildcards.sample}
        resultsDir={params.out_folder}
        resolution={params.resolution}
        refrag={input.refrag}
        lowerBound=$((resolution * 2))
        refragL=50 #$((seqLength * 2))
        refragU=500

        echo "Start step 3 - categorize read pairs!"
        python s3_categorizePairs.py \
            -f ${{refrag}} \
            -r {input.bam} \
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
        unimulti=expand("{output}/{{sample}}/s3_{{resolution}}/{{sample}}.uniMulti",
            output=config["output"])
    output:
        expand("{output}/{{sample}}/s4_{{resolution}}/{{sample}}.validPairs.MULTI",
            output=config["output"]),
        expand("{output}/{{sample}}/s4_{{resolution}}/{{sample}}.validPairs.MultiRMdupList",
            output=config["output"]),
        expand("{output}/{{sample}}/s4_{{resolution}}/{{sample}}.validPairs.UNI",
            output=config["output"]),
        expand("{output}/{{sample}}/s4_{{resolution}}/{{sample}}.validPairs.UNI.binPair",
            output=config["output"])
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}",
            output=config["output"], sample=wildcards.sample),
        resolution=lambda wildcards: wildcards.resolution,
        chr_list=config["chr_list"]
    conda:
        "conda_environments/mhic_env.yaml"
    shell:
        """
        name={wildcards.sample}
        bin=$(pwd)/bin
        validP={input.pairs}
        validI=$(echo {input.pairs} | sed s/s3/s4/)
        splitByChrom=0
        saveSplitContact=0
        chrList=({params.chr_list})

        echo "Start step 4 - duplicates removal and binning!"
        bash s4_bin.sh \
            "$validP" \
            "$validI" \
            "$bin" \
            "$splitByChrom" \
            "$saveSplitContact" \
            "${{chrList[@]}}"
        """

# merge_replicates
rule merge_multi:
    input:
        lambda wildcards: expand("{output}/{sample}_{lib}/s4_{resolution}/{sample}_{lib}.validPairs.MULTI",
            output=config["output"], sample=wildcards.merged_sample, lib=replicates[wildcards.merged_sample],
            resolution=wildcards.resolution)
    output:
        expand("{output}/{{merged_sample}}_merged/s4_{{resolution}}/{{merged_sample}}.validPairs.MULTI",
            output=config["output"])
    shell:
        "cat {input} > {output}"

rule merge_duplist:
    input:
        lambda wildcards: expand("{output}/{sample}_{lib}/s4_{resolution}/{sample}_{lib}.validPairs.MultiRMdupList",
            output=config["output"], sample=wildcards.merged_sample, lib=replicates[wildcards.merged_sample],
            resolution=wildcards.resolution)
    output:
        expand("{output}/{{merged_sample}}_merged/s4_{{resolution}}/{{merged_sample}}.validPairs.MultiRMdupList",
            output=config["output"])
    shell:
        "cat {input} > {output}"

rule merge_uni:
    input:
        lambda wildcards: expand("{output}/{sample}_{lib}/s4_{resolution}/{sample}_{lib}.validPairs.UNI",
            output=config["output"], sample=wildcards.merged_sample, lib=replicates[wildcards.merged_sample],
            resolution=wildcards.resolution)
    output:
        expand("{output}/{{merged_sample}}_merged/s4_{{resolution}}/{{merged_sample}}.validPairs.UNI",
            output=config["output"])
    shell:
        "cat {input} > {output}"

rule merge_binpair:
    input:
        lambda wildcards: expand("{output}/{sample}_{lib}/s4_{resolution}/{sample}_{lib}.validPairs.UNI.binPair",
            output=config["output"], sample=wildcards.merged_sample, lib=replicates[wildcards.merged_sample],
            resolution=wildcards.resolution)
    output:
        expand("{output}/{{merged_sample}}_merged/s4_{{resolution}}/{{merged_sample}}.validPairs.UNI.binPair",
            output=config["output"])
    shell:
        "cat {input} > {output}"



rule gem_mappability:
    input:
        expand("{output}/genome_files/{genome}_without_unitigs.fa",
            output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/genome_files/{genome}_without_unitigs_mappability.bg",
            output=config["output"], genome=config["genomeName"])
    params:
        kmer=76
    conda:
        "conda_environments/gem_env.yaml"
    threads: 64
    shell:
        """
        PREFIX=$(echo {output} | sed s/.bg//)
        gem-indexer -T {threads} -c dna -i {input} -o ${{PREFIX}}
        gem-mappability -T {threads} -I ${{PREFIX}}.gem -l {params.kmer} -o ${{PREFIX}}
        gem-2-bed mappability -I ${{PREFIX}}.gem -i ${{PREFIX}}.mappability -o ${{PREFIX}}
        """

rule ice_mappability:
    input:
        fai=expand("{output}/genome_files/{genome}_without_unitigs.fa.fai",
            output=config["output"], genome=config["genomeName"]),
        bed=expand("{output}/genome_files/{genome}_without_unitigs_mappability.bg",
            output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/{{sample}}_merged/s4_{{resolution}}/{genome}_without_unitigs_mappability.ice",
            output=config["output"], genome=config["genomeName"])
    params:
        resolution=lambda wildcards: wildcards.resolution
    conda:
        "conda_environments/r_env.yaml"
    script:
        "bin/bed2mappability.R"

## ***************************************
## step 4 - Normalization and binning.
## ***************************************
rule norm_bin:
    input:
        sizes=expand("{output}/genome_files/{genome}.sizes",
            output=config["output"], genome=config["genomeName"]),
        mappability=expand("{output}/{{sample}}_merged/s4_{{resolution}}/{genome}_without_unitigs_mappability.ice",
            output=config["output"], genome=config["genomeName"]),
        multi=expand("{output}/{{sample}}_merged/s4_{{resolution}}/{{sample}}.validPairs.MULTI",
            output=config["output"]),
        duplist=expand("{output}/{{sample}}_merged/s4_{{resolution}}/{{sample}}.validPairs.MultiRMdupList",
            output=config["output"]),
        uni=expand("{output}/{{sample}}_merged/s4_{{resolution}}/{{sample}}.validPairs.UNI",
            output=config["output"]),
        binpair=expand("{output}/{{sample}}_merged/s4_{{resolution}}/{{sample}}.validPairs.UNI.binPair",
            output=config["output"])
    output:
        expand("{output}/{{sample}}_merged/s4_{{resolution}}/{{sample}}.validPairs.binPair.Marginal",
            output=config["output"]),
        expand("{output}/{{sample}}_merged/s4_{{resolution}}/{{sample}}.validPairs.binPairCount.uni",
            output=config["output"]),
        expand("{output}/{{sample}}_merged/s4_{{resolution}}/{{sample}}.validPairs.MULTI.binPair.multi",
            output=config["output"])
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}_merged",
            output=config["output"], sample=wildcards.sample),
        resolution=lambda wildcards: wildcards.resolution,
        chr_list=config["chr_list"]
    conda:
        "conda_environments/mhic_env.yaml"
    shell:
        """
        name={wildcards.sample}
        resultsDir={params.out_folder}
        resolution={params.resolution}
        bin=$(pwd)/bin
        validP="ignored"
        validI=$(echo {input.multi} | sed s/.MULTI//)
        minCount=1 #min contact counts allowed

        normMethod="ICE" #1. "ICE" 2. "KR" 3."None"
        ICEmappFile={input.mappability} ## mappability file for ICE method
        ICEminMap=0.25 ## min mappability threshold for ICE method
        ICEmaxIter=100 ## maximum number of iteration for ICE method
        KRchromSizeFile="" ## chromosome size file for KR method
        KRsparsePerc=5 ## remove *% of sparse regions for KR method
        splitByChrom=0
        saveSplitContact=0
        chrList=({params.chr_list})


        echo "Start step 4 - duplicates removal and binning!"
        bash s4.2_normalization.sh \
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
        pairs=expand("{output}/{{sample}}_merged/s4_{{resolution}}/{{sample}}.validPairs.binPair.Marginal",
            output=config["output"]),
        sizes=expand("{output}/genome_files/{genome}.sizes",
            output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/{{sample}}_merged/s5_{{resolution}}/s5_prior.mhic",
            output=config["output"])
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}_merged",
            output=config["output"], sample=wildcards.sample),
        resolution=lambda wildcards: wildcards.resolution
    conda:
        "conda_environments/mhic_env.yaml"
    shell:
        """
        name={wildcards.sample}
        resultsDir={params.out_folder}
        resolution={params.resolution}
        validI=$(echo {input.pairs} | sed 's/.binPair.Marginal//')
        splineBin=200
        lower=$((resolution * 2))
        priorName="uniPrior"
        normMethod="ICE" #"ICE" ##"KR"
        chromSizeFile={input.sizes}
        contactFile=$validI.binPairCount.uni

        echo "Starts step 5 - prior construction based on uni-reads only!"
        python bin/createFitHiCFragments-fixedsize.py \
            --chrLens "$chromSizeFile" \
            --resolution "$resolution" \
            --outFile "$resultsDir/s5_${{resolution}}/$name.$resolution.uni.fragments.mHiC"

        python s5_prior.py \
            -f $resultsDir/s5_${{resolution}}/$name.$resolution.uni.fragments.mHiC \
            -i $contactFile \
            -o $resultsDir/s5_${{resolution}} \
            -t $validI.binPairCount.uni.ICEnorm.bias \
            -b $splineBin \
            -L $lower \
            -r $resolution \
            -p 1
        """

## ************************************************************************************
## step 6 - Generative model to assign probability to multi-reads potential alignments.
## ************************************************************************************
rule assign_multi_reads:
    input:
        multi=expand("{output}/{{sample}}_merged/s4_{{resolution}}/{{sample}}.validPairs.MULTI.binPair.multi",
            output=config["output"]),
        uni=expand("{output}/{{sample}}_merged/s4_{{resolution}}/{{sample}}.validPairs.binPairCount.uni",
            output=config["output"]),
        prior=expand("{output}/{{sample}}_merged/s5_{{resolution}}/s5_prior.mhic",
            output=config["output"]),
        sizes=expand("{output}/genome_files/{genome}.sizes",
            output=config["output"], genome=config["genomeName"])
    output:
        expand("{output}/{{sample}}_merged/s6_{{resolution}}/{{sample}}.validPairs.binPair.multi.mHiC",
            output=config["output"])
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}_merged",
            output=config["output"], sample=wildcards.sample),
        resolution=lambda wildcards: wildcards.resolution
    conda:
        "conda_environments/mhic_env.yaml"
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
        awk -v OFS="=" '{{print $2, $3, $4, $5}}' $multi | sort -u >$multiKeys
        python s6_em.py -p $prior -u $uni -m $multi -mk $multiKeys -t $threshold -o "$resultsDir/s6_${{resolution}}" -f $filename
        """

## ************************************************************************************
## step 6 - Post-mHiC processing.
## ************************************************************************************
rule post_processing:
    input:
        multi=expand("{output}/{{sample}}_merged/s4_{{resolution}}/{{sample}}.validPairs.MULTI.binPair.multi",
            output=config["output"]),
        uni=expand("{output}/{{sample}}_merged/s4_{{resolution}}/{{sample}}.validPairs.binPairCount.uni",
            output=config["output"]),
        s6=expand("{output}/{{sample}}_merged/s6_{{resolution}}/{{sample}}.validPairs.binPair.multi.mHiC",
            output=config["output"])
    output:
        expand("{output}/{{sample}}_merged/s6_{{resolution}}_{{threshold}}/{{sample}}.validPairs.binPair.uniMulti",
            output=config["output"])
    params:
        out_folder=lambda wildcards: expand("{output}/{sample}",
            output=config["output"], sample=wildcards.sample),
        resolution=lambda wildcards: wildcards.resolution,
        threshold=lambda wildcards: wildcards.threshold
    shell:
        """
        name={wildcards.sample}
        resultsDir={params.out_folder}
        resolution={params.resolution}
        filterT={params.threshold} # Filter the mHi-C outcome by the posterior probability
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

rule normalize_for_ploidy:
    input:
        expand("{output}/{{sample}}_merged/s6_{{resolution}}_{{threshold}}/{{sample}}.validPairs.binPair.uniMulti",
            output=config["output"])
    output:
        expand("{output}/{{sample}}_merged/s7_{{resolution}}_{{threshold}}_4C/{{sample}}.ploidy_normalized.uniMulti",
            output=config["output"])
    conda:
        "conda_environments/r_env.yaml"
    script:
        "bin/normalize_ploidy.R"

rule convert_to_homer:
    input:
        expand("{output}/{{sample}}_merged/s7_{{resolution}}_{{threshold}}_4C/{{sample}}.ploidy_normalized.uniMulti",
            output=config["output"])
    output:
        expand("{output}/{{sample}}_merged/s7_{{resolution}}_{{threshold}}_4C/{{sample}}.interactions.homer",
            output=config["output"])
    conda:
        "conda_environments/hicsd_env.yaml"
    shell:
        "hicsd mHiC2homer -m {input} -o {output}"

rule virtual_4c:
    input:
        matrix=expand("{output}/{{sample}}_merged/s7_{{resolution}}_{{threshold}}_4C/{{sample}}.interactions.homer",
            output=config["output"]),
        gff=expand("{gff_folder}/{{gff_file}}.gff", gff_folder=config["gff_folder"])
    output:
        expand("{output}/{{sample}}_merged/s7_{{resolution}}_{{threshold}}_4C/Virtual_4C/{{gff_file}}.wig",
            output=config["output"])
    params:
        resolution=lambda wildcards: wildcards.resolution,
    conda:
        "conda_environments/hicsd_env.yaml"
    shell:
        """
        hicsd \
         virtual_4C \
         --matrix_file {input.matrix} \
         --gff_file {input.gff} \
         --track_name $(echo $(basename {input.gff}) | sed "s/.gff//") \
         --bin_size {params.resolution} \
         --output_file {output}
        """

## ************************************************************************************
## step 7 - Ranking annotations
## ************************************************************************************
rule rank_annotations:
    input:
        matrix=expand("{output}/{{sample}}_merged/s7_{{resolution}}_{{threshold}}_4C/{{sample}}.ploidy_normalized.uniMulti",
            output=config["output"]),
        gff=expand("{gff_folder}/{{gff_file}}.gff", gff_folder=config["gff_folder"])
    output:
        expand("{output}/{{sample}}_merged/s8_{{resolution}}_{{threshold}}_4C/{{sample}}.{{gff_file}}.{{annotation}}.ranking.gff",
            output=config["output"])
    shell:
        "python bin/rank_annotations.py -v {input.gff} {input.matrix} {{resolution}} {{genomeAnnotation}} {output}"
