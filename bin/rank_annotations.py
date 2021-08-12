__author__ = "Markus Schmidt"
__version__ = "1.0.0"
__email__ = "Markus.Schmidt@lmu.de"

import argparse

def bisect_left(a, x, lo=0, hi=None, key=lambda x: x):
    """Return the index where to insert item x in list a, assuming a is sorted.
    The return value i is such that all e in a[:i] have e < x, and all e in
    a[i:] have e >= x.  So if x already appears in the list, a.insert(x) will
    insert just before the leftmost x already there.
    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
    """
    if lo < 0:
        raise ValueError('lo must be non-negative')
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        # Use __lt__ to match the logic in list.sort() and in heapq
        if key(a[mid]) < x: lo = mid+1
        else: hi = mid
    return lo

def as_tsv(*args):
    ret = ""
    for arg in args:
        ret += str(arg)
        ret += "\t"
    return ret[:-1]

class Annotation:
    def __init__(self, chrom, db_name, annotation_type, from_pos, to_pos, _1, strand, _2, extras, *opt):
        self.chrom = chrom
        self.type = annotation_type
        self.from_pos = int(from_pos)
        self.to_pos = int(to_pos)
        self.strand = strand
        self.extras = extras
        self.opt = (opt)
        self.val = 0
        self.bin_id = 0

    def key(self):
        return self.chrom

    def sort_pos(self):
        return self.from_pos

    def id(self):
        es = self.extras.split(";")
        for e in es:
            key, val = e.split("=")
            if key == "ID":
                return val
        return "idless_" + str(self.type)

    def __str__(self):
        return as_tsv(
                      self.chrom,
                      self.type,
                      self.from_pos,
                      self.to_pos,
                      self.strand,
                      "interactions=" + str(self.val) + ";bin_id=" + str(self.bin_id) + ";" + self.extras,
                      *self.opt
                    )

def overlaps(chr_1, from_1, to_1, chr_2, from_2, to_2):
    if chr_1 != chr_2:
        return False
    return to_1 > from_2 and to_2 > from_1

class Annotations:
    def __init__(self, annotation_file, annotation_filter, bin_size):
        self.data = {}
        self.filter = annotation_filter
        self.encountered_annos = set()

        with open(annotation_file, "r") as in_file_1:
            for line in in_file_1:
                if line[0] == "#":
                    continue
                # parse file colum
                anno = Annotation(*line.split())
                if anno.type in annotation_filter or "all" in annotation_filter:
                    self._add(anno)
                    self.encountered_annos.add(anno.type)
        self._prep(bin_size)

    def _add(self, anno):
        if not anno.key() in self.data:
            self.data[anno.key()] = []
        self.data[anno.key()].append(anno)

    def _prep(self, bin_size):
        bin_id_max = 0
        for val in self.data.values():
            bin_id_max_chr = bin_id_max
            val.sort(key=lambda x: x.sort_pos())

            # register overlapping annotations
            for anno in val:
                from_id = (anno.from_pos // bin_size) + bin_id_max
                to_id = (anno.to_pos // bin_size) + bin_id_max
                bin_id_max_chr = max(bin_id_max_chr, from_id, to_id)
                if from_id == to_id:
                    anno.bin_id = str(from_id)
                else:
                    anno.bin_id = str(from_id) + "-" + str(to_id) 
            bin_id_max = bin_id_max_chr + 1


    def interaction(self, chr_, from_pos, to_pos, val):
        if chr_ not in self.data:
            return
        
        idx = bisect_left(self.data[chr_], to_pos, key=lambda x: x.sort_pos()) - 1
        while idx >= 0 and self.data[chr_][idx].to_pos > from_pos:
            self.data[chr_][idx].val += val
            idx -= 1

    def iter_types(self):
        if "all" in self.filter:
            for anno in self.encountered_annos:
                yield anno
            yield "all"
        else:
            for anno in self.filter:
                yield anno

    def rank(self):
        for anno_type in self.iter_types():
            ret = []
            for key, val in self.data.items():
                for anno in val:
                    if anno.type == anno_type or anno_type is "all":
                        ret.append(anno)
            ret.sort(key=lambda x: (-x.val, x.chrom, x.from_pos))
            yield anno_type, ret

def parse_interactions(in_filenames):
    for in_filename in in_filenames:
        with open(in_filename, "r") as in_file_1:
            for line in in_file_1:
                # parse file columns
                chr_1, pos_1, chr_2, pos_2, val = line.split()

                yield in_filename, chr_1, int(pos_1), chr_2, int(pos_2), int(val)


def main(interaction_file, bin_size, annotation_file, annotation_filter, viewpoint_file, count, out_filename, args):
    annos = Annotations(annotation_file, annotation_filter, bin_size)

    if not viewpoint_file is None:
        with open(viewpoint_file, "r") as in_file_1:
            for line in in_file_1:
                if line[0] == "#":
                    continue
                view_chr, _, _, view_from, view_to, *_ = line.split()
                view_from = int(view_from)
                view_to = int(view_to)


    for file_name, chr_1, pos_1, chr_2, pos_2, val in parse_interactions(interaction_file):
        if not viewpoint_file is None:
            if not overlaps(chr_1, pos_1, pos_1 + bin_size, view_chr, view_from, view_to) and \
               not overlaps(chr_2, pos_2, pos_2 + bin_size, view_chr, view_from, view_to):
                continue

        annos.interaction(chr_1, pos_1 - bin_size//2, pos_1 + bin_size//2, val if not count else 1)
        annos.interaction(chr_2, pos_2 - bin_size//2, pos_2 + bin_size//2, val if not count else 1)

    with open(out_filename, "w") as out_file:
        out_file.write("##parameters:\n")
        for arg in vars(args):
            out_file.write("##" + str(arg) + "=" + str(getattr(args, arg)) + "\n")
        if not viewpoint_file is None:
            out_file.write("##viewpoint: " + view_chr + ":" + str(view_from) + "-" + str(view_to) + "\n")
        for anno_type, anno_list in annos.rank():
            out_file.write("##ranking for " + str(anno_type) + "\n")
            out_file.write("#")
            out_file.write(as_tsv("contig", "type", "start", "end",
                                  "strand", "bin_id", "extra"))
            out_file.write("\n")
            for anno in anno_list:
                out_file.write(str(anno))
                out_file.write("\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Rank annotations by Hi-C interactions.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('interactions', type=str,
                        help='File containing the Hi-C interactions. Use a comma seperated list for adding multiple files.')
    parser.add_argument('bin_size', type=int,
                        help='bin size used in the interactions file.')
    parser.add_argument('annotations', type=str,
                        help='File containing the genome annotations.')
    parser.add_argument('-f', "--filter", type=str, default="gene",
                        help='Comma seperated list of annotation types to rank. Parameter can be set to \'all\' for ranking all types of annotation')
    parser.add_argument('-v', "--viewpoint", type=str,
                        help='Genomic viewpoint for the interactions as gff file.')
    parser.add_argument('-c', "--count", action="store_true",
                        help='Count the number of bins that interactions occur with instead of summing up all interactions.')
    parser.add_argument('out_filename', type=str,
                        help='Prefix of the output file.')

    args = parser.parse_args()
    main(args.interactions.strip().split(","), args.bin_size, args.annotations, args.filter.strip().split(","),
         args.viewpoint, args.count, args.out_filename, args)