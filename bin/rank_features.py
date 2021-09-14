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
        self.num_bins = 1

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

    def get_val(self):
        return self.val / self.num_bins

    def __str__(self):
        return as_tsv(
                      self.chrom,
                      "rank_features.py",
                      self.type,
                      self.from_pos,
                      self.to_pos,
                      self.get_val(),
                      self.strand,
                      ".",
                      "interactions=" + str(self.get_val()) + ";bin_id=" + str(self.bin_id) + ";" + \
                            self.extras,
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
        self.preamble = ""

        with open(annotation_file, "r") as in_file_1:
            for line in in_file_1:
                if line[0] == "#":
                    self.preamble += line
                    continue
                # parse file colum
                anno = Annotation(*line[:-1].split("\t"))
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
                anno.num_bins = 1 + to_id - from_id
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
            ret.sort(key=lambda x: (-x.get_val(), x.chrom, x.from_pos))
            yield anno_type, ret

def parse_interactions(in_filenames):
    for in_filename in in_filenames:
        with open(in_filename, "r") as in_file_1:
            curr_chr = None
            bin_size = 0
            for line in in_file_1:
                if len(line) <= 1:
                    continue
                split_str = line[:-1].split(" ")
                front = split_str[0]
                if front == "track":
                    continue
                elif front == "variableStep":
                    curr_chr = split_str[1].split("=")[1]
                    bin_size = int(split_str[2].split("=")[1])
                else:
                    yield in_filename, curr_chr, int(front) - bin_size//2, int(front) + bin_size//2, float(split_str[1])


def main(interaction_file, bin_size, annotation_file, annotation_filter, count, out_filename, args):
    annos = Annotations(annotation_file, annotation_filter, bin_size)

    for file_name, chr_1, start, end, val in parse_interactions(interaction_file):
        annos.interaction(chr_1, start, end, val if not count else 1)

    with open(out_filename, "w") as out_file:
        out_file.write(annos.preamble)
        out_file.write("##rank_features parameters:\n")
        for arg in vars(args):
            out_file.write("##rank_features " + str(arg) + "=" + str(getattr(args, arg)) + "\n")
        for anno_type, anno_list in annos.rank():
            out_file.write("##rank_features feature=" + str(anno_type) + "\n")
            out_file.write("#")
            out_file.write(as_tsv("contig", "source", "feature", "start", "end", "score",
                                  "strand", "phase", "attributes"))
            out_file.write("\n")
            for anno in anno_list:
                out_file.write(str(anno))
                out_file.write("\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Rank features by Hi-C interactions.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('interactions', type=str,
                        help='File containing the Hi-C interactions. Use a comma seperated list for adding multiple files.')
    parser.add_argument('bin_size', type=int,
                        help='bin size used in the interactions file.')
    parser.add_argument('annotations', type=str,
                        help='File containing the genome annotations.')
    parser.add_argument('-f', "--filter", type=str, default="gene",
                        help='Comma seperated list of features to rank. Parameter can be set to \'all\' for ranking all types of features')
    parser.add_argument('-c', "--count", action="store_true",
                        help='Count the number of bins that interactions occur with instead of summing up all interactions.')
    parser.add_argument('out_filename', type=str,
                        help='Prefix of the output file.')

    args = parser.parse_args()
    main(args.interactions.strip().split(","), args.bin_size, args.annotations, args.filter.strip().split(","),
         args.count, args.out_filename, args)