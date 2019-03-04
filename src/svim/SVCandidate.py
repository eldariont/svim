class Candidate:
    """Candidate class for structural variant candidates. Candidates reflect the final SV types and can be merged from signatures of several reads.
    """
    def __init__(self, source_contig, source_start, source_end, members, score, std_span, std_pos):
        self.source_contig = source_contig
        self.source_start = source_start
        self.source_end = source_end
        
        self.members = members
        self.score = score
        self.std_span = std_span
        self.std_pos = std_pos
        self.type = "unk"


    def get_source(self):
        return (self.source_contig, self.source_start, self.source_end)


    def get_key(self):
        contig, start, end = self.get_source()
        return (self.type, contig, (start + end) // 2)


    def mean_distance_to(self, candidate2):
        """Return distance between means of two candidates."""
        this_contig, this_start, this_end = self.get_source()
        other_contig, other_start, other_end = candidate2.get_source()
        if self.type == candidate2.type and this_contig == other_contig:
            return abs(((this_start +this_end) // 2) - ((other_start + other_end) // 2))
        else:
            return float("inf")


    def get_bed_entry(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(self.source_contig, self.source_start, self.source_end, "{0};{1};{2}".format(self.type, self.std_span, self.std_pos), self.score, ".", "["+"][".join([ev.as_string("|") for ev in self.members])+"]")


    def get_vcf_entry(self):
        raise NotImplementedError


class CandidateDeletion(Candidate):
    def __init__(self, source_contig, source_start, source_end, members, score, std_span, std_pos):
        self.source_contig = source_contig
        #0-based start of the deletion (first deleted base)
        self.source_start = source_start
        #0-based end of the deletion (one past the last deleted base)
        self.source_end = source_end

        self.members = members
        self.score = score
        self.std_span = std_span
        self.std_pos = std_pos
        self.type = "del"


    def get_vcf_entry(self):
        contig, start, end = self.get_source()
        svtype = "DEL"
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(contig, start, ".", "N", "<" + svtype + ">", int(self.score), "q30" if self.score < 30 else "PASS", "SVTYPE={0};END={1};SVLEN={2};STD_SPAN={3};STD_POS={4}".format(svtype, end, start - end, self.std_span, self.std_pos), "GT", "./.")


class CandidateInversion(Candidate):
    def __init__(self, source_contig, source_start, source_end, members, score, std_span, std_pos):
        self.source_contig = source_contig
        #0-based start of the inversion (first inverted base)
        self.source_start = source_start
        #0-based end of the inversion (one past the last inverted base)
        self.source_end = source_end

        self.members = members
        self.score = score
        self.std_span = std_span
        self.std_pos = std_pos
        self.type = "inv"


    def get_vcf_entry(self):
        contig, start, end = self.get_source()
        svtype = "INV"
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(contig, start+1, ".", "N", "<" + svtype + ">", int(self.score), "q20" if self.score < 20 else "PASS", "SVTYPE={0};END={1};STD_SPAN={2};STD_POS={3}".format(svtype, end, self.std_span, self.std_pos), "GT", "./.")


class CandidateNovelInsertion(Candidate):
    def __init__(self, dest_contig, dest_start, dest_end, members, score, std_span, std_pos):
        self.dest_contig = dest_contig
        #0-based start of the insertion (base after the insertion)
        self.dest_start = dest_start
        #0-based start of the insertion (base after the insertion) + length of the insertion
        self.dest_end = dest_end

        self.members = members
        self.score = score
        self.std_span = std_span
        self.std_pos = std_pos
        self.type = "nov_ins"

    def get_destination(self):
        return (self.dest_contig, self.dest_start, self.dest_end)

    def get_bed_entry(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(self.dest_contig, self.dest_start, self.dest_end, "{0};{1};{2}".format(self.type, self.std_span, self.std_pos), self.score, ".", "["+"][".join([ev.as_string("|") for ev in self.members])+"]")

    def get_vcf_entry(self):
        contig, start, end = self.get_destination()
        svtype = "INS:NOVEL"
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(contig, start, ".", "N", "<" + svtype + ">", int(self.score), "q30" if self.score < 30 else "PASS", "SVTYPE={0};END={1};SVLEN={2};STD_SPAN={3};STD_POS={4}".format(svtype, start, end - start, self.std_span, self.std_pos), "GT", "./.")


class CandidateDuplicationTandem(Candidate):
    def __init__(self, source_contig, source_start, source_end, copies, members, score, std_span, std_pos):
        self.source_contig = source_contig
        #0-based start of the region (first copied base)
        self.source_start = source_start
        #0-based end of the region (one past the last copied base)
        self.source_end = source_end
        
        self.copies = copies

        self.members = members
        self.score = score
        self.std_span = std_span
        self.std_pos = std_pos
        self.type = "dup_tan"


    def get_destination(self):
        source_contig, source_start, source_end = self.get_source()
        return (source_contig, source_end, source_end + self.copies * (source_end - source_start))


    def get_bed_entries(self, sep="\t"):
        source_contig, source_start, source_end = self.get_source()
        dest_contig, dest_start, dest_end = self.get_destination()
        source_entry = sep.join(["{0}", "{1}", "{2}", "{3}", "{4}", "{5}", "{6}"]).format(source_contig, source_start,
                                                                                     source_end,
                                                                                     "tan_dup_source;>{0}:{1}-{2};{3};{4}".format(
                                                                                         dest_contig, dest_start,
                                                                                         dest_end, self.std_span, self.std_pos), self.score, ".",
                                                                                     "[" + "][".join(
                                                                                         [ev.as_string("|") for ev in
                                                                                          self.members]) + "]")
        dest_entry = sep.join(["{0}", "{1}", "{2}", "{3}", "{4}", "{5}", "{6}"]).format(dest_contig, dest_start, dest_end,
                                                                                   "tan_dup_dest;<{0}:{1}-{2};{3};{4}".format(
                                                                                       source_contig, source_start,
                                                                                       source_end, self.std_span, self.std_pos), self.score, ".",
                                                                                   "[" + "][".join(
                                                                                       [ev.as_string("|") for ev in
                                                                                        self.members]) + "]")
        return (source_entry, dest_entry)


    def mean_distance_to(self, candidate2):
        """Return distance between means of two candidates."""
        this_source_contig, this_source_start, this_source_end = self.get_source()
        this_dest_contig, this_dest_start, this_dest_end = self.get_destination()
        other_source_contig, other_source_start, other_source_end = candidate2.get_source()
        other_dest_contig, other_dest_start, other_dest_end = candidate2.get_destination()
        if self.type == candidate2.type and this_source_contig == other_source_contig and this_dest_contig == other_dest_contig:
            return abs(((this_source_start + this_source_end) // 2) - ((other_source_start + other_source_end) // 2)) + \
                   abs(((this_dest_start + this_dest_end) // 2) - ((other_dest_start + other_dest_end) // 2))
        else:
            return float("inf")


    def get_vcf_entry(self):
        contig = self.source_contig
        start = self.source_end
        end = self.source_end + self.copies * (self.source_end - self.source_start)
        svtype = "DUP:TANDEM"
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(contig, start, ".", "N", "<" + svtype + ">", int(self.score), "q30" if self.score < 30 else "PASS", "SVTYPE={0};END={1};SVLEN={2};STD_SPAN={3};STD_POS={4}".format(svtype, start, end - start, self.std_span, self.std_pos), "GT", "./.")


class CandidateDuplicationInterspersed(Candidate):
    def __init__(self, source_contig, source_start, source_end, dest_contig, dest_start, dest_end, members, score, std_span, std_pos, cutpaste=False):
        self.source_contig = source_contig
        #0-based start of the region (first copied base)
        self.source_start = source_start
        #0-based end of the region (one past the last copied base)
        self.source_end = source_end

        self.dest_contig = dest_contig
        #0-based start of the insertion (base after the insertion)
        self.dest_start = dest_start
        #0-based end of the insertion (base after the insertion) + length of the insertion
        self.dest_end = dest_end

        self.members = members
        self.score = score
        self.std_span = std_span
        self.std_pos = std_pos
        self.cutpaste= cutpaste
        self.type = "dup_int"


    def get_destination(self):
        return (self.dest_contig, self.dest_start, self.dest_end)


    def get_bed_entries(self, sep="\t"):
        source_contig, source_start, source_end = self.get_source()
        dest_contig, dest_start, dest_end = self.get_destination()
        source_entry = sep.join(["{0}", "{1}", "{2}", "{3}", "{4}", "{5}", "{6}"]).format(source_contig, source_start,
                                                                                     source_end,
                                                                                     "int_dup_source;>{0}:{1}-{2};{3};{4}".format(
                                                                                         dest_contig, dest_start,
                                                                                         dest_end, self.std_span, self.std_pos), self.score, "origin potentially deleted" if self.cutpaste else ".",
                                                                                     "[" + "][".join(
                                                                                         [ev.as_string("|") for ev in
                                                                                          self.members]) + "]")
        dest_entry = sep.join(["{0}", "{1}", "{2}", "{3}", "{4}", "{5}", "{6}"]).format(dest_contig, dest_start, dest_end,
                                                                                   "int_dup_dest;<{0}:{1}-{2};{3};{4}".format(
                                                                                       source_contig, source_start,
                                                                                       source_end, self.std_span, self.std_pos), self.score, "origin potentially deleted" if self.cutpaste else ".",
                                                                                   "[" + "][".join(
                                                                                       [ev.as_string("|") for ev in
                                                                                        self.members]) + "]")
        return (source_entry, dest_entry)


    def mean_distance_to(self, candidate2):
        """Return distance between means of two candidates."""
        this_source_contig, this_source_start, this_source_end = self.get_source()
        this_dest_contig, this_dest_start, this_dest_end = self.get_destination()
        other_source_contig, other_source_start, other_source_end = candidate2.get_source()
        other_dest_contig, other_dest_start, other_dest_end = candidate2.get_destination()
        if self.type == candidate2.type and this_source_contig == other_source_contig and this_dest_contig == other_dest_contig:
            return abs(((this_source_start + this_source_end) // 2) - ((other_source_start + other_source_end) // 2)) + \
                   abs(((this_dest_start + this_dest_end) // 2) - ((other_dest_start + other_dest_end) // 2))
        else:
            return float("inf")


    def get_vcf_entry(self):
        contig, start, end = self.get_destination()
        svtype = "DUP:INT"
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(contig, start, ".", "N", "<" + svtype + ">", int(self.score), "q30" if self.score < 30 else "PASS", "SVTYPE={0};{1}END={2};SVLEN={3};STD_SPAN={4};STD_POS={5}".format(svtype, "CUTPASTE;" if self.cutpaste else "", start, end - start, self.std_span, self.std_pos), "GT", "./.")
