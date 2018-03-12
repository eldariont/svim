class Candidate:
    """Candidate class for structural variant candidates. Candidates reflect the final SV types and can be merged from evidences of several reads.
    """
    def __init__(self, source_contig, source_start, source_end, members, score):
        self.source_contig = source_contig
        self.source_start = source_start
        self.source_end = source_end
        
        self.members = members
        self.score = score
        self.type = "unk"


    def get_source(self):
        return (self.source_contig, self.source_start, self.source_end)


    def get_key(self):
        contig, start, end = self.get_source()
        return (self.type, contig, (start + end) / 2)


    def mean_distance_to(self, candidate2):
        """Return distance between means of two candidates."""
        this_contig, this_start, this_end = self.get_source()
        other_contig, other_start, other_end = candidate2.get_source()
        if self.type == candidate2.type and this_contig == other_contig:
            return abs(((this_start +this_end) / 2) - ((other_start + other_end) / 2))
        else:
            return float("inf")


    def gowda_diday_distance(self, candidate2, largest_indel_size):
        """Return Gowda-Diday distance between two evidences."""
        this_contig, this_start, this_end = self.get_source()
        other_contig, other_start, other_end = candidate2.get_source()
        # non-intersecting
        if this_end <= other_start or other_end <= this_start:
            return float("inf")
        dist_pos = abs(this_start - other_start) / float(largest_indel_size)
        span1 = abs(this_end - this_start)
        span2 = abs(other_end - other_start)
        span_total = abs(max(this_end, other_end) - min(this_start, other_start))
        dist_span = abs(span1 - span2) / float(span_total)
        inter = min(this_end, other_end) - max(this_start, other_start)
        dist_content = (span1 + span2 - 2 * inter) / float(span_total)
        return dist_pos + dist_span + dist_content


    def get_bed_entry(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.source_contig, self.source_start, self.source_end, "{0}".format(self.type), self.score, "["+"][".join([ev.as_string("|") for ev in self.members])+"]")


    def get_vcf_entry(self):
        raise NotImplementedError


class CandidateDeletion(Candidate):
    def get_vcf_entry(self):
        contig, start, end = self.get_source()
        svtype = "DEL"
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(contig, start+1, ".", "N", "<" + svtype + ">", ".", "PASS", "SVTYPE={0};END={1};SVLEN={2}".format(svtype, end, end - start))


class CandidateInversion(Candidate):
    def __init__(self, source_contig, source_start, source_end, members, score):
        self.source_contig = source_contig
        self.source_start = source_start
        self.source_end = source_end

        self.members = members
        self.score = score
        self.type = "inv"


    def get_vcf_entry(self):
        contig, start, end = self.get_source()
        svtype = "INV"
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(contig, start+1, ".", "N", "<" + svtype + ">", ".", "PASS", "SVTYPE={0};END={1};SVLEN={2}".format(svtype, end, end - start))


class CandidateInsertion(Candidate):
    def __init__(self, source_contig, source_start, source_end, dest_contig, dest_start, dest_end, members, score):
        self.source_contig = source_contig
        self.source_start = source_start
        self.source_end = source_end

        self.dest_contig = dest_contig
        self.dest_start = dest_start
        self.dest_end = dest_end

        self.members = members
        self.score = score
        self.type = "ins"


    def get_destination(self):
        return (self.dest_contig, self.dest_start, self.dest_end)

    
    def get_bed_entries(self, sep="\t"):
        source_contig, source_start, source_end = self.get_source()
        dest_contig, dest_start, dest_end = self.get_destination()
        source_entry = sep.join(["{0}","{1}","{2}","{3}","{4}","{5}",]).format(source_contig, source_start, source_end,
                                                                       "ins_source;>{0}:{1}-{2}".format(dest_contig, dest_start, dest_end), self.score, "["+"][".join([ev.as_string("|") for ev in self.members])+"]")
        dest_entry = sep.join(["{0}","{1}","{2}","{3}","{4}","{5}",]).format(dest_contig, dest_start, dest_end,
                                                                             "ins_dest;<{0}:{1}-{2}".format(source_contig, source_start, source_end), self.score, "["+"][".join([ev.as_string("|") for ev in self.members])+"]")
        return (source_entry, dest_entry)


    def mean_distance_to(self, candidate2):
        """Return distance between means of two candidates."""
        this_source_contig, this_source_start, this_source_end = self.get_source()
        this_dest_contig, this_dest_start, this_dest_end = self.get_destination()
        other_source_contig, other_source_start, other_source_end = candidate2.get_source()
        other_dest_contig, other_dest_start, other_dest_end = candidate2.get_destination()
        if self.type == candidate2.type and this_source_contig == other_source_contig and this_dest_contig == other_dest_contig:
            return abs(((this_source_start + this_source_end) / 2) - ((other_source_start + other_source_end) / 2)) + \
                   abs(((this_dest_start + this_dest_end) / 2) - ((other_dest_start + other_dest_end) / 2))
        else:
            return float("inf")


    def get_vcf_entry(self):
        contig, start, end = self.get_destination()
        svtype = "INS"
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(contig, start+1, ".", "N", "<" + svtype + ">", ".", "PASS", "SVTYPE={0};END={1};SVLEN={2}".format(svtype, end, end - start))


class CandidateDuplicationTandem(Candidate):
    def __init__(self, source_contig, source_start, source_end, copies, members, score):
        self.source_contig = source_contig
        self.source_start = source_start
        self.source_end = source_end
        
        self.copies = copies

        self.members = members
        self.score = score
        self.type = "dup_tan"


    def get_destination(self):
        source_contig, source_start, source_end = self.get_source()
        return (source_contig, source_end, source_end + self.copies * (source_end - source_start))


    def get_bed_entries(self, sep="\t"):
        source_contig, source_start, source_end = self.get_source()
        dest_contig, dest_start, dest_end = self.get_destination()
        source_entry = sep.join(["{0}", "{1}", "{2}", "{3}", "{4}", "{5}", ]).format(source_contig, source_start,
                                                                                     source_end,
                                                                                     "int_dup_source;>{0}:{1}-{2}".format(
                                                                                         dest_contig, dest_start,
                                                                                         dest_end), self.score,
                                                                                     "[" + "][".join(
                                                                                         [ev.as_string("|") for ev in
                                                                                          self.members]) + "]")
        dest_entry = sep.join(["{0}", "{1}", "{2}", "{3}", "{4}", "{5}", ]).format(dest_contig, dest_start, dest_end,
                                                                                   "int_dup_dest;<{0}:{1}-{2}".format(
                                                                                       source_contig, source_start,
                                                                                       source_end), self.score,
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
            return abs(((this_source_start + this_source_end) / 2) - ((other_source_start + other_source_end) / 2)) + \
                   abs(((this_dest_start + this_dest_end) / 2) - ((other_dest_start + other_dest_end) / 2))
        else:
            return float("inf")


    def get_vcf_entry(self):
        contig = self.source_contig
        start = self.source_end
        end = self.source_end + self.copies * (self.source_end - self.source_start)
        svtype = "DUP:TANDEM"
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(contig, start+1, ".", "N", "<" + svtype + ">", ".", "PASS", "SVTYPE={0};END={1};SVLEN={2}".format(svtype, end, end - start))


class CandidateDuplicationInterspersed(Candidate):
    def __init__(self, source_contig, source_start, source_end, dest_contig, dest_start, dest_end, members, score):
        self.source_contig = source_contig
        self.source_start = source_start
        self.source_end = source_end

        self.dest_contig = dest_contig
        self.dest_start = dest_start
        self.dest_end = dest_end

        self.members = members
        self.score = score
        self.type = "dup_int"


    def get_destination(self):
        return (self.dest_contig, self.dest_start, self.dest_end)


    def get_bed_entries(self, sep="\t"):
        source_contig, source_start, source_end = self.get_source()
        dest_contig, dest_start, dest_end = self.get_destination()
        source_entry = sep.join(["{0}", "{1}", "{2}", "{3}", "{4}", "{5}", ]).format(source_contig, source_start,
                                                                                     source_end,
                                                                                     "int_dup_source;>{0}:{1}-{2}".format(
                                                                                         dest_contig, dest_start,
                                                                                         dest_end), self.score,
                                                                                     "[" + "][".join(
                                                                                         [ev.as_string("|") for ev in
                                                                                          self.members]) + "]")
        dest_entry = sep.join(["{0}", "{1}", "{2}", "{3}", "{4}", "{5}", ]).format(dest_contig, dest_start, dest_end,
                                                                                   "int_dup_dest;<{0}:{1}-{2}".format(
                                                                                       source_contig, source_start,
                                                                                       source_end), self.score,
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
            return abs(((this_source_start + this_source_end) / 2) - ((other_source_start + other_source_end) / 2)) + \
                   abs(((this_dest_start + this_dest_end) / 2) - ((other_dest_start + other_dest_end) / 2))
        else:
            return float("inf")


    def get_vcf_entry(self):
        contig, start, end = self.get_destination()
        svtype = "DUP:INT"
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(contig, start+1, ".", "N", "<" + svtype + ">", ".", "PASS", "SVTYPE={0};END={1};SVLEN={2}".format(svtype, end, end - start))
