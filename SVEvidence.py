import logging

class Evidence:
    """Evidence class for basic evidences of structural variants. An evidence is always detected from a single read.
    """
    def __init__(self, contig, start, end, evidence, read):
        self.contig = contig
        self.start = start
        self.end = end
        self.evidence = evidence
        self.read = read
        self.type = "unk"
        if self.end < self.start:
            logging.warning("Evidence with invalid coordinates (end < start): " + self.as_string())


    def get_source(self):
        return (self.contig, self.start, self.end)


    def get_key(self):
        contig, start, end = self.get_source()
        return (self.type, contig, (start + end) // 2)


    def mean_distance_to(self, evidence2):
        """Return distance between means of two evidences."""
        this_contig, this_start, this_end = self.get_source()
        other_contig, other_start, other_end = evidence2.get_source()
        if self.type == evidence2.type and this_contig == other_contig:
            return abs(((this_start +this_end) // 2) - ((other_start + other_end) // 2))
        else:
            return float("inf")


    def gowda_diday_distance(self, evidence2, largest_indel_size):
        """Return Gowda-Diday distance between two evidences."""
        this_contig, this_start, this_end = self.get_source()
        other_contig, other_start, other_end = evidence2.get_source()
        # non-intersecting
        if this_contig != other_contig or this_end <= other_start or other_end <= this_start:
            return float("inf")
        dist_pos = abs(this_start - other_start) / largest_indel_size
        span1 = abs(this_end - this_start)
        span2 = abs(other_end - other_start)
        span_total = abs(max(this_end, other_end) - min(this_start, other_start))
        dist_span = abs(span1 - span2) / span_total
        inter = min(this_end, other_end) - max(this_start, other_start)
        dist_content = (span1 + span2 - 2 * inter) / span_total
        return dist_pos + dist_span + dist_content


    def span_loc_distance(self, evidence2, distance_normalizer):
        this_contig, this_start, this_end = self.get_source()
        other_contig, other_start, other_end = evidence2.get_source()
        if this_contig != other_contig:
            return float("inf")
        #Component 1: difference in spans
        this_span = this_end - this_start
        other_span = other_end - other_start
        dist_span = abs(this_span - other_span) / max(this_span, other_span)
        #Component 2: difference in locations
        this_center = (this_start + this_end) // 2
        other_center = (other_start + other_end) // 2
        dist_loc = min(abs(this_start - other_start), abs(this_end - other_end), abs(this_center - other_center)) / distance_normalizer
        return dist_span + dist_loc

    def as_string(self, sep="\t"):
        contig, start, end = self.get_source()
        return sep.join(["{0}","{1}","{2}","{3}","{4}"]).format(contig, start, end, "{0};{1}".format(self.type, self.evidence), self.read)


class EvidenceDeletion(Evidence):
    """SV Evidence: a region (contig:start-end) has been deleted and is not present in sample"""
    def __init__(self, contig, start, end, evidence, read):
        self.contig = contig
        self.start = start
        self.end = end
        self.evidence = evidence
        self.read = read
        self.type = "del"


class EvidenceInsertion(Evidence):
    """SV Evidence: a region of length end-start has been inserted at contig:start"""
    def __init__(self, contig, start, end, evidence, read):
        self.contig = contig
        self.start = start
        self.end = end
        self.evidence = evidence
        self.read = read
        self.type = "ins"


class EvidenceInversion(Evidence):
    """SV Evidence: a region (contig:start-end) has been inverted in the sample"""
    def __init__(self, contig, start, end, evidence, read, direction):
        self.contig = contig
        self.start = start
        self.end = end
        self.evidence = evidence
        self.read = read
        self.type = "inv"
        self.direction = direction


    def as_string(self, sep="\t"):
        contig, start, end = self.get_source()
        return sep.join(["{0}","{1}","{2}","{3}","{4}"]).format(contig, start, end, "{0};{1};{2}".format(self.type, self.direction, self.evidence), self.read)


class EvidenceInsertionFrom(Evidence):
    """SV Evidence: a region (contig:start-end) has been inserted at contig2:pos in the sample"""
    def __init__(self, contig1, start, end, contig2, pos, evidence, read):
        self.contig1 = contig1
        self.start = start
        self.end = end

        self.contig2 = contig2
        self.pos = pos

        self.evidence = evidence
        self.read = read
        self.type = "ins_dup"


    def get_source(self):
        return (self.contig1, self.start, self.end)


    def get_destination(self):
        source_contig, source_start, source_end = self.get_source()
        return (self.contig2, self.pos, self.pos + (source_end - source_start))


    def get_key(self):
        source_contig, source_start, source_end = self.get_source()
        dest_contig, dest_start, dest_end = self.get_destination()
        return (self.type, source_contig, dest_contig, dest_start + (source_start + source_end) // 2)


    def mean_distance_to(self, evidence2):
        """Return distance between means of two evidences."""
        this_source_contig, this_source_start, this_source_end = self.get_source()
        this_dest_contig, this_dest_start, this_dest_end = self.get_destination()
        other_source_contig, other_source_start, other_source_end = evidence2.get_source()
        other_dest_contig, other_dest_start, other_dest_end = evidence2.get_destination()
        if self.type == evidence2.type and this_source_contig == other_source_contig and this_dest_contig == other_dest_contig:
            return abs(((this_source_start + this_source_end) // 2) - ((other_source_start + other_source_end) // 2)) + \
                   abs(((this_dest_start + this_dest_end) // 2) - ((other_dest_start + other_dest_end) // 2))
        else:
            return float("inf")


    def as_string(self, sep="\t"):
        source_contig, source_start, source_end = self.get_source()
        dest_contig, dest_start, dest_end = self.get_destination()
        return sep.join(["{0}:{1}-{2}","{3}:{4}-{5}","{6}", "{7}"]).format(source_contig, source_start, source_end,
                                                                           dest_contig, dest_start, dest_end,
                                                                           "{0};{1}".format(self.type, self.evidence), self.read)


class EvidenceDuplicationTandem(Evidence):
    """SV Evidence: a region (contig:start-end) has been tandemly duplicated"""

    def __init__(self, contig, start, end, copies, evidence, read):
        self.contig = contig
        self.start = start
        self.end = end

        self.copies = copies

        self.evidence = evidence
        self.read = read
        self.type = "dup"


    def get_destination(self):
        source_contig, source_start, source_end = self.get_source()
        return (source_contig, source_end, source_end + self.copies * (source_end - source_start))


    def as_string(self, sep="\t"):
        source_contig, source_start, source_end = self.get_source()
        dest_contig, dest_start, dest_end = self.get_destination()
        return sep.join(["{0}:{1}-{2}","{3}:{4}-{5}","{6}", "{7}"]).format(source_contig, source_start, source_end,
                                                                           dest_contig, dest_start, dest_end,
                                                                           "{0};{1};{2}".format(self.type, self.evidence, self.copies), self.read)


class EvidenceTranslocation(Evidence):
    """SV Evidence: two positions (contig1:pos1 and contig2:pos2) are connected in the sample"""
    def __init__(self, contig1, pos1, direction1, contig2, pos2, direction2, evidence, read):
        self.contig1 = contig1
        self.pos1 = pos1
        self.direction1 = direction1
        self.contig2 = contig2
        self.pos2 = pos2
        self.direction2 = direction2
        self.evidence = evidence
        self.read = read
        self.type = "tra"


    def get_source(self):
        return (self.contig1, self.pos1, self.pos1 + 1)


    def get_destination(self):
        return (self.contig2, self.pos2, self.pos2 + 1)


    def as_string(self, sep="\t"):
        source_contig, source_start, source_end = self.get_source()
        dest_contig, dest_start, dest_end = self.get_destination()
        return sep.join(["{0}:{1}-{2}","{3}:{4}-{5}","{6}", "{7}"]).format(source_contig, source_start, source_end,
                                                                           dest_contig, dest_start, dest_end,
                                                                           "{0};{1}".format(self.type, self.evidence), self.read)


    def get_key(self):
        return (self.type, self.contig1, self.pos1)


    def mean_distance_to(self, evidence2):
        """Return distance between means of two evidences."""
        if self.type == evidence2.type and self.contig1 == evidence2.contig1:
            return abs(self.pos1 - evidence2.pos1)
        else:
            return float("inf")


class EvidenceClusterUniLocal(Evidence):
    def __init__(self, contig, start, end, score, size, members, type, std_span, std_pos):
        self.contig = contig
        self.start = start
        self.end = end
        
        self.score = score
        self.std_span = std_span
        self.std_pos = std_pos
        self.size = size
        self.members = members
        self.type = type


    def get_bed_entry(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.contig, self.start, self.end, "{0};{1};{2};{3}".format(self.type, self.size, self.std_span, self.std_pos), self.score, "["+"][".join([ev.as_string("|") for ev in self.members])+"]")


    def get_vcf_entry(self):
        if self.type == "del":
            svtype = "DEL"
        elif self.type == "ins":
            svtype = "INS"
        elif self.type == "inv":
            svtype = "INV"
        else:
            return
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(self.contig, self.start+1, ".", "N", "<" + svtype + ">", ".", "PASS", "SVTYPE={0};END={1};SVLEN={2};STD_SPAN={3};STD_POS={4}".format(svtype, self.end, self.end - self.start, self.std_span, self.std_pos))


    def get_length(self):
        return self.end - self.start

class EvidenceClusterBiLocal(Evidence):
    def __init__(self, source_contig, source_start, source_end, dest_contig, dest_start, dest_end, score, size, members, type, std_span, std_pos):
        self.source_contig = source_contig
        self.source_start = source_start
        self.source_end = source_end
        self.dest_contig = dest_contig
        self.dest_start = dest_start
        self.dest_end = dest_end
        
        self.score = score
        self.std_span = std_span
        self.std_pos = std_pos
        self.size = size
        self.members = members
        self.type = type


    def get_source(self):
        return (self.source_contig, self.source_start, self.source_end)


    def get_destination(self):
        return (self.dest_contig, self.dest_start, self.dest_end)


    def get_bed_entries(self):
        source_entry = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.source_contig, self.source_start, self.source_end, "{0}_source;{1}:{2}-{3};{4};{5};{6}".format(self.type, self.dest_contig, self.dest_start, self.dest_end, self.size, self.std_span, self.std_pos), self.score, "["+"][".join([ev.as_string("|") for ev in self.members])+"]")
        dest_entry = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.dest_contig, self.dest_start, self.dest_end, "{0}_dest;{1}:{2}-{3};{4}".format(self.type, self.source_contig, self.source_start, self.source_end, self.size), self.score, "["+"][".join([ev.as_string("|") for ev in self.members])+"]")
        return (source_entry, dest_entry)


    def get_vcf_entry(self):
        if self.type == "dup":
            svtype = "DUP:TANDEM"
        else:
            return
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(self.source_contig, self.source_start+1, ".", "N", "<" + svtype + ">", ".", "PASS", "SVTYPE={0};END={1};SVLEN={2};STD_SPAN={3};STD_POS={4}".format(svtype, self.source_end, self.source_end - self.source_start, self.std_span, self.std_pos))


    def get_source_length(self):
        return self.source_end - self.source_start


    def get_destination_length(self):
        return self.dest_end - self.dest_start