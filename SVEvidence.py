class Evidence:
    def __init__(self, contig, start, end, evidence, read):
        self.contig = contig
        self.start = start
        self.end = end
        self.evidence = evidence
        self.read = read
        self.type = "unk"


    def get_key(self):
        return (self.type, self.contig, (self.start + self.end) / 2)


    def get_bed_contig(self):
        return self.contig


    def get_bed_start(self):
        return self.start


    def get_bed_end(self):
        return self.end


    def get_bed_name(self):
        return "{0};{1}".format(self.type, self.evidence)


    def get_bed_entry(self, sep="\t"):
        return sep.join(["{0}","{1}","{2}","{3}","{4}","{5}",]).format(self.get_bed_contig(), self.get_bed_start(), self.get_bed_end(), self.get_bed_name(), 0, self.read)


    def mean_distance_to(self, evidence2):
        """Return distance between means of two evidences."""
        if self.type == evidence2.type and self.contig == evidence2.contig:
            return abs(((self.start + self.end) / 2) - ((evidence2.start + evidence2.end) / 2))
        else:
            return float("inf")


    def gowda_diday_distance(self, evidence2, largest_indel_size):
        """Return Gowda-Diday distance between two evidences."""
        # non-intersecting
        if self.end <= evidence2.start or evidence2.end <= self.start:
            return float("inf")
        dist_pos = abs(self.start - evidence2.start) / float(largest_indel_size)
        span1 = abs(self.end - self.start)
        span2 = abs(evidence2.end - evidence2.start)
        span_total = abs(max(self.end, evidence2.end) - min(self.start, evidence2.start))
        dist_span = abs(span1 - span2) / float(span_total)
        inter = min(self.end, evidence2.end) - max(self.start, evidence2.start)
        dist_content = (span1 + span2 - 2 * inter) / float(span_total)
        return dist_pos + dist_span + dist_content



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


    def get_key(self):
        return (self.type, self.contig1, self.contig2, self.pos + (self.start + self.end) / 2)


    def get_bed_contig(self):
        return self.contig2


    def get_bed_start(self):
        return self.pos


    def get_bed_end(self):
        return self.pos + (self.end - self.start)


    def get_bed_name(self):
        return "{0};{1}:{2}-{3};{4}".format(self.type, self.contig1, self.start, self.end, self.evidence)


    def mean_distance_to(self, evidence2):
        """Return distance between means of two evidences."""
        if self.type == evidence2.type and self.contig1 == evidence2.contig1 and self.contig2 == evidence2.contig2:
            return abs(((self.start + self.end) / 2) - ((evidence2.start + evidence2.end) / 2)) + abs(self.pos - evidence2.pos)
        else:
            return float("inf")


class EvidenceInversion(Evidence):
    """SV Evidence: a region (contig:start-end) has been inverted in the sample"""
    def __init__(self, contig, start, end, evidence, read):
        self.contig = contig
        self.start = start
        self.end = end
        self.evidence = evidence
        self.read = read
        self.type = "inv"


class EvidenceTranslocation(Evidence):
    """SV Evidence: two positions (contig1:pos1 and contig2:pos2) are connected in the sample"""
    def __init__(self, contig1, pos1, contig2, pos2, evidence, read):
        self.contig1 = contig1
        self.pos1 = pos1
        self.contig2 = contig2
        self.pos2 = pos2
        self.evidence = evidence
        self.read = read
        self.type = "tra"


    def get_bed_contig(self):
        return self.contig1


    def get_bed_start(self):
        return self.pos1


    def get_bed_end(self):
        return self.pos1 + 1


    def get_bed_name(self):
        return "{0};{1}:{2};{3}".format(self.type, self.contig2, self.pos2, self.evidence)


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


    def get_key(self):
        return (self.type, self.copies, self.contig, (self.start + self.end) / 2)


    def get_bed_name(self):
        return "{0};{1};{2}".format(self.type, self.copies,self.evidence)


    def mean_distance_to(self, evidence2):
        """Return distance between means of two evidences."""
        if self.type == evidence2.type and self.contig == evidence2.contig and self.copies == evidence2.copies:
            return abs(((self.start + self.end) / 2) - ((evidence2.start + evidence2.end) / 2))
        else:
            return float("inf")


class EvidenceCluster(Evidence):
    def __init__(self, contig, start, end, score, size, members, type):
        self.contig = contig
        self.start = start
        self.end = end
        
        self.score = score
        self.size = size
        self.members = members
        self.type = type
    
    def get_bed_entry(self):
        if self.type == "dup":
            return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.contig, self.start, self.end, "{0};{1};{2}".format(self.type, self.members[0].copies, self.size), self.score, "["+"][".join([ev.get_bed_entry("|") for ev in self.members])+"]")
        else:
            return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.contig, self.start, self.end, "{0};{1}".format(self.type, self.size), self.score, "["+"][".join([ev.get_bed_entry("|") for ev in self.members])+"]")
