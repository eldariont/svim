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

    
    def get_bed_entry(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.contig, self.start, self.end, self.type + ";" + self.evidence, 0, self.read)


    def mean_distance_to(self, evidence2):
        """Return distance between means of two evidences."""
        if self.contig == evidence2.contig and self.type == evidence2.type:
            return abs(((self.start + self.end) / 2) - ((evidence2.start + evidence2.end) / 2))
        else:
            return float("inf")


class EvidenceDeletion(Evidence):
    def __init__(self, contig, start, end, evidence, read):
        self.contig = contig
        self.start = start
        self.end = end
        self.evidence = evidence
        self.read = read
        self.type = "del"


class EvidenceInsertion(Evidence):
    def __init__(self, contig, start, end, evidence, read):
        self.contig = contig
        self.start = start
        self.end = end
        self.evidence = evidence
        self.read = read
        self.type = "ins"
    

class EvidenceInversion(Evidence):
    def __init__(self, contig, start, end, evidence, read):
        self.contig = contig
        self.start = start
        self.end = end
        self.evidence = evidence
        self.read = read
        self.type = "inv"


class EvidenceTranslocation(Evidence):
    def __init__(self, contig1, pos1, contig2, pos2, evidence, read):
        self.contig1 = contig1
        self.pos1 = pos1
        self.contig2 = contig2
        self.pos2 = pos2
        self.evidence = evidence
        self.read = read
        self.type = "tra"

    def get_bed_entry(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.contig1, self.pos1, self.pos1+1, "{0};{1}:{2};{3}".format(self.type, self.contig2, self.pos2, self.evidence), 0, self.read)


class EvidenceDuplicationTandem(Evidence):
    def __init__(self, contig, start, end, copies, evidence, read):
        self.contig = contig
        self.start = start
        self.end = end
        
        self.copies = copies

        self.evidence = evidence
        self.read = read
        self.type = "dup"


    def get_bed_entry(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.contig, self.start, self.end, "{0};{1};{2}".format(self.type, self.copies,self.evidence), 0, self.read)


    def mean_distance_to(self, evidence2):
        """Return distance between means of two evidences."""
        if self.contig == evidence2.contig and self.type == evidence2.type and self.copies == evidence2.copies:
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
            return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.contig, self.start, self.end, "{0};{1};{2}".format(self.type, self.members[0].copies, self.size), self.score, ";".join([ev.read for ev in self.members]))
        else:
            return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.contig, self.start, self.end, "{0};{1}".format(self.type, self.size), self.score, ";".join([ev.read for ev in self.members]))
