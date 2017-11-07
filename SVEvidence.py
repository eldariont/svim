class SVEvidence:
    def __init__(self, contig1, start, end, type, evidence, read, contig2 = None):
        self.contig1 = contig1
        if contig2 == None:
            self.contig2 = contig1
        else:
            self.contig2 = contig2
        self.start = start
        self.end = end
        self.type = type
        self.evidence = evidence
        self.read = read

    def as_tuple(self):
        return (self.contig1, self.start, self.end, self.type, self.contig2, self.evidence, self.read)