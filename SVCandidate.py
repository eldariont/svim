class CandidateInsertion:
    def __init__(self, sourceContig, sourceStart, sourceEnd, destContig, destStart, destEnd, evidence, read):
        self.sourceContig = sourceContig
        self.sourceStart = sourceStart
        self.sourceEnd = sourceEnd

        self.destContig = destContig
        self.destStart = destStart
        self.destEnd = destEnd

        self.evidence = evidence
        self.read = read


class CandidateDeletion:
    def __init__(self, sourceContig, sourceStart, sourceEnd, evidence, read):
        self.sourceContig = sourceContig
        self.sourceStart = sourceStart
        self.sourceEnd = sourceEnd
        
        self.evidence = evidence
        self.read = read


class CandidateInversion:
    def __init__(self, sourceContig, sourceStart, sourceEnd, evidence, read):
        self.sourceContig = sourceContig
        self.sourceStart = sourceStart
        self.sourceEnd = sourceEnd
        
        self.evidence = evidence
        self.read = read


class CandidateDuplicationTandem:
    def __init__(self, sourceContig, sourceStart, sourceEnd, copies, evidence, read):
        self.sourceContig = sourceContig
        self.sourceStart = sourceStart
        self.sourceEnd = sourceEnd
        
        self.copies = copies

        self.evidence = evidence
        self.read = read


class CandidateDuplicationInterspersed:
    def __init__(self, sourceContig, sourceStart, sourceEnd, destContig, destStart, destEnd, evidence, read):
        self.sourceContig = sourceContig
        self.sourceStart = sourceStart
        self.sourceEnd = sourceEnd

        self.destContig = destContig
        self.destStart = destStart
        self.destEnd = destEnd

        self.evidence = evidence
        self.read = read