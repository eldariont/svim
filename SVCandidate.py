class Candidate:
    def __init__(self, source_contig, source_start, source_end, members, score):
        self.source_contig = source_contig
        self.source_start = source_start
        self.source_end = source_end
        
        self.members = members
        self.score = score


    def get_source(self):
        return (self.source_contig, self.source_start, self.source_end)


class CandidateDeletion(Candidate):
    pass


class CandidateInversion(Candidate):
    pass


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


class CandidateDuplicationTandem(Candidate):
    def __init__(self, source_contig, source_start, source_end, copies, members, score):
        self.source_contig = source_contig
        self.source_start = source_start
        self.source_end = source_end
        
        self.copies = copies

        self.members = members
        self.score = score


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