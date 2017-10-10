class callPacParams:
    def __init__(self):
        # Read tail mapping
        self.tail_span = 1000
        self.tail_min_mapq = 30
        self.tail_min_deviation = -0.1
        self.tail_max_deviation = 0.2

        # Counting
        self.count_win_size = 50
        self.count_k = 7
        self.count_band = 0.5

        # Find stretches
        self.stretch_threshold = 7
        self.stretch_tolerance = 2
        self.stretch_min_length = 3

        # Find best path through segments
        self.path_constant_gap_cost = 0
        self.path_convex_gap_cost = 3
        self.path_tolerance = 2

        # Alignment
        self.align_costs = (3, -12, -12)
