class callPacParams:
    def __init__(self):
        # Read tail mapping
        self.tail_span = 1000
        self.tail_min_mapq = 30
        self.tail_min_deviation = -0.02
        self.tail_max_deviation = 0.1

        # Counting
        self.count_win_size = 100
        self.count_k = 13
        self.count_band = 0.5

        # Find stretches
        self.stretch_threshold = 4
        self.stretch_tolerance = 2
        self.stretch_min_length = 3

        # Find best path through segments
        self.path_constant_gap_cost = 0
        self.path_linear_gap_cost = 0.01
        self.path_convex_gap_cost = 0
        self.path_root_gap_cost = 1
        self.path_tolerance = 2

        # Alignment
        self.align_costs = (3, -12, -12)

    def set_with_options(self, options):
        # Read tail mapping
        self.tail_span = options.tail_span
        self.tail_min_mapq =  options.tail_min_mapq
        self.tail_min_deviation = options.tail_min_deviation
        self.tail_max_deviation = options.tail_max_deviation

        # Counting
        self.count_win_size = options.count_win_size
        self.count_k = options.count_k
        self.count_band = options.count_band

        # Find stretches
        self.stretch_threshold = options.stretch_threshold
        self.stretch_tolerance = options.stretch_tolerance
        self.stretch_min_length = options.stretch_min_length

        # Find best path through segments
        self.path_constant_gap_cost = options.path_constant_gap_cost
        self.path_linear_gap_cost = options.path_linear_gap_cost
        self.path_convex_gap_cost = options.path_convex_gap_cost
        self.path_root_gap_cost = options.path_root_gap_cost
        self.path_tolerance = options.path_tolerance

        # Alignment
        self.align_costs = (options.align_costs_match, options.align_costs_mismatch, options.align_costs_gap)
