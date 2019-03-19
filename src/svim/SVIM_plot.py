import matplotlib
import logging
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_sv_lengths(deletion_candidates, inversion_candidates, int_duplication_candidates, tan_dup_candidates, novel_insertion_candidates, options):
    len_dict_5 = dict()
    len_dict_10 = dict()
    logging.info("Generate SV length plots..")
    len_dict_5["DEL"] = [v.get_source()[2] - v.get_source()[1] for v in deletion_candidates if v.score >= 5]
    len_dict_5["INV"] = [v.get_source()[2] - v.get_source()[1] for v in inversion_candidates if v.score >= 5]
    len_dict_5["DUP_INT"] = [v.get_destination()[2] - v.get_destination()[1] for v in int_duplication_candidates if v.score >= 5]
    len_dict_5["DUP_TAN"] = [v.get_destination()[2] - v.get_destination()[1] for v in tan_dup_candidates if v.score >= 5]
    len_dict_5["INS"] = [v.get_destination()[2] - v.get_destination()[1] for v in novel_insertion_candidates if v.score >= 5]
    make_plot(dict_of_lengths=len_dict_5, output=options.working_dir + "/sv-lengths-q5.png")
    len_dict_10["DEL"] = [v.get_source()[2] - v.get_source()[1] for v in deletion_candidates if v.score >= 10]
    len_dict_10["INV"] = [v.get_source()[2] - v.get_source()[1] for v in inversion_candidates if v.score >= 10]
    len_dict_10["DUP_INT"] = [v.get_destination()[2] - v.get_destination()[1] for v in int_duplication_candidates if v.score >= 10]
    len_dict_10["DUP_TAN"] = [v.get_destination()[2] - v.get_destination()[1] for v in tan_dup_candidates if v.score >= 10]
    len_dict_10["INS"] = [v.get_destination()[2] - v.get_destination()[1] for v in novel_insertion_candidates if v.score >= 10]
    make_plot(dict_of_lengths=len_dict_10, output=options.working_dir + "/sv-lengths-q10.png")


def make_plot(dict_of_lengths, output):
    """Makes two stacked bar charts
    Plotting two bar charts of number of SVs by length split by SV type
    Use a consistent colouring scheme for those in "standard_order" to
    make comparison reasonable

    First bar chart is up to 2kb with bins of 10bp
    Second bar chart is up to 20kb, with bins of 100bp
     and uses log scaling on the y-axis
    """
    standard_order = ['DEL', 'INS', 'INV', 'DUP_INT', 'DUP_TAN']
    names, lengths = zip(
            *sorted([(svtype, lengths) for svtype, lengths in dict_of_lengths.items()],
                    key=lambda x: standard_order.index(x[0])))
    plt.subplot(2, 1, 1)
    plt.hist(x=lengths,
             bins=[i for i in range(0, 2000, 10)],
             stacked=True,
             histtype='bar',
             label=names)
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")

    plt.subplot(2, 1, 2)
    plt.hist(x=lengths,
             bins=[i for i in range(0, 20000, 100)],
             stacked=True,
             histtype='bar',
             label=names,
             log=True)
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")
    plt.tight_layout()
    plt.savefig(output)
    plt.clf()