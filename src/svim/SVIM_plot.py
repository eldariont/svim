import matplotlib
import logging
import random
import math
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_sv_lengths(deletion_candidates, inversion_candidates, int_duplication_candidates, tan_dup_candidates, novel_insertion_candidates, options):
    len_dict_5 = dict()
    len_dict_10 = dict()
    len_dict_5["DEL"] = [v.get_source()[2] - v.get_source()[1] for v in deletion_candidates if v.score >= 5]
    len_dict_5["INV"] = [v.get_source()[2] - v.get_source()[1] for v in inversion_candidates if v.score >= 5]
    len_dict_5["DUP_INT"] = [v.get_destination()[2] - v.get_destination()[1] for v in int_duplication_candidates if v.score >= 5]
    len_dict_5["DUP_TAN"] = [v.get_destination()[2] - v.get_destination()[1] for v in tan_dup_candidates if v.score >= 5]
    len_dict_5["INS"] = [v.get_destination()[2] - v.get_destination()[1] for v in novel_insertion_candidates if v.score >= 5]
    draw_sv_length_plot(dict_of_lengths=len_dict_5, output=options.working_dir + "/sv-lengths-q5.png")
    len_dict_10["DEL"] = [v.get_source()[2] - v.get_source()[1] for v in deletion_candidates if v.score >= 10]
    len_dict_10["INV"] = [v.get_source()[2] - v.get_source()[1] for v in inversion_candidates if v.score >= 10]
    len_dict_10["DUP_INT"] = [v.get_destination()[2] - v.get_destination()[1] for v in int_duplication_candidates if v.score >= 10]
    len_dict_10["DUP_TAN"] = [v.get_destination()[2] - v.get_destination()[1] for v in tan_dup_candidates if v.score >= 10]
    len_dict_10["INS"] = [v.get_destination()[2] - v.get_destination()[1] for v in novel_insertion_candidates if v.score >= 10]
    draw_sv_length_plot(dict_of_lengths=len_dict_10, output=options.working_dir + "/sv-lengths-q10.png")


def draw_sv_length_plot(dict_of_lengths, output):
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


def plot_sv_alleles(candidates, options):
    refs_11 = [candidate.ref_reads for candidate in candidates if candidate.genotype == '1/1' and candidate.score >= 5 and candidate.ref_reads != None and candidate.alt_reads != None]
    alts_11 = [candidate.alt_reads for candidate in candidates if candidate.genotype == '1/1' and candidate.score >= 5 and candidate.ref_reads != None and candidate.alt_reads != None]
    refs_10 = [candidate.ref_reads for candidate in candidates if candidate.genotype == '0/1' and candidate.score >= 5 and candidate.ref_reads != None and candidate.alt_reads != None]
    alts_10 = [candidate.alt_reads for candidate in candidates if candidate.genotype == '0/1' and candidate.score >= 5 and candidate.ref_reads != None and candidate.alt_reads != None]
    refs_00 = [candidate.ref_reads for candidate in candidates if candidate.genotype == '0/0' and candidate.score >= 5 and candidate.ref_reads != None and candidate.alt_reads != None]
    alts_00 = [candidate.alt_reads for candidate in candidates if candidate.genotype == '0/0' and candidate.score >= 5 and candidate.ref_reads != None and candidate.alt_reads != None]
    refs_nn = [candidate.ref_reads for candidate in candidates if candidate.genotype == './.' and candidate.score >= 5 and candidate.ref_reads != None and candidate.alt_reads != None]
    alts_nn = [candidate.alt_reads for candidate in candidates if candidate.genotype == './.' and candidate.score >= 5 and candidate.ref_reads != None and candidate.alt_reads != None]

    draw_allele_plot(refs_11, alts_11, refs_10, alts_10, refs_00, alts_00, refs_nn, alts_nn, output=options.working_dir + "/sv-genotypes-q5.png")


def draw_allele_plot(refs_11, alts_11, refs_10, alts_10, refs_00, alts_00, refs_nn, alts_nn, output):
    """Makes a scatter plot of allele support
    """
    num_points = len(refs_11)+len(refs_10)+len(refs_00)+len(refs_nn)
    point_alpha = 10 / math.sqrt(max(100, num_points))
    plt.scatter(x=[ref+random.uniform(-0.5, 0.5) for ref in refs_11],
                y=[alt+random.uniform(-0.5, 0.5) for alt in alts_11],
                c='tab:red',
                alpha=point_alpha,
                label='1/1',
                edgecolors='none')
    plt.scatter(x=[ref+random.uniform(-0.5, 0.5) for ref in refs_10],
                y=[alt+random.uniform(-0.5, 0.5) for alt in alts_10],
                c='tab:purple',
                alpha=point_alpha,
                label='0/1',
                edgecolors='none')
    plt.scatter(x=[ref+random.uniform(-0.5, 0.5) for ref in refs_00],
                y=[alt+random.uniform(-0.5, 0.5) for alt in alts_00],
                c='tab:blue',
                alpha=point_alpha,
                label='0/0',
                edgecolors='none')
    plt.scatter(x=[ref+random.uniform(-0.5, 0.5) for ref in refs_nn],
                y=[alt+random.uniform(-0.5, 0.5) for alt in alts_nn],
                c='tab:brown',
                alpha=point_alpha,
                label='./.',
                edgecolors='none')
    axes = plt.gca()
    axes.set_xlim([0,60])
    axes.set_ylim([0,60])
    plt.xlabel('Reference allele support')
    plt.ylabel('Variant allele support')
    leg = plt.legend(frameon=True,
               fontsize="medium")
    for lh in leg.legendHandles:
        lh.set_alpha(1.0)

    plt.tight_layout()
    plt.savefig(output)
    plt.clf()