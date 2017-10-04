import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from scipy import stats
from time import time
import sys
import math
import argparse
from cCounting import c_count

def parseArguments(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""""")
    parser.add_argument('fasta', type=argparse.FileType('r'))
    parser.add_argument('--debug', '-d', action='store_true')
    return parser.parse_args()


def read_fasta(filehandle):
    sequences = dict()
    seq = ""
    for line in filehandle:
        sline = line.strip()
        if sline.startswith('>'):
            if len(seq) > 0:
                sequences[name] = seq
                seq = ""
            name = sline[1:]
        else:
            seq += sline
    if len(seq) > 0:
        sequences[name] = seq
        seq = ""
    return sequences


def find_stretches(values, threshold = 7, tolerance = 2, min_length = 5):
    stretches = []
    begin = -1
    last_good = -1
    for index, value in enumerate(values):
        if value > threshold:
            if begin == -1:
                begin = index
            last_good = index
        else:
            if begin != -1 and index - last_good > tolerance:
                if last_good - begin >= min_length - 1:
                    stretches.append((begin, last_good))
                begin = -1
                last_good = -1
    if begin != -1 and last_good - begin >= min_length - 1:
        stretches.append((begin, last_good))
    return stretches


def convert_segments(segments):
    converted_segments = []
    for segment in segments:
        stretches = []
        for offset, start, end in segment:
            if offset > 0:
                stretches.append(((start, start + offset), (end, end + offset)))
            else:
                stretches.append(((start - offset, start), (end - offset, end)))
        start = sorted(stretches, key=lambda stretch: stretch[0][0] + stretch[0][1])[0][0]
        end = sorted(stretches, key=lambda stretch: stretch[1][0] + stretch[1][1], reverse=True)[0][1]
        new_segment = {'start': start, 'end': end, 'stretches': stretches}
        converted_segments.append(new_segment)
    return converted_segments


def distance(start1, end1, start2, end2, convexA = 0, convexB = 3,shift_tolerance = 2):
    if start2[0] + shift_tolerance >= end1[0] and start2[1] + shift_tolerance >= end1[1]:
        if start2[0] - end1[0] <= 0:
            dist = 0
        else:
            dist = convexA + (convexB * math.log(start2[0] - end1[0]))
        if start2[1] - end1[1] <= 0:
            dist += 0
        else:
            dist += convexA + (convexB * math.log(start2[1] - end1[1]))
        return dist
    else:
        return float("inf")


def find_best_path(segments, matrix_end, debug = False):
    #print segments
    G = nx.DiGraph()
    G.add_nodes_from(["start", "end"])
    G.add_edge("start", "end", weight = distance((0,0), (0,0), matrix_end, matrix_end))
    for i1, s1 in enumerate(segments):
        for i2, s2 in enumerate(segments):
            if i1 != i2:
                dist = distance(s1['start'], s1['end'], s2['start'], s2['end'])
                if dist < float('inf'):
                    G.add_edge(i1, i2, weight = dist)
        dist = distance((0,0), (0,0), s1['start'], s1['end'])
        if dist < float('inf'):
            G.add_edge("start", i1, weight = dist)
        dist = distance(s1['start'], s1['end'], matrix_end, matrix_end)
        if dist < float('inf'):
            G.add_edge(i1, "end", weight = dist)
    if debug:
        print "Edge-list:", nx.to_edgelist(G)
    shortest_path = nx.shortest_path(G, "start", "end", weight = "weight")
    start_segment = {'start': (0,0), 'end': (0,0), 'stretches': []}
    end_segment = {'start': matrix_end, 'end': matrix_end, 'stretches': []}
    return [start_segment] + [segments[node] for node in shortest_path[1:-1]] + [end_segment]


def plot_array(array, rows, cols, d, figure = 1):
    plt.figure(figure)
    plt.subplot(rows, cols, d)
    plt.imshow(array)


def find_svs(ref, read, winSize = 50, k = 7, debug = False, times = False):
    band = 0.5
    start_time = time()

    lengths = (len(ref), len(read))
    rows = len(ref) / winSize
    cols = len(read) / winSize
    last_row_size = len(ref) % winSize
    last_col_size = len(read) % winSize
    if last_row_size > (winSize / 3):
        rows += 1
    if last_col_size > (winSize / 3):
        cols += 1

    s2 = time()
    counts = np.zeros((rows, cols), dtype=int)

    #Prepare kmer set for each read bucket
    ykmers = []
    for ybucket in xrange(cols):
        bucketkmers = set()
        for i in xrange(winSize):
            if (ybucket*winSize + i + k) <= len(read):
                bucketkmers.add(read[(ybucket*winSize + i) : (ybucket*winSize + i + k)])
        ykmers.append(bucketkmers)
    
    for xbucket in xrange(rows):
        for i in xrange(winSize):
            for ybucket in xrange(cols):
                if (xbucket*winSize + i + k) <= len(ref):
                    if ref[(xbucket*winSize + i) : (xbucket*winSize + i + k)] in ykmers[ybucket]:
                        counts[xbucket, ybucket] += 1
    
    if last_row_size > (winSize / 3):
        for col in xrange(len(read) / winSize):
            counts[rows - 1, col] = counts[rows - 1, col] * (winSize / last_row_size)
    if last_col_size > (winSize / 3):
        for row in xrange(len(ref) / winSize):
            counts[row, cols - 1] = counts[row, cols - 1] * (winSize / last_col_size)
    if last_row_size > (winSize / 3) and last_col_size > (winSize / 3):
        counts[rows - 1, cols - 1] = counts[rows - 1, cols - 1] * (winSize * winSize) / (last_row_size * last_col_size)
    
    if times:
        print "The size of the last row/column was {0}/{1}bps.".format(last_row_size, last_col_size)
        print "Counting finished ({0} s)".format(time()-s2)
    
    #s2 = time()
    #counts2 = c_count(ref, read, (rows, cols), winSize, k)
    #print "Counting finished ({0} s)".format(time()-s2)

    #if (counts == counts2).all():
        #print "Same"
    #else:
        #print counts - counts2

    s3 = time()
    counts2 = np.nan_to_num(stats.zscore(counts, axis=1)  + stats.zscore(counts, axis=0))
    if times:
        print "Z-Score computation finished ({0} s)".format( time()-s3)
    #print counts2

    if rows > cols:
        posLim = int(band * cols)
        negLim = - (rows - cols) - int(band * cols)
    else:
        posLim = int(band * rows) + (cols - rows)
        negLim = -int(band * rows)

    s4 = time()
    counts3 = np.zeros((rows, cols), dtype=int)
    completed_segments = []
    active_segments = []
    for offset in xrange(negLim, posLim):
        #Find stretches
        values = np.diagonal(counts2, offset)
        stretches = find_stretches(values)
        for start, end in stretches:
            for i in xrange(start, end+1):
                if offset >= 0:
                    counts3[i, i+offset] = 1
                else:
                    counts3[i-offset, i] = 1

        #Combine stretches to segment
        new_active_segments = []
        new_completed_segments = completed_segments + active_segments[:]
        for start, end in stretches:
            found_matching_segment = False
            for segment in active_segments:
                last_offset, last_start, last_end = segment[-1]
                if last_offset >= 0 and offset >= 0:
                    if start - 1 <= last_end and end >= last_start:
                        current_segment = segment[:]
                        current_segment.append((offset, start, end))
                        new_active_segments.append(current_segment)
                        if segment in new_completed_segments:
                            new_completed_segments.remove(segment)
                        found_matching_segment = True
                elif last_offset <= 0 and offset <= 0:
                    if start <= last_end and end + 1 >= last_start:
                        current_segment = segment[:]
                        current_segment.append((offset, start, end))
                        new_active_segments.append(current_segment)
                        if segment in new_completed_segments:
                            new_completed_segments.remove(segment)
                        found_matching_segment = True
            if not found_matching_segment:
                new_active_segments.append([(offset, start, end)])
        active_segments = new_active_segments
        completed_segments = new_completed_segments

    completed_segments.extend(active_segments)
    if times:
        print "Line finding finished ({0} s)".format(time()-s4)

    final_segments = convert_segments(completed_segments)
    if debug:
        print "Final segments", final_segments
    best_path = find_best_path(final_segments, (rows, cols), debug)
    #print best_path
    
    sv_results = []
    for index in xrange(len(best_path) - 1):
        deletion_gap = best_path[index+1]['start'][0] - best_path[index]['end'][0]
        insertion_gap = best_path[index+1]['start'][1] - best_path[index]['end'][1]
        
        if insertion_gap > 1:
            #print "Insertion found:", best_path[index]['end'][1], best_path[index+1]['start'][1]
            insertion_length = best_path[index+1]['start'][1] - best_path[index]['end'][1]
            sv_results.append( ('ins', best_path[index]['end'][0], best_path[index]['end'][0] + insertion_length) )
        if deletion_gap > 1:
            #print "Deletion found:", best_path[index]['end'][0], best_path[index+1]['start'][0]
            sv_results.append( ('del', best_path[index]['end'][0], best_path[index+1]['start'][0]) )

    #np.putmask(counts3, counts2<5, 0)

    total_time = time() - start_time
    if times:
        print "Total time: {0}s".format(total_time)
    if debug:
        return sv_results, counts, counts2, counts3
    
    return sv_results


def main():
    options = parseArguments(sys.argv)

    #Read fasta file
    s1 = time()
    sequences = read_fasta(options.fasta)
    #print "Reading finished ({0} s)".format(time()-s1)

    num_pairs = len(sequences) / 2
    for i, d in enumerate(xrange(16, 24, 1)):
        print "Read", d
        for key in sequences.keys():
            if key.startswith("ref" + str(d)):
                ref = sequences[key]
            if key.startswith("read" + str(d)):
                read = sequences[key]
        sv_results, counts, zscores, stretches = find_svs(ref, read, debug=True)
        plot_array(counts, 4, 2, i+1, 1)
        plot_array(zscores, 4, 2, i+1, 2)
        plot_array(stretches, 4, 2, i+1, 3)
        print sv_results
        print("")
    plt.show()
        
    
if __name__ == "__main__":
    sys.exit(main())
