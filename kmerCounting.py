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


def find_svs(fasta, winSize = 50, k = 7, debug = False):
    s1 = time()

    sequences = read_fasta(fasta)

    print "Reading finished ({0} s)".format(time()-s1)

    band = 0.5

    counting_time = 0
    start_time = time()

    num_pairs = len(sequences) / 2
    #for d in xrange(1, num_pairs+1):
    for d in xrange(1, 17):
        ref = sequences["ref" + str(d)]
        read = sequences["read" + str(d)]

        lengths = (len(ref), len(read))
        rows = len(ref) / winSize
        cols = len(read) / winSize
        #if len(ref) % winSize > 0:
            #rows += 1
        #if len(read) % winSize > 0:
            #cols += 1
        buckets = (rows, cols)

        s2 = time()
        counts = np.zeros(buckets, dtype=int)

        #Prepare kmer set for each read bucket
        ykmers = []
        for ybucket in xrange(buckets[1]):
            bucketkmers = set()
            for i in xrange(winSize):
                if (ybucket*winSize + i + k) <= len(read):
                    bucketkmers.add(read[(ybucket*winSize + i) : (ybucket*winSize + i + k)])
            ykmers.append(bucketkmers)
        
        for xbucket in xrange(buckets[0]):
            for i in xrange(winSize):
                for ybucket in xrange(buckets[1]):
                    if (xbucket*winSize + i + k) <= len(ref):
                        if ref[(xbucket*winSize + i) : (xbucket*winSize + i + k)] in ykmers[ybucket]:
                            counts[xbucket, ybucket] += 1

        #counts = c_count(ref, read, buckets, winSize, k)
        
        print "Counting finished for pair {0} ({1} s)".format(d, time()-s2)
        counting_time += time()-s2

        s3 = time()
        counts2 = np.nan_to_num(stats.zscore(counts, axis=1)  + stats.zscore(counts, axis=0))
        print "Z-Score computation finished for pair {0} ({1} s)".format(d, time()-s3)
        #print counts2

        if buckets[0] > buckets[1]:
            posLim = int(band * buckets[1])
            negLim = - (buckets[0] - buckets[1]) - int(band * buckets[1])
        else:
            posLim = int(band * buckets[0]) + (buckets[1] - buckets[0])
            negLim = -int(band * buckets[0])

        s4 = time()
        counts3 = np.zeros(buckets, dtype=int)
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
        print "Line finding finished for pair {0} ({1} s)".format(d, time()-s4)

        final_segments = convert_segments(completed_segments)
        if debug:
            print "Final segments", final_segments
        best_path = find_best_path(final_segments, buckets, debug)
        #print best_path
        for index in xrange(len(best_path) - 1):
            deletion_gap = best_path[index+1]['start'][0] - best_path[index]['end'][0]
            insertion_gap = best_path[index+1]['start'][1] - best_path[index]['end'][1]
            
            if insertion_gap > 1:
                print "Insertion found:", best_path[index]['end'][1], best_path[index+1]['start'][1]
            if deletion_gap > 1:
                print "Deletion found:", best_path[index]['end'][0], best_path[index+1]['start'][0]

        #np.putmask(counts3, counts2<5, 0)
        plt.figure(1)
        plt.subplot(4, 4, d)
        plt.imshow(counts3)
        plt.figure(2)
        plt.subplot(4, 4, d)
        plt.imshow(counts2)

    plt.show()
    total_time = time() - start_time
    print "Total time: {0}s; counting time: {1}s".format(total_time, counting_time)


def main():
    options = parseArguments(sys.argv)
    find_svs(options.fasta, debug = options.debug)
    
if __name__ == "__main__":
    sys.exit(main())
