#!/usr/bin/env python

from optparse import OptionParser
import subprocess, os
import dna

############################################################
# compostbin.py
#
# CBCB version of CompostBin for unsupervised clustering
# of metagenomic sequences.  Original CompostBin code is
# available at http://sites.google.com/site/souravc/compostbin
# and citable as "Sourav Chatterji, Ichitaro Yamazaki,
# Zhaojun Bai and Jonathan Eisen, CompostBin: A DNA
# composition-based algorithm for binning environmental
# shotgun reads , to appear in RECOMB 2008."
#
# Author: David Kelley
############################################################

cb_bin = "/fs/szasmg/dakelley/classes/metagenomics/software/Scimm/CBCBCompostBin"

# add matlab bin to matlabpath...FAILS IF JOB IS RUN IN BACKGROUND?
if not os.environ.has_key('MATLABPATH'):
    os.environ['MATLABPATH'] = cb_bin
elif os.environ['MATLABPATH'].find(cb_bin) == -1:
    os.environ['MATLABPATH'] = cb_bin + ':' + os.environ['MATLABPATH']

############################################################
# Read
############################################################
class Read:
    def __init__(self, header, seq, k):
        self.header = header
        self.kmers = dna.canonical_kmers(dna.count_kmers(k,seq.upper(),True))


############################################################
# main
############################################################
def main():
    parser = OptionParser()
    parser.add_option('-r', dest='readsf', help='Fasta file of reads')
    parser.add_option('-c', dest='num_clust', type='int', help='Number of clusters')
    parser.add_option('-k', dest='k', type='int', help='k-mers to count')
    parser.add_option('--constrained', dest='constraintsf', help='File of constrained reads') # constrained reads needs more testing
    (options, args) = parser.parse_args()
    
    # count sequences
    S = 0
    for line in open(options.readsf):
        S += (line[0] == '>')

    # get constraints
    if options.constraintsf:
        constraints = [line.rstrip() for line in open(options.constraintsf).readlines()]

    # count kmers
    count_kmers(options.readsf, options.k, 'kmers.dat')
    
    # initialize to single cluster, REMEMBER THIS IS 1-BASED
    clusters = [range(1,S+1)]

    # while we have too few clusters
    while(len(clusters) < options.num_clust):
        # find min cut
        mincut = -1
        for c in range(len(clusters)):
            # tell matlab script which kmers to consider
            iout = open('indexes.dat','w')
            print >> iout, '\n'.join([str(i) for i in clusters[c]])
            iout.close()

            # tell matlab script constraints
            if options.constraintsf:
                cout = open('constraints.dat','w')
                print >> cout, '\n'.join([str(constraints[i-1]) for i in clusters[c]])
                cout.close();
            elif os.path.isfile('constraints.dat'):
                print 'Renaming constraints.dat as constraints.tmp to avoid confusion'
                os.rename('constraints.dat', 'constraints.tmp')                

            # matlab partition
            matlab_attempts = 0
            while matlab_attempts != -1 and matlab_attempts < 3:
                p = subprocess.Popen('matlab -nodisplay -nosplash -nodesktop -wait -r partition', shell=True)
                sts = os.waitpid(p.pid, 0)[1]

                if os.path.isfile('partition.txt'):
                    matlab_attempts = -1
                else:
                    print >> sys.stderr, 'Matlab failed, trying again.'
                    matlab_attempts += 1
            if matlab_attempts != -1:
                exit(1)

            # compare ncut value
            ncut = float(open('ncut.txt').readline().rstrip())
            if mincut == -1 or ncut < mincut:
                minc = c
                mincut = ncut
                os.rename('partition.txt','min_part.txt')
                
        # create new clusters with the min cut
        # in the original slot and a new one
        tmplist = []
        clusters.append([])
        i = 0
        for line in open('min_part.txt'):
            if int(line.rstrip()):
                tmplist.append(clusters[minc][i])
            else:
                clusters[-1].append(clusters[minc][i])
            i += 1
        clusters[minc] = tmplist

    # output clusters
    out_clust(options.readsf, clusters)

    os.remove('min_part.txt')


############################################################
# count_kmers
#
# Count kmers in the reads file, and output vectors
# for clustering
############################################################
def count_kmers(readsf, k, output_file):
    # load reads and count kmers
    datf = open(output_file, 'w')

    header = ''
    for line in open(readsf):
        if line[0] == '>':
            # add last (if not first)
            if header:
                # count kmers
                r = Read(header, seq, k)
                print_kmers(r.kmers, datf)
                
            header = line[1:].rstrip()
            seq = ''
        else:
            seq += line.rstrip()

    # finish last
    r = Read(header, seq, k)
    print_kmers(r.kmers, datf)
        
    datf.close()


############################################################
# print_kmers
#
# Print the kmer count vector to file f
############################################################
def print_kmers(kmer_counts, f):
    kmers = sorted(kmer_counts.keys())
    kmers_sum = sum(kmer_counts.values())
    print >> f, '\t'.join([str(kmer_counts[kmer] / kmers_sum) for kmer in kmers])


############################################################
# out_clust
#
# Output sequence headers and cluster mappings to
# 'partition.txt'
############################################################
def out_clust(readsf, clusters):
    # map indexes to sequence headers
    seqs = []
    for line in open(readsf):
        if line[0] == '>':
            seqs.append(line[1:].rstrip())

    # output header,cluster tuples
    out = open('partition.txt','w')
    for c in range(len(clusters)):
        for i in clusters[c]:
            print >> out, '%d\t%s' % (c,seqs[i-1])
    out.close()
            

############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
