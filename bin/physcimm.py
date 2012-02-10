#!/usr/bin/env python
from optparse import OptionParser, SUPPRESS_HELP
import os, glob, subprocess, math, random, sys
import scimm, util, dna

############################################################
# physcimm.py
#
# Sequence Clustering with Interpolated Markov Models using
# classification with Phymm to initialize.
#
# Author: David Kelley
############################################################

bin_dir = os.path.abspath(os.path.dirname(sys.argv[0]))


############################################################
# main
############################################################
def main():
    parser = OptionParser()

    # generic options
    parser.add_option('-s', dest='readsf', help='Fasta file of sequences')
    parser.add_option('-p', dest='proc', type='int', default=2, help='Number of processes to run')

    # phymm options
    parser.add_option('--taxlevel', dest='taxlevel', default='family', help='Taxonomic level at which to cluster reads with Phymm [Default=%default]')
    parser.add_option('--minbp_pct', dest='minbp_pct', type='float', default=.01, help='Minimum proportion of bp assigned to a class to become a cluster [Default=%default]')
    parser.add_option('-n','--numreads', dest='numreads', type='int', default=3000, help='Number of reads to sample from the data set to classify with Phymm [Default=%default]')
    parser.add_option('-r','--phymm_results', dest='phymm_results_file', help='Phymm results file to be used rather than running Phymm from scratch.')

    # my testing options
    # help='Use a soft assignment of reads to clusters [Default=%default]'
    parser.add_option('--em',dest='soft_assign', action='store_true', default=False, help=SUPPRESS_HELP)
    # help='Ask Phymm to ignore the IMMs in the given file'
    parser.add_option('-i', dest='ignore', help=SUPPRESS_HELP)
    # help='Run Phymm and initialize clusters only'
    parser.add_option('--init', dest='init', action='store_true', default=False, help=SUPPRESS_HELP)
    # help='Run my version of Phymm w/o Blast and w/ chromosomes only'
    parser.add_option('--bc', dest='bc', action='store_true', default=False, help=SUPPRESS_HELP)

    (options, args) = parser.parse_args()

    # check data
    data_integrity(options.readsf)
    
    # make robust to directory changes
    options.readsf = os.path.abspath(options.readsf)
    if options.ignore:
        options.ignore = os.path.abspath(options.ignore)    

    if options.soft_assign:
        em = '--em'
    else:
        em = ''

    if options.phymm_results_file:
        if not os.path.isfile('sample.fa') and not os.path.islink('sample.fa'):
            print >> sys.stderr, 'Assuming Phymm results file includes all reads'
            os.symlink(options.readsf, 'sample.fa')
        phymm_results_file = options.phymm_results_file

    else:
        # randomly sample reads
        total_reads = 0
        for line in open(options.readsf):
            if line[0] == '>':
                total_reads += 1
        if options.numreads and options.numreads < total_reads:
            dna.fasta_rand_big(options.numreads, options.readsf, 'sample.fa')
        else:
            os.symlink(options.readsf, 'sample.fa')

        # classify
        if options.bc:
            bc_str = '-b -c'
            phymm_results_file = 'results.01.phymm_sample_fa.txt'
        else:
            bc_str = ''
            phymm_results_file = 'results.03.phymmBL_sample_fa.txt'
        p = subprocess.Popen('%s/phymm_par.py -p %d %s sample.fa' % (bin_dir, options.proc,bc_str), shell=True)
        os.waitpid(p.pid, 0)

    # determine minimum bp for cluster
    total_bp = 0
    for line in open('sample.fa'):
        if line[0] != '>':
            total_bp += len(line.rstrip())
    minbp = options.minbp_pct*total_bp

    # initialize clusters
    class_k = init_clusters(options.readsf, phymm_results_file, options.taxlevel, minbp, options.soft_assign)

    # run IMM clustering
    if not options.init:
        p = subprocess.Popen('%s/imm_cluster.py -k %d -r %s -p %d -s %s &> immc.log' % (bin_dir, class_k, options.readsf, options.proc, em), shell=True)
        os.waitpid(p.pid, 0)
        


############################################################
# data_integrity
#
# Check for uniqueness of headers.
############################################################
def data_integrity(readsf):
    reads = {}
    for line in open(readsf):
        if line[0] == '>':
            r = line[1:].split()[0]
            if reads.has_key(r):
                print 'Sorry, Phymm only considers fasta headers up to the first whitespace.  Please make these unique in your file.  E.g. %s is not unique' % r
                exit()
            reads[r] = True


############################################################
# init_clusters.py
#
# Convert Phymm output to an initial partitioning of
# reads for imm_cluster
############################################################
def init_clusters(readsf, phymm_results_file, taxlevel, minbp, soft_assign):
    if open(phymm_results_file).readline().find('CONF') == -1:
        class2index = {'strain':1, 'species':1, 'genus':3, 'family':4, 'order':5, 'class':6, 'phylum':7}
    else:
        class2index = {'strain':1, 'species':1, 'genus':3, 'family':5, 'order':7, 'class':9, 'phylum':11}
    col = class2index[taxlevel.lower()]

    # get read sizes and map between Phymm headers and true headers
    readbp = {}
    phymm2true = {}
    for line in open(readsf):
        if line[0] == '>':
            header = line[1:].rstrip()
            readbp[header] = 0
            phymm2true[header.split()[0]] = header
        else:
            readbp[header] += len(line.rstrip())    

    # fill clusters with reads
    clusters = {}
    clustbp = {}
    for line in open(phymm_results_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        if a[0] != 'QUERY_ID' and a[col]: # some species are missing classifications
            if phymm2true.has_key(a[0]): # results file is allowed to have extra sequences
                header = phymm2true[a[0]]

                if clusters.has_key(a[col]):
                    clusters[a[col]].append(header)
                    clustbp[a[col]] += readbp[header]
                else:
                    clusters[a[col]] = [header]
                    clustbp[a[col]] = readbp[header]

    # extra cluster for deleted classes
    clusters['extra'] = []
    clustbp['extra'] = 0

    # number clusters, forget clusters w/ < minreads
    clust_nums = {}
    for c in clusters.keys():
        # don't consider deleting extra
        if c == 'extra':
            continue

        if clustbp[c] < minbp:
            clusters['extra'] += clusters[c]
            clustbp['extra'] += clustbp[c]
            del clusters[c]
        else:
            clust_nums[c] = len(clust_nums)

    # check extra size (for now, allowing anything)
    if clustbp['extra']:
        clust_nums['extra'] = len(clust_nums)
    else:
        del clusters['extra']

    # map reads -> clusters
    read2cluster = {}
    for c in clusters:
        for r in clusters[c]:
            read2cluster[r] = clust_nums[c]

    # open files
    init_files = []
    build_files = []
    for i in range(len(clusters)):
        init_files.append(open('cluster-%d.fa' % i, 'w'))
        if soft_assign:
            build_files.append(open('cluster-%d.build.fa' % i, 'w'))
        
    # read fasta to cluster-*.fa
    for line in open(readsf):
        if line[0] == '>':
            r = line[1:].rstrip()
            if read2cluster.has_key(r):
                c = read2cluster[r]
                init_files[c].write(line)
                if soft_assign:
                    build_files[c].write('>1.0;%s' % line[1:])
            else:
                c = -1

        elif c != - 1:
            init_files[c].write(line)
            if soft_assign:
                build_files[c].write(line)

    # close files
    for i in range(len(clusters)):
        init_files[i].close()
        if soft_assign:
            build_files[i].close()
    
    return len(clusters)

    
############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
