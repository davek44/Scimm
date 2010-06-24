#!/usr/bin/env python
from optparse import OptionParser
import os, glob, subprocess, math, random
import scimm, util, dna

############################################################
# physcimm.py
#
# Sequence Clustering with Interpolated Markov Models using
# classification with Phymm to initialize.
#
# Author: David Kelley
############################################################

phymmdir = '/fs/szasmg/dakelley/classes/metagenomics/software/Scimm/phymm'

############################################################
# main
############################################################
def main():
    parser = OptionParser()

    # generic options
    parser.add_option('-s','-r', dest='readsf', help='Fasta file of sequences')
    parser.add_option('-p', dest='proc', type='int', default=2, help='Number of processes to run')
    parser.add_option('--em',dest='soft_assign', action='store_true', default=False, help='Use a soft assignment of reads to clusters')

    # phymm options
    parser.add_option('--taxlevel', dest='taxlevel', default='family', help='Taxonomic level at which to cluster reads with Phymm')
    parser.add_option('--minbp_pct', dest='minbp_pct', type='float', default=.01, help='Minimum proportion of bp assigned to a class to become a cluster')
    parser.add_option('-n','--numreads', dest='numreads', type='int', default=3000, help='Number of reads to sample from the data set to classify with Phymm')
    parser.add_option('-i', dest='ignore', help='Ask Phymm to ignore the IMMs in the given file')
    parser.add_option('--init', dest='init', action='store_true', default=False, help='Just initialize the clusters with Phymm; do not run cluster with IMMs')
    #parser.add_option('--nophymm', dest='nophymm', action='store_true', default=False, help='Phymm results have already been computed, and are in results.txt')

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

    # move to phymm directory
    origdir = os.getcwd()
    os.chdir(phymmdir)

    # choose a random id for the phymm dir
    pid_unchosen = True
    while pid_unchosen:
        pid = random.randint(0,1000000)
        if not os.path.isfile('sample%d.fa'%pid):
            pid_unchosen = False

    # randomly sample reads
    total_reads = 0
    for line in open(options.readsf):
        if line[0] == '>':
            total_reads += 1
    if options.numreads and options.numreads < total_reads:
        dna.fasta_rand_big(options.numreads, options.readsf, 'sample%d.fa' % pid)
    else:
        os.system('ln -s %s sample%d.fa' % (options.readsf,pid))
        options.numreads = total_reads

    # classify
    phymm_parallel(pid, options.proc, origdir, options.ignore)

    # determine minimum bp for cluster
    total_bp = 0
    for line in open('sample%d.fa'%pid):
        if line[0] != '>':
            total_bp += len(line.rstrip())
    minbp = options.minbp_pct*total_bp

    # clean up phymm
    os.system('rm *sample%d*' % pid)

    # move back to regular directory
    os.chdir(origdir)    

    # initialize clusters
    class_k = init_clusters(options.readsf, options.taxlevel, minbp, options.soft_assign)    

    # run IMM clustering
    if not options.init:
        os.system('%s/imm_cluster.py -k %d -r %s -p %d -s %s &> immc.log' % (scimm.scimm_bin,class_k, options.readsf, options.proc, em))


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
                print 'Sorry, Phymm only considers fasta headers up to the first whitespace.  Please make these unique in your file.  E.f. %s is not unique' % r
                exit()
            reads[r] = True


############################################################
# phymm_parallel
#
# Break up file of reads and run Phymm in parallel. Output
# is placed in results.txt
############################################################
def phymm_parallel(pid, proc, origdir, ignoref):
    # open tmp files
    tmpfiles = []
    for i in range(proc):
        tmpfiles.append(open('sample%d.fa.%d' % (pid,i), 'w'))

    # distribute sequences among temp files
    seq_count = -1
    for line in open('sample%d.fa' % pid):
        if line[0] == '>':
            seq_count += 1
        tmpfiles[seq_count % proc].write(line)

    # close temp files
    for i in range(proc):
        tmpfiles[i].close()

    # launch Phymm
    cmds = []
    for i in range(proc):
        if ignoref:
            cmds.append('./scoreReads.pl sample%d.fa.%d -i %s 2> /dev/null' % (pid,i,ignoref))
        else:
            cmds.append('./scoreReads.pl sample%d.fa.%d 2> /dev/null' % (pid,i))
    util.exec_par(cmds, proc)

    # collect output
    open('%s/results.txt' % origdir,'w')
    for i in range(proc):
        os.system('cat results.01.phymm_sample%d_fa_%d.txt >> %s/results.txt' % (pid,i,origdir))


############################################################
# init_clusters.py
#
# Convert Phymm output to an initial partitioning of
# reads for imm_cluster
############################################################
def init_clusters(readsf, taxlevel, minbp, soft_assign):
    class2index = {'species':1, 'genus':3, 'family':4, 'order':5, 'class':6, 'phylum':7}
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
    for line in open('results.txt'):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        if a[0] != 'QUERY_ID' and a[col]: # some species are missing classifications
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
