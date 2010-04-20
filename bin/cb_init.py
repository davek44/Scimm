#!/usr/bin/env python

from optparse import OptionParser
import os, glob, subprocess
import scimm, dna

############################################################
# cb_init.py
#
# Sample a subset of reads and run compostbin and an 
# initial run of imm_cluster
############################################################

def main():
    parser = OptionParser()

    parser.add_option('-r', dest='readsf', help='Fasta file of reads')
    parser.add_option('-n', dest='numreads', type='int', help='Number of reads to sample for MCMC')
    parser.add_option('-k', dest='clusters', type='int', help='Number of clusters')
    parser.add_option('-m', dest='mers', type='int', help='Mers to count')
    parser.add_option('-p', dest='proc', type='int', help='Number of processes to run')
    parser.add_option('--em', dest='soft_assign', action='store_true', default=False, help='Use a soft assignment of reads to clusters')

    (options, args) = parser.parse_args()

    if options.soft_assign:
        em = '--em'
    else:
        em = ''

    # randomly sample reads
    total_reads = 0
    for line in open(options.readsf):
        if line[0] == '>':
            total_reads += 1
    if options.numreads and options.numreads < total_reads:
        dna.fasta_rand(options.numreads, options.readsf, 'sample.fa')
    else:
        os.system('ln -s %s sample.fa' % options.readsf)

    # CompostBin
    os.system('%s/compostbin.py -r sample.fa -c %d -k %d &> cb.log' % (scimm.scimm_bin, options.clusters, options.mers))

    # initialize clusters
    init_clusters(options.readsf, options.clusters, options.soft_assign)

    # run seed_only
    os.system('%s/imm_cluster.py -k %d -r %s -p %d -s --seed_only %s >> cb.log' % (scimm.scimm_bin, options.clusters, options.readsf, options.proc, em))
    

############################################################
# init_clusters
#
# Convert CompostBin output to an initial partitioning of
# reads for imm_cluster
############################################################
def init_clusters(readsf, clusters, soft_assign):
    read_clusters = {}
    for line in open('partition.txt'):
        (c,r) = line.split('\t')
        read_clusters[r.rstrip()] = int(c)

    # open files
    init_files = []
    build_files = []
    for c in range(clusters):
        init_files.append(open('cluster-%d.fa' % c, 'w'))
        if soft_assign:
            build_files.append(open('cluster-%d.build.fa' % c, 'w'))
        
    # read fasta to cluster-*.fa
    for line in open(readsf):
        if line[0] == '>':
            r = line[1:].rstrip()
            if read_clusters.has_key(r):
                c = read_clusters[r]
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
    for c in range(clusters):
        init_files[c].close()
        if soft_assign:
            build_files[c].close()

############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
