#!/usr/bin/env python

from optparse import OptionParser
import os, glob, util, subprocess, sys, pdb
import imm_cluster, dna

############################################################
# lb_init.py
#
# Sample a subset of reads and run LikelyBin and an initial
# run of imm_cluster
############################################################

bin_dir = os.path.abspath(os.path.dirname(sys.argv[0]))

############################################################
# main
############################################################
def main():
    parser = OptionParser()

    parser.add_option('-r', dest='readsf', help='Fasta file of reads')
    parser.add_option('-n', dest='numreads', type='int', help='Number of reads to sample for MCMC')
    parser.add_option('-k', dest='k', type='int', help='Number of clusters')
    parser.add_option('-o', dest='order', type='int', help='Order of Markov model')
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
        dna.fasta_rand_big(options.numreads, options.readsf, 'sample.fa')
    else:
        if os.path.isfile('sample.fa') or os.path.islink('sample.fa'):
            os.remove('sample.fa')
        os.symlink(options.readsf, 'sample.fa')

    # LikelyBin
    p = subprocess.Popen('%s/mcmc.pl sample.fa -num_sources %d -chain_order %d -num_threads %d &> lb.log' % (bin_dir,options.k, options.order, options.proc), shell=True)
    os.waitpid(p.pid, 0)

    if os.path.isfile('sample.fa.binning.allprobs') and os.path.getsize('sample.fa.binning.allprobs') > 0:

        # initialize clusters
        init_clusters(options.readsf, options.soft_assign)
        
        # check for k clusters
        new_k = drop_empty(options.k, options.soft_assign)
    
        # run seed_only
        p = subprocess.Popen('%s/imm_cluster.py -k %d -r %s -p %d -s --seed_only %s >> lb.log' % (bin_dir, new_k, options.readsf, options.proc, em), shell=True)
        os.waitpid(p.pid, 0)


############################################################
# init_clusters.py
#
# Convert LikelyBin output to an initial partitioning of
# reads for imm_cluster
############################################################
def init_clusters(readsf, soft_assign):
    # load_mates
    mates = {}
    #if matesf:   ... just in case I need this later ...
    #    for line in open(options.mates_file):
    #        (lr,rr) = line.split()
    #        mates[lr] = rr
    #        mates[rr] = lr
            
    read_likes = {}
    for line in open('sample.fa.binning.allprobs'):
        a = line.split('\t')
        r = a[0].strip()
        read_likes[r] = [float(x) for x in a[1:]]
        k = len(a[1:])

    # Note that I'm assuming here that the cluster priors
    # are incorporated into the printed likelihoods.  I can't
    # be sure but I'd rather be wrong and have used a uniform
    # prior than reestimate the priors myself and be wrong and
    # double count them.  However, this means I am double
    # counting the prior for mated reads.

    # assign to clusters
    hard_clusters = {}
    soft_clusters = {}
    for r in read_likes:
        if not hard_clusters.has_key(r):  # mate may have been done
            # if mated, combine likelihood with 
            if mates.has_key(r):
                m = mates[r]
                clust_likes = [read_likes[r][i] + read_likes[m][i] for i in range(k)]
            else:
                m = r   # it works
                clust_likes = read_likes[r]

            # hard assignment
            (like_max, clust) = util.max_i(clust_likes)
            hard_clusters[r] = clust
            hard_clusters[m] = clust

            # soft assignment
            if soft_assign:
                sum_score = clust_likes[0]
                for i in range(1,k):
                    sum_score = imm_cluster.log_add(sum_score, clust_likes[i])

                soft_clusters[r] = []
                for i in range(k):
                    prob = math.exp(clust_likes[i] - sum_score)
                    if r != m:
                        soft_clusters[m] = []
                    if prob > imm_cluster.soft_assign_t:
                        soft_clusters[r].append((i,prob))
                        if r != m:
                            soft_clusters[m].append((i,prob))

    chunk_size = 50
    chunk_i = 0
    while chunk_i*chunk_size < k:
        # open files
        init_files = {}
        build_files = {}
        for c in range(chunk_i*chunk_size, min(k, (chunk_i+1)*chunk_size)):
            init_files[c] = open('cluster-%d.fa' % c, 'w')
            if soft_assign:
                build_files[c] = open('cluster-%d.build.fa' % c, 'w')

        # read fasta to cluster-*.fa
        for line in open(readsf):
            if line[0] == '>':
                r = line[1:].strip()  # front spaces are removed by LikelyBin
                if hard_clusters.has_key(r):
                    hc = hard_clusters[r]
                    if init_files.has_key(hc):
                        init_files[hc].write(line)
                        if soft_assign:
                            for (sc,p) in soft_clusters[r]:
                                build_files[sc].write('>%f;%s' % (p,line[1:]))
                else:
                    hc = -1

            elif hc != -1:
                if init_files.has_key(hc):
                    init_files[hc].write(line)
                    if soft_assign:
                        for (sc,p) in soft_clusters[r]:
                            build_files[sc].write(line)

        # close files
        for c in init_files:
            init_files[c].close()
            if soft_assign:
                build_files[c].close()

        # increment
        chunk_i += 1


############################################################
# drop_empty
#
# If a cluster file is empty, just get rid of it, and
# reorder the files
############################################################
def drop_empty(k, soft_assign):
    # find all empty cluster files
    nonempty = []
    empty = []
    for i in range(k):
        if soft_assign:
            f = 'cluster-%d.build.fa' % i
        else:
            f = 'cluster-%d.fa' % i

        if os.path.getsize(f) > 0:
            nonempty.append(i)
        else:
            empty.append(i)

    # move largest index of nonempty to index of empty
    for i in range(len(empty)):
        if soft_assign:
            ef = 'cluster-%d.build.fa' % empty[i]
            nef = 'cluster-%d.build.fa' % nonempty[-1]
        else:
            ef = 'cluster-%d.fa' % empty[i]
            nef = 'cluster-%d.fa' % nonempty[-1]
    
        if empty[i] < nonempty[-1]:
            os.rename(nef,ef)
            nonempty[-1] = empty[i]
            nonempty.sort()
        else:
            os.remove(ef)

    return k - len(empty)


def ratio_sort(x, y):
    if x['ratio'] < y['ratio']:
        return -1
    elif x['ratio'] > y['ratio']:
        return 1
    else:
        return 0
    

if __name__ == '__main__':
    main()
    #pdb.runcall(init_clusters, 'testScimm.fasta', False)
