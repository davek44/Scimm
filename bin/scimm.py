#!/usr/bin/env python

from optparse import OptionParser
import os, glob, subprocess, math
import imm_cluster, util

############################################################
# scimm.py
#
# Sequence Clustering with Interpolated Markov Models
#
# Author: David Kelley
############################################################

scimm_bin = "/fs/szasmg/dakelley/classes/metagenomics/software/Scimm/bin"
if 'PYTHONPATH' in os.environ:
    os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + ':' + scimm_bin
else:
    os.environ['PYTHONPATH'] = scimm_bin

############################################################
# main
############################################################
def main():
    parser = OptionParser()

    # generic options
    parser.add_option('-s','-r', dest='readsf', help='Fasta file of sequences')
    parser.add_option('-k', dest='k', type='int', help='Number of clusters')
    parser.add_option('-p', dest='proc', type='int', default=2, help='Number of processes to run [Default=%default]')
    parser.add_option('--em',dest='soft_assign', action='store_true', default=False, help='Use a soft assignment of reads to clusters [Default=%default]')

    # likelybin options
    parser.add_option('--ls', dest='lb_starts', type='int', default=1, help='Number of random LikelyBin starts [Default=%default]')
    parser.add_option('--ln', dest='lb_numreads', type='int', default=3000, help='Number of reads to sample for LikelyBin [Default=%default]')
    parser.add_option('--lt', dest='lb_threads', type='int', default=2, help='Number of LikelyBin threads per start, and CPUs for imm_cluster [Default=%default]')
    parser.add_option('--lo', dest='lb_order', type='int', default=3, help='Order of LikelyBin Markov model [Default=%default]')

    # compostbin options
    parser.add_option('--cs', dest='cb_starts', type='int', default=1, help='Number of random CompostBin starts [Default=%default]')
    parser.add_option('--cn', dest='cb_numreads', type='int', default=3000, help='Number of reads to sample for CompostBin [Default=%default]')
    parser.add_option('--ct', dest='cb_threads', type='int', default=1, help='Number of CPUs for imm_cluster [Default=%default]')
    parser.add_option('--co','--cm', dest='cb_mers', type='int', default=4, help='mers to count in CompostBin [Default=%default]')

    (options, args) = parser.parse_args()

    options.readsf = os.path.abspath(options.readsf)

    total_starts = options.lb_starts + options.cb_starts

    if options.soft_assign:
        em = '--em'
    else:
        em = ''

    # run initial samples
    i = 0    
    while i < total_starts:
        p = []
        j = 0
        while j < options.proc:
            # LikelyBin
            if i < options.lb_starts:
                # double check processes
                if j + options.lb_threads <= options.proc:
                    # make a temp dir to compute in and cd to it
                    temp_dir('tmp.start%d' % i)                    
                    p.append(subprocess.Popen('%s/lb_init.py -r %s -n %d -k %d -o %d -p %d %s' % (scimm_bin, options.readsf, options.lb_numreads, options.k, options.lb_order, options.lb_threads, em), shell=True))
                    os.chdir('..')
                    i += 1
                elif j == 0:
                    print 'Cannot use more lb threads than processes'
                    exit()
            
                j += options.lb_threads  # even if not true, just move things along

            # CompostBin
            else:
                # double check processes
                if j + options.cb_threads <= options.proc:
                    # make a temp dir to compute in and cd to it
                    temp_dir('tmp.start%d' % i)
                    p.append(subprocess.Popen('%s/cb_init.py -r %s -n %d -k %d -m %d -p %d %s' % (scimm_bin, options.readsf, options.cb_numreads, options.k, options.cb_mers, options.cb_threads, em), shell=True))
                    os.chdir('..')
                    i += 1
                elif j == 0:
                    print 'Cannot use more cb threads than processes'
                    exit()
            
                j += options.lb_threads  # even if not true, just move things along

        # wait for processes to finish
        for j in range(len(p)):
            os.waitpid(p[j].pid, 0)

    # choose best start
    #maxlike_clusters(total_starts, options.readsf, options.k, options.soft_assign)    
    minentropy_clusters(total_starts, options.readsf, options.k, options.soft_assign)    

    # delete others
    #for i in range(options.starts):
    #    os.system('rm -r tmp.start%d' % i)

    # in case k changed
    new_k = determine_k(options.soft_assign, options.k)

    # run imm clustering completely
    os.system('%s/imm_cluster.py -k %d -r %s -p %d -i --trained %s &> immc.log' % (scimm_bin, new_k, options.readsf, options.proc, em))


############################################################
# temp_dir
#
# Create and change to a temporary directory to do initial
# runs within
############################################################
def temp_dir(tmpdir):
    if os.path.isdir(tmpdir):
        os.chdir(tmpdir)
        if len(glob.glob('*')) > 0:
            os.system('rm *')
    else:
        os.mkdir(tmpdir)
        os.chdir(tmpdir)


############################################################
# maxlike_clusters
#
# Copy the clustering with maximum likelihood to the main
# directory 
############################################################
def maxlike_clusters(total_starts, readsf, k, soft_assign):
    like = [0]*total_starts
    for i in range(total_starts):
        os.chdir('tmp.start%d' % i)
        if len(glob.glob('cluster-*.fa')) > 0:
            # determine likelihood
            like[i] = scimm_like(readsf, k, soft_assign)
        else:
            # something failed
            like[i] = ''
        os.chdir('..')

    # find max likelihood initial partitioning
    max_like = min(like)   # '' is greater than numbers
    for i in range(len(like)):        
        if like[i] != '' and like[i] >= max_like:
            max_like = like[i]
            max_clust = i

    # get files from max
    os.system('cp tmp.start%d/cluster-*.fa tmp.start%d/icm-*scores.tmp .' % (max_clust,max_clust))


############################################################
# scimm_like
#
# Calculate the likelihood of the given clustering and IMM
############################################################
def scimm_like(readsf, k, soft_assign):
    new_k = determine_k(soft_assign, k)
    priors = imm_cluster.update_priors([1.0/new_k]*new_k, readsf, {}, {}, soft_assign)
    (likelihood,read_probs) = imm_cluster.get_read_probs(priors, {}, {}, soft_assign)
    return likelihood


############################################################
# minentropy_clusters
#
# Copy the clustering with minimum entropy to the main
# directory.
############################################################
def minentropy_clusters(total_starts, readsf, k, soft_assign):
    entropy = [0]*total_starts
    for i in range(total_starts):
        os.chdir('tmp.start%d' % i)
        if len(glob.glob('cluster-*.fa')) > 0:
            # determine likelihood
            entropy[i] = get_entropy(readsf, k, soft_assign)
        else:
            # something failed
            entropy[i] = ''
        os.chdir('..')

    # find min entropy partitioning ('' is greater than numbers)
    (min_entropy, min_clust) = util.min_i(entropy)

    # get files from min
    os.system('cp tmp.start%d/cluster-*.fa tmp.start%d/icm-*scores.tmp .' % (min_clust,min_clust))


############################################################
# get_entropy
#
# Return the entropy of the clusters in the current
# directory.
############################################################
def get_entropy(readsf, k, soft_assign):
    new_k = determine_k(soft_assign, k)
    priors = imm_cluster.update_priors([1.0/new_k]*new_k, readsf, {}, {}, soft_assign)
    (like, read_probs) = imm_cluster.get_read_probs(priors, {}, {}, soft_assign)

    entropy = 0.0
    for r in read_probs:
        for c in range(len(read_probs[r])):
            if read_probs[r][c] > 0:
                entropy += -read_probs[r][c]*math.log(read_probs[r][c])

    return entropy

############################################################
# determine_k
#
# In case, I'm letting k change within LikelyBin
############################################################
def determine_k(soft_assign, k):
    new_k = 0
    for i in range(k):
        if soft_assign:
            f = 'cluster-%d.build.fa' % i
        else:
            f = 'cluster-%d.fa' % i

        if os.path.isfile(f) and os.path.getsize(f) > 0:
            new_k += 1

        # or if job is done
        elif os.path.isfile(f+'.headers') and os.path.getsize(f+'.headers') > 0:
            new_k += 1

    return new_k

    
############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
