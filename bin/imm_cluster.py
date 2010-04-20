#!/usr/bin/env python

from __future__ import division
from optparse import OptionParser
import sys, os, glob, random, math, util, pdb
import scimm

############################################################
# imm_cluster.py
#
# Perform IMM clustering on a set of seuqneces.
############################################################

max_iter = 200
use_priors = True
like_decrease_t = 5
rsments_t = .0005
soft_assign_t = .01

############################################################
# main
############################################################
def main():
    parser = OptionParser()

    parser.add_option('-k', dest='k', type='int', help='Number of clusters')
    parser.add_option('-d','--read_dir', dest='reads_dir', help='Directory with read files')
    parser.add_option('-r', '--read_file', dest='reads_file', help='Fasta file with reads')
    parser.add_option('-p','--par', dest='par', type='int', default=1, help='Parallelize with PAR cpus')
    parser.add_option('-m','--mates', dest='mates_file', help='Read mate pairs')
    parser.add_option('-c','--constraints', dest='constraints_file', help='Read cluster constraints')
    parser.add_option('-i','--initial', dest='initial_done', action='store_true', default=False, help='Initial partition is given')
    parser.add_option('-s','--seed', dest='seed', action='store_true', default=False, help='Incomplete initial partition is given')
    parser.add_option('--seed_only', dest='seed_only', action='store_true', default=False, help='Perform a single iteration of the algorithm using a seeded initialization')
    parser.add_option('--trained', dest='trained', action='store_true', default=False, help='The models are already trained for the first iteration (e.g. by --seed_only')
    parser.add_option('--em', dest='soft_assign', action='store_true', default=False, help='Use a soft assignment of sequences to clusters and use expectation maximization')

    (options, args) = parser.parse_args()

    if not options.k:
        parser.error('Must define k')
    if not options.reads_file and not options.reads_dir:
        parser.error('Must provide reads')

    num_reads = 0
    for line in open(options.reads_file):
        if line[0] == '>':
            num_reads += 1

    # in case k shrinks
    k = options.k

    # load mate pairs
    mates = load_mates(options.mates_file)

    # load constraints
    constraints = load_constraints(options.constraints_file)
    
    # partition reads into k files
    if options.seed:
        if options.constraints_file:
            constraint_seed(options.reads_file, k, mates, constraints, options.soft_assign)
        (like,priors) = seed_partition(options.reads_file, k, mates, constraints, options.soft_assign, options.par)
        k = filter_empty(k, priors, constraints)
        print 'Iter 0:\t%d' % int(like)

        if options.seed_only:
            # train on all sequences for fair likelihood comparisons
            train_imm(k, options.soft_assign, options.par)
            # score each read with each IMM
            score_reads(k, options.reads_file, options.par)
            exit()
    else:
        priors = [1.0/k]*k            

    if options.initial_done or options.seed:
        verify_constraints(k, constraints)
        
    else:
        random_partition(options.reads_file, options.reads_dir, k, mates, options.soft_assign)
        #correct_partition_cb(options.reads_file, k, options.soft_assign)

    # progress data
    prog = Progress(k)
    good_prog = True
    rsments = num_reads
    iter = 0
    
    while iter < max_iter and rsments >= num_reads*rsments_t and good_prog:
        iter += 1
                
        if iter > 1 or not options.trained:
            # train an IMM on each cluster
            train_imm(k, options.soft_assign, options.par)
        
            # score each read with each IMM
            score_reads(k, options.reads_file, options.par)

        # reassign reads to max scoring IMM
        (rsments,like,priors) = reassign_reads(options.reads_file, priors, mates, constraints, options.soft_assign, False)
        k = filter_empty(k, priors, constraints)

        print 'Iter %d:\t%d\t%d reassignments' % (iter,int(like),rsments)

        good_prog = prog.assess(like)

    # take max
    for i in range(k):
        os.system('mv cluster-%d.max cluster-%d.fa' % (i,i))

    #os.system('rm *.tmp cluster-*.build.fa')

############################################################
# train_imm
#
# Train an IMM on each cluster.  The reads are in 
# cluster-#.fa, and the IMM will be in cluster-#.icm
############################################################
def train_imm(k, soft_assign, par):
    cmds = []
    for i in range(k):
        if soft_assign:
            cmds.append('%s/em_build-icm -p 1 cluster-%d.icm < cluster-%d.build.fa' % (scimm.scimm_bin,i,i))
        else:
            cmds.append('%s/build-icm -p 1 cluster-%d.icm < cluster-%d.fa' % (scimm.scimm_bin,i,i))

    util.exec_par(cmds, par)

############################################################
# score_reads
#
# For each IMM, run on each cluster of reads outputting the
# scores to a temp file cluster-#.icm-#.scores.tmp
############################################################
def score_reads(k, readsf, par):
    cmds = []
    for c in range(k):
        cmds.append('%s/simple-score -N cluster-%d.icm < %s > icm-%d.scores.tmp 2>/dev/null' % (scimm.scimm_bin,c,readsf,c))
    
    util.exec_par(cmds, par)

############################################################
# reassign_reads
#
# For each cluster of reads, check their IMM scores and
# assign each read to a new cluster.  Also, calculate the
# likelihood of the reads (in their current clusters)
# given the current model.
############################################################
def reassign_reads(readsf, priors, mates, constraints, soft_assign, initial_seed):
    k = len(priors)

    if use_priors:
        priors = update_priors(priors, readsf, mates, constraints, soft_assign)

    (likelihood, read_probs) = get_read_probs(priors, mates, constraints, soft_assign)

    # open files
    read_files = []
    build_files = []
    for i in range(k):
        read_files.append(open('cluster-%d.tmp' % i,'w'))
        if soft_assign:
            build_files.append(open('cluster-%d.build.fa' % i, 'w'))
    
    if initial_seed:
        myk = 1
    else:
        myk = k

    rsments = 0
    for c in range(myk):
        # reassign
        for line in open('cluster-%d.fa' % c):
            if line[0] == '>':
                r = line.split()[0][1:]
                if not read_probs.has_key(r):
                    print 'ERROR: missing read %s scores' % r
                    exit()
                    
                elif constraints.has_key(r):
                    if constraints[r] != c:
                        print 'Found a constrained read in the wrong cluster'
                    max_icm = constraints[r]

                else:
                    (max_prob, max_icm) = util.max_i(read_probs[r])

                    # count reassignments
                    if max_icm != c:
                        rsments += 1

            # print line to files
            read_files[max_icm].write(line)
            if soft_assign:
                for i in range(k):
                    if read_probs[r][i] > soft_assign_t:
                        if line[0] == '>':
                            build_files[i].write('>%f;%s' % (read_probs[r][i],line[1:]))
                        else:
                            build_files[i].write(line)
            
    # close files
    for i in range(k):
        read_files[i].close()
        if soft_assign:
            build_files[i].close()

    # move tmp
    for i in range(k):
        os.system('mv cluster-%d.tmp cluster-%d.fa' % (i,i))

    return (rsments,likelihood,priors)

############################################################
# log_add
#
# Safely compute log(e^l_i + e^l_j)
############################################################
def log_add(l_i, l_j):
    if l_i > l_j:
        return l_i + math.log(1 + math.exp((l_j - l_i)))
    else:
        return l_j + math.log(1 + math.exp((l_i - l_j)))
    

############################################################
# update_priors
#
# Calculate the proportion of sequence in each cluster
# to be used as the prior probability of a read being
# from that cluster.
############################################################
def update_priors(prev_priors, readsf, mates, constraints, soft_assign):    
    k = len(prev_priors)

    # get read probs
    (l,read_probs) = get_read_probs(prev_priors, mates, constraints, soft_assign)

    # count expected bp
    exp_bp = [0]*k
    for line in open(readsf):
        if line[0] == '>':
            r = line[1:].rstrip()
        else:
            for c in range(k):
                exp_bp[c] += read_probs[r][c]*len(line.rstrip())

    # normalize
    total_bp = sum(exp_bp)
    priors = [bp/total_bp for bp in exp_bp]

    return priors

############################################################
# get_read_probs
#
# Given a set of priors, return the overall likelihood
# and a dict of read probabilities
############################################################
def get_read_probs(priors, mates, constraints, soft_assign):
    k = len(priors)

    read_likes = {}
    for c in range(k):
        for line in open('icm-%d.scores.tmp' % c):
            (r,s) = line.split()
            if not read_likes.has_key(r):
                read_likes[r] = [0]*k
            read_likes[r][c] = float(s)
            
    read_probs = {}
    likelihood = 0.0
    for r in read_likes:
        # if constrained, set prob accordingly
        if constraints.has_key(r):
            read_probs[r] = [0]*k
            read_probs[r][constraints[r]] = 1.0

            # I don't care about their likelihood

        else:
            # combine mate likelihoods and priors
            if mates.has_key(r):
                r1 = read_likes[r]
                r2 = read_likes[mates[r]['mate']]
                #read_scores = [r1[x]+r2[x]+priors[x] for x in range(k)]
                read_scores = [r1[x]+r2[x]+math.log(priors[x]) for x in range(k)]
            else:
                r1 = read_likes[r]
                #read_scores = [r1[x]+priors[x] for x in range(k)]
                read_scores = [r1[x]+math.log(priors[x]) for x in range(k)]

            # determine probabilities of assignments
            sum_score = read_scores[0]
            for i in range(1,k):
                sum_score = log_add(sum_score,read_scores[i])
            read_probs[r] = []
            for i in range(k):
                read_probs[r].append(math.exp(read_scores[i] - sum_score))

            # update likelihood, accounting for mates being assigned twice
            if mates.has_key(r):
                if soft_assign:
                    likelihood += sum_score/2
                else:
                    likelihood += max(read_scores)/2
            else:
                if soft_assign:
                    likelihood += sum_score
                else:
                    likelihood += max(read_scores)

    return (likelihood, read_probs)


############################################################
# random_partition
#
# The input is a directory name where there will be a bunch
# of read fasta files.  Split them up randomly into k files
# called cluster-#.fa
############################################################
def random_partition(reads_file, reads_dir, k, mates, soft_assign):
    # open files
    rand_read_files = []
    build_files = []
    for c in range(k):
        rand_read_files.append(open('cluster-%d.fa' % c,'w'))
        if soft_assign:
            build_files.append(open('cluster-%d.build.fa' % c,'w'))

    # check for read files
    if reads_file:
        read_files = [reads_file]
    elif reads_dir:        
        read_files = glob.glob(reads_dir+'/*.fa')
        
    if not read_files:
        print reads_dir
        print 'No read files.'
        exit()
        
    # split reads into new files
    for readf in read_files:
        for line in open(readf):
            if line[0] == '>':
                # keep mates together
                header = line.split()[0][1:]
                if mates.has_key(header):
                    m = mates[header]

                    # assign to mates cluster
                    if m['cluster'] != -1:
                        rf = m['cluster']
                    else:
                        # or to random (and save)
                        rf = random.randint(0,k-1)
                        m['cluster'] = rf
                        mates[m['mate']]['cluster'] = rf
                        
                # or to random
                else:
                    rf = random.randint(0,k-1)
                    
            rand_read_files[rf].write(line)
            if soft_assign:
                if line[0] == '>':
                    build_files[rf].write('>1.0;%s' % line[1:])
                else:
                    build_files[rf].write(line)
        
    # close files
    for c in range(k):
        rand_read_files[c].close()
        if soft_assign:
            build_files[c].close()

    # get back mate 'cluster' memory
    for r in mates:
        del mates[r]['cluster']
        

############################################################
# correct_partition_cb
#
# The input is a directory name where there will be a bunch
# of read fasta files.  Just move them as is to cluster-#.fa
############################################################
def correct_partition_cb(reads_file, k, soft_assign):
    # open files
    rand_read_files = []
    build_files = []
    for c in range(k):
        rand_read_files.append(open('cluster-%d.fa' % c,'w'))
        if soft_assign:
            build_files.append(open('cluster-%d.build.fa' % c,'w'))
        
    # split reads into new files
    for line in open(reads_file):
        if line[0] == '>':
            species = line.split('|')[0]
            sp = int(species[species.find('_')+1:])

        # print
        rand_read_files[sp].write(line)
        if soft_assign:
            if line[0] == '>':
                build_files[sp].write('>1.0;%s' % line[1:])
            else:
                build_files[sp].write(line)
                


############################################################
# correct_partition_simXC
#
# The input is a directory name where there will be a bunch
# of read fasta files.  Just move them as is to cluster-#.fa
#
# Doesn't work properly with EM
############################################################
def correct_partition_simXC(reads_dir, k):
    # check for read files
    read_files = glob.glob(reads_dir+'/*.fa')
    if not read_files:
        print 'No read files.'
        exit()
    # split reads into new files
    cf = 0
    for readf in read_files:
        os.system('cp %s cluster-%d.fa' % (readf,cf))
        cf += 1
        

############################################################
# seed_partition
#
# To initialize the algorithm, train IMM's on incomplete
# read clusters and do a maximization step to partition
# the remainder of the reads. (Actually the seed reads
# can be moved as well which I think is ok.)
############################################################
def seed_partition(readsf, k, mates, constraints, soft_assign, par):
    # train IMMs
    train_imm(k, soft_assign, par)

    # score all reads
    score_reads(k, readsf, par)

    # check scores and partition
    os.system('cp %s cluster-0.fa' % readsf)
    (rsments, likelihood, priors) = reassign_reads(readsf, [1.0/k]*k, mates, constraints, soft_assign, True)

    return(likelihood,priors)

############################################################
# constraint_seed
#
# Seed the algorithm's initial partitioning with the
# constraint reads.
############################################################
def constraint_seed(reads_file, k, mates, constraints, soft_assign):
    # open files
    seed_read_files = []
    build_files = []
    seed_reads = [0]*k    
    for i in range(k):
        seed_read_files.append(open('cluster-%d.fa' % i,'w'))
        if soft_assign:
            build_files.append(open('cluster-%d.build.fa' % i,'w'))

    # parse reads file, printing constrained reads
    for line in open(reads_file):
        if line[0] == '>':
            header = line[1:].rstrip()
            if constraints.has_key(header):
                c = constraints[header]
                seed_reads[c] += 1
            else:
                c = -1
                
        if c != -1:
            seed_read_files[c].write(line)
            if soft_assign:
                if line[0] == '>':
                    build_files[c].write('>1.0;%s' % line[1:])
                else:
                    build_files[c].write(line)

    # close files
    for c in range(k):
        seed_read_files[c].close()
        if soft_assign:
            build_files[c].close()
    
    # check seeds all > 0
    for c in range(k):
        if seed_reads[c] == 0:
            print 'Cluster %d has no seed reads.  Try a different initialization method.' % c
            exit()        



############################################################
# verify_constraints
#
# My constraints map reads to certain cluster #'s so I may
# need to swap some files if the initial clusters are
# given.
############################################################
def verify_constraints(k, constraints):
    cluster_map = {}    
    for c in range(k):
        for line in open('cluster-%d.fa' % c):
            if line[0] == '>':
                header = line[1:].rstrip()

                # if read is constrainted
                if constraints.has_key(header):
                    # either initialize cluster map
                    if not cluster_map.has_key(c):
                        cluster_map[c] = constraints[header]

                    # or verify that it matches
                    else:
                        if cluster_map[c] != constraints[header]:
                            print 'Inconsistent constraints: cluster %d' % c
                            exit()

    # handle unconstrained clusters
    # by finding clusters that aren't mapped to
    open_clusters = []
    for c in range(k):
        if c not in cluster_map.values():
            open_clusters.append(c)

    # and mapping unconstrained clusters to them
    i = 0
    for c in range(k):
        if not cluster_map.has_key(c):
            cluster_map[c] = open_clusters[i]
            i += 1

    # move clusters to their matching constraint number
    for c in range(k):
        os.system('mv cluster-%d.fa cluster-%d.tmp.fa' % (c,cluster_map[c]))
    for c in range(k):
        os.system('mv cluster-%d.tmp.fa cluster-%d.fa' % (c,c))
                  
                            

############################################################
# filter_empty
#
# Filter out empty cluster files
############################################################
def filter_empty(k, priors, constraints):
    # find empty clusters
    empty = []
    for c in range(k):
        if os.path.getsize('cluster-%d.fa' % c) == 0:
            empty.append(c)

    while empty:
        # choose cluster to delete
        c = empty[0]
        empty = empty[1:]

        # rename all following clusters by 1, update priors
        for i in range(c+1,k):
            os.system('mv cluster-%d.fa cluster-%d.fa' % (i,i-1))
            priors[i-1] = priors[i]

        # update constraints
        for r in constraints:
            if constraints[r] > c:
                constraints[r] -= 1

        # update empty
        for i in range(len(empty)):
            empty[i] -= 1

        # update k
        k -= 1

    # re-normalize priors
    sp = sum(priors)
    priors = [p/sp for p in priors]

    return k

############################################################
# load_mates
#
# Load mates pairs from file
############################################################
def load_mates(matef):
    mates = {}
    if matef:
        for line in open(matef):
            (r1,r2) = line.split()
            mates[r1] = {'mate':r2, 'cluster':-1, 'scores':[]}
            mates[r2] = {'mate':r1, 'cluster':-1, 'scores':[]}
    return mates

############################################################
# load_constraints
#
# Load constraints from file
############################################################
def load_constraints(constrainf):
    constraints = {}
    if constrainf:
        for line in open(constrainf):
            (r,c) = line.split()
            constraints[r] = int(c)
    return constraints


############################################################
# Progress
#
# A class to help manage convergence progress
############################################################
class Progress:
    def __init__(self, k):
        self.last_like = ''
        self.max_like = ''
        self.like_decr = 0
        self.k = k

    def assess(self,like):
        # compare to last likelihood
        if self.last_like and like > self.last_like:
            self.like_decr = 0            
            
        elif self.last_like and like < self.last_like:
            self.like_decr += 1

        self.last_like = like

        # compare to max likelihood
        if not self.max_like or like > self.max_like:
            # save current as max
            for i in range(self.k):
                os.system('cp cluster-%d.fa cluster-%d.max' % (i,i))

        if self.like_decr >= like_decrease_t:
            return False
        else:
            return True
        
############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
