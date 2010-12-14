#!/usr/bin/env python
import os, subprocess, stat

############################################################
# install_scimm.py
#
# Compile source and set paths for Scimm.
#
# Author: David Kelley
############################################################

############################################################
# main
############################################################
def main():
    installdir = os.getcwd()

    # IMM code
    os.chdir('glimmer3.02/src')
    p = subprocess.Popen('make clean; make', shell=True)
    os.waitpid(p.pid,0)
    os.chdir('../..')
    if not os.path.isfile('bin/simple-score'):
        os.symlink('../glimmer3.02/bin/simple-score','bin/simple-score')
    if not os.path.isfile('bin/build-icm'):
        os.symlink('../glimmer3.02/bin/build-icm', 'bin/build-icm')

    # LikelyBin
    os.chdir('likelybin-0.1')
    p = subprocess.Popen('make clean; make', shell=True)
    os.waitpid(p.pid,0)
    os.chdir('..')
    if not os.path.isfile('bin/mcmc.pl'):
        os.symlink('../likelybin-0.1/mcmc.pl','bin/mcmc.pl')

    # CBCBCompostBin
    p = subprocess.Popen('sed \'s,cb_bin = "[a-zA-Z/]*",cb_bin = "%s/CBCBCompostBin",\' CBCBCompostBin/compostbin.py > cb.tmp' % installdir, shell=True)
    os.waitpid(p.pid,0)
    os.rename('cb.tmp', 'CBCBCompostBin/compostbin.py')
    p = subprocess.Popen('chmod ug+x CBCBCompostBin/compostbin.py', shell=True)
    os.waitpid(p.pid,0)
    if not os.path.isfile('bin/compostbin.py'):
        os.symlink('../CBCBCompostBin/compostbin.py', 'bin/compostbin.py')

    # Scimm
    p = subprocess.Popen('sed \'s,scimm_bin = "[a-zA-Z/]*",scimm_bin = "%s/bin",\' bin/scimm.py > sc.tmp' % installdir, shell=True)
    os.waitpid(p.pid,0)
    os.rename('sc.tmp', 'bin/scimm.py')
    p = subprocess.Popen('chmod ug+x bin/scimm.py', shell=True)
    os.waitpid(p.pid,0)
    

############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
