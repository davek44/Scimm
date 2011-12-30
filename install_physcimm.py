#!/usr/bin/env python
import os, sys, subprocess

############################################################
# install_physcimm.py
#
# Compile source and set paths for Physcimm.
#
# Author: David Kelley
############################################################

prior_phymm_dir = ''

############################################################
# main
############################################################
def main():
    status = raw_input('Please note that this will install Phymm,\nwhich will likely take ~24 hours and use\n~50 GB of space. If you already have Phymm\ninstalled, edit this script and set\n"prior_phymm_dir" to the package directory\npath.\nContinue? [y/n] ')
    if status == 'n' or status == 'N':
        exit()

    installdir = os.getcwd()

    # compile IMM code
    os.chdir('glimmer3.02/src')
    p = subprocess.Popen('make clean; make', shell=True)
    os.waitpid(p.pid,0)

    os.chdir('../..')

    # set IMM links
    if not os.path.isfile('bin/simple-score'):
        os.symlink('../glimmer3.02/bin/simple-score','bin/simple-score')
    if not os.path.isfile('bin/build-icm'):
        os.symlink('../glimmer3.02/bin/build-icm', 'bin/build-icm')
    
    # set scimm bin variable
    p = subprocess.Popen('sed \'s,scimm_bin = ".*",scimm_bin = "%s/bin",\' bin/scimm.py > sc.tmp' % installdir, shell=True)
    os.waitpid(p.pid, 0)
    os.rename('sc.tmp', 'bin/scimm.py')
    p = subprocess.Popen('chmod ug+x bin/scimm.py', shell=True)
    os.waitpid(p.pid,0)

    # set physcimm bin variable
    p = subprocess.Popen('sed \'s,phymmdir = ".*",phymmdir = "%s/phymm",\' bin/physcimm.py > ph.tmp' % installdir, shell=True)
    os.waitpid(p.pid, 0)
    os.rename('ph.tmp', 'bin/physcimm.py')
    p = subprocess.Popen('chmod ug+x bin/physcimm.py', shell=True)
    os.waitpid(p.pid,0)

    if prior_phymm_dir:
        if os.path.islink(prior_phymm_dir):
            os.remove(prior_phymm_dir)
        os.symlink(prior_phymm_dir,'phymm')
    else:
        os.mkdir('phymm')
        os.chdir('phymm')
        p = subprocess.Popen('curl -o phymmInstaller.tar.gz http://www.cbcb.umd.edu/software/phymm/phymmInstaller.tar.gz', shell=True)
        os.waitpid(p.pid, 0)
        p = subprocess.Popen('tar -xzvf phymmInstaller.tar.gz', shell=True)
        os.waitpid(p.pid, 0)
        p = subprocess.Popen('./phymmSetup.pl', shell=True)
        os.waitpid(p.pid, 0)


############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
