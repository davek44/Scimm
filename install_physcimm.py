#!/usr/bin/env python
import os, sys

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

    # IMM code
    os.chdir('glimmer3.02/src')
    os.system('make clean; make')
    os.chdir('../..')
    if not os.path.isfile('bin/simple-score'):
        os.chdir('bin')
        os.system('ln -s ../glimmer3.02/bin/simple-score')
        os.chdir('..')
    if not os.path.isfile('bin/build-icm'):
        os.chdir('bin')
        os.system('ln -s ../glimmer3.02/bin/build-icm')
        os.chdir('..')
    
    # Scimm
    os.system('sed -i .bak \'s,scimm_bin = "[a-zA-Z/]*",scimm_bin = "%s/bin",\' bin/scimm.py' % installdir) 
    os.system('rm bin/scimm.py.bak')

    # Phymm
    os.system('sed -i .bak \'s,phymmdir = "[a-zA-Z/]*",phymmdir = "%s/phymm",\' bin/physcimm.py' % installdir) 
    os.system('rm bin/physcimm.py.bak')
    if prior_phymm_dir:
        os.system('ln -s %s phymm' % prior_phymm_dir)
    else:
        os.mkdir('phymm')
        os.chdir('phymm')
        os.system('curl -o phymmInstaller.tar.gz http://www.cbcb.umd.edu/software/phymm/phymmInstaller.tar.gz')
        os.system('tar -xzvf phymmInstaller.tar.gz')
        os.system('./phymmSetup.pl')
        os.chdir('..')


############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
