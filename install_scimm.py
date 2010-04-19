#!/fs/sz-user-supported/Linux-i686/bin/python2.5
import os

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

    # LikelyBin
    os.chdir('likelybin-0.1')
    os.system('make clean; make')
    os.chdir('..')
    if not os.path.isfile('bin/mcmc.pl'):
        os.chdir('bin')
        os.system('ln -s ../likelybin-0.1/mcmc.pl')
        os.chdir('..')

    # CBCBCompostBin
    os.system('sed -i \'s,cb_bin = "[a-zA-Z/]\+",cb_bin = "%s/CBCBCompostBin",\' CBCBCompostBin/compostbin.py' % installdir)
    if not os.path.isfile('bin/compostbin.py'):
        os.chdir('bin')
        os.system('ln -s ../CBCBCompostBin/compostbin.py')
        os.chdir('..')

    # Scimm
    os.system('sed -i \'s,scimm_bin ="[a-zA-Z/]\+",scimm_bin = "%s/bin",\' bin/scimm.py' % installdir) 
    

############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
