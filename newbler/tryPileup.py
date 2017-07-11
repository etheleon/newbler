#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import pprint
import multiprocessing as mp

from .pileup import Alignment
#from newbler.newbler import Newbler
pp = pprint.PrettyPrinter(indent = 4)

parser = argparse.ArgumentParser(description="Generates Pileup and forms second", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--root', metavar='ROOT', dest='root', type=str, default = os.getcwd(),
help="""Path to root directory, Directory structure such the following:
<your project root directory>
└── out
    ├── newbler
    │   └── K0000X
    │       └── input
    │           └── K0000X.1.fq (binned reads from Diamond+blast2lca)
    │           └── K0000X.2.fq (binned reads from Diamond+blast2lca)
    ├── pileup
    │   └── K0000X
    │       └── input
    │           └── K0000X-contig00001
    │           └── K0000X-contig00002
    ├── preNewbler
    │   └── K0000X
    │       └── K0000X (fastQ file)
    ├── pAss01
    ├── pAss03
    ├── pAss05
    ├── pAss10
    └── pAss11
""")
parser.add_argument('--pool', metavar='pool', dest="pool",type=int, default = 1, help="The number individual newbler sessions")
parser.add_argument('--cpu', metavar='CPU', dest="cpu",type=int, default = 1, help="total number of threads")
parser.add_argument('--subset', metavar='N', dest='subset', type=int, nargs=2, help="Run the script for a subset of KOs")

args = parser.parse_args()
pp.pprint(args)

kos = os.listdir("%s/out/newbler" % args.root)
if args.subset is not None:
    kos = kos[args.subset[0]:args.subset[1]]
#KOS = ['K00927']
print(kos)
def preAssembly(root, koid):
    """
    calls object to generate output
    """
    alignment = Alignment(root, koid)
    alignment.doPile()
    alignment.getReadsFromPileUP() #needs out/preNewbler

##################################################
def callback(response):
    print("Done - Error:", response)

def err_call(response):
    print("Done - error:", response)

POOL = mp.Pool(processes=args.cpu)
for ko in kos:
    POOL.apply_async(preAssembly, args=(args.root, ko), callback=callback, error_callback=err_call)
POOL.close()
POOL.join()
print("Outputs:")
print("\tPileup: %s/out/pileup/" % args.root)
print("\tExtracted reads: %s/out/preNewbler" % args.root)
##################################################

#def runAssembly(root, koid):
#    """
#    runsAssembly
#    """
#    alignment = Alignment(root, koid)
#    try:
#        os.mkdir("out/assemble/%s"%koid)
#        print("dir out/assemble already exists")
#    except OSError as e:
#        if e.errno != errno.EEXIST:
#            raise  # raises the error again
#    alignment.roundTwoAssembly()


#def runAssembly(root, koid):
#    """
#    calls object to generate output
#    """
#    #alignment = Alignment(root, koid)
#    newbler = Newbler(root, ko, cpu)
#
#POOL = mp.Pool(processes=2)
#print("Assembling now")
#for ko in KOS:
#    POOL.apply_async(runAssembly, args=(ROOT, ko), callback=callback, error_callback=err_call)
#POOL.close()
#POOL.join()
#print("Done Assembling")
