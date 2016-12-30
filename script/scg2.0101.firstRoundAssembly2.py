#!/usr/bin/python

import argparse
from newbler.newbler import Newbler

parser = argparse.ArgumentParser(description='Gene Centric Assembly',
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('start', metavar='k', type=int, nargs='?',
                    help='the starting nth KO')
parser.add_argument('end', metavar='p', type=int, nargs='?',
                    help='the ending nth KO')
parser.add_argument('cpu', metavar='t', type=str, nargs='?',
                    help='the number of cpu cores to use')
parser.add_argument('--newbler', default="out/newbler",help='''
binned KO reads default: ./out/newbler
EXAMPLE:
    out/newbler
    ├── K00001
    │   └── input
    │       ├── K00001.1.fq
    │       ├── K00001.2.fq
'''
)
args = parser.parse_args()

kos = os.listdir("%s" % args.newbler)
for ko in kos[args.start:args.end]:
    newbler = Newbler(args.newbler, ko, args.cpu)
    newbler.geneCentricAssembly()
