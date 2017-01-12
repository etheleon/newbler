#!/usr/bin/env python

import  os
import  errno
import  multiprocessing as mp

from newbler.pileup import Alignment

ROOT      = "./"
DIRECTORY = "./out/newbler"
KOS       = os.listdir(DIRECTORY)

#KOS = ['K00927']

def extract_MDR_Reads(root, koid):
    """
    calls object to generate output
    """
    alignment = Alignment(root, koid)
    try:
        os.mkdir("out/pileup")
        print("dir out/pileup already exists")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise  # raises the error again
    try:
        os.mkdir("out/preNewbler/%s"%koid)
        print("dir out/preNewbler already exists")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise  # raises the error again
    try:
        os.mkdir("out/assemble/%s"%koid)
        print("dir out/assemble already exists")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise  # raises the error again
    alignment.doPile()
    alignment.getReadsFromPileUP()
    #alignment.roundTwoAssembly()

def runAssembly(root, koid):
    """
    runsAssembly
    """
    alignment = Alignment(root, koid)
    try:
        os.mkdir("out/assemble/%s"%koid)
        print("dir out/assemble already exists")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise  # raises the error again
    alignment.roundTwoAssembly()


def runAssembly(root, koid):
    """
    calls object to generate output
    """
    alignment = Alignment(root, koid)

def err_call(response, ko):
    """
    error callback function
    """
    print("Done - error:", response)

def callback(response, ko):
    """"
    callback function
    """
    print("Done - Error:", response)

POOL = mp.Pool(processes=20)
for ko in KOS:
    print("Extracting %s MDR reads for 2nd Assembly now" % ko)
    POOL.apply_async(extract_MDR_Reads, args=(ROOT, ko), callback=callback, error_callback=err_call)
POOL.close()
POOL.join()
print("Done extracting")

POOL = mp.Pool(processes=2)
print("Assembling now")
for ko in KOS:
    POOL.apply_async(runAssembly, args=(ROOT, ko), callback=callback, error_callback=err_call)
POOL.close()
POOL.join()
print("Done Assembling")
