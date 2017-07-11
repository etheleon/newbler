#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import shutil
import subprocess

class Newbler:
    '''
    For running newbler, overlap group assembler, for all KOs.
    This is to be used together with a specific
    docker image: etheleon/python3
    '''

    def __init__(self, root, ko, cpu, assm = "/home/uesu/Downloads/newbler/opt/454/apps/mapper/bin/runAssembly"):
        self.assm = assm
        self.root = root #the directory which contains the KOs
        self.ko = ko
        self.cpu = cpu
        self.info = {}

    def __cleanup(self):
        '''
        removes content in the KO's directory save for the
        fastQ files
        '''
        files = os.listdir("%s/%s" % (self.root,self.ko))
        offending = ["%s/%s/%s"%(self.root, self.ko,file) for file in files]
        for nonsense in list(filter(lambda file: not re.search('input', file), offending)):
            if os.path.isfile(nonsense):
                os.remove
            else:
                shutil.rmtree(nonsense)

    def __checkNewblerIsDone(self):
        '''
        newbler fails to run without error, this checks
        454NewblerProgress for completion
        '''
        try:
            filePath = "%s/%s/454NewblerProgress.txt" % (self.root, self.ko)
            with open(filePath) as progressFile:
                content = progressFile.read()
                if re.search("Assembly computation succeeded", content):
                    return True
                else:
                    return False
        except IOError as e:
            print("I/O error({0}): {1}: {2}".format(e.errno, e.strerror, filePath))
            return False

    def __check(self, fq, num):
        fqExists = os.path.isfile(fq)
        if fqExists:
            num_lines = sum(1 for line in open(fq))
            return fqExists and num_lines > 0
        else:
            print("fastQ read%s does not exists for %s"%(num, self.ko))
            return False


    def geneCentricAssembly(self, debug=False, MDR=True, timeoutlimit=7200):
        '''
        Running assembler: NEWBLER first time to generate gene centric assemblies
        '''
        #times out after 2 hours
        if MDR:
            inputFile = "%s/%s/%s" % (self.root, self.ko, self.ko)
        else:
            inputFile = "%s/%s/input/%s" % (self.root, self.ko, self.ko)

        for i in ("1","2"):
            self.info['fq'+i] = {
                'filePath' : "%s.%s.fq"%(inputFile, i)
            }
            self.info['fq'+i]['status'] = self.__check(self.info['fq'+i]['filePath'], i)

        haveBothReads = (self.info['fq1']['status'] and self.info['fq2']['status'])
        onlyHaveRead1 = (self.info['fq1']['status'] and  not self.info['fq2']['status'])
        onlyHaveRead2 = (not self.info['fq1']['status'] and self.info['fq2']['status'])

        headCMD = self.assm + " -cpu " + self.cpu + " -force -m -urt -rip -o %s/%s" % (self.root, self.ko)
        if haveBothReads:
            cmd = headCMD + " %s.1.fq %s.2.fq" % (inputFile, inputFile)
        elif onlyHaveRead1:
            cmd = headCMD + " %s.1.fq" % inputFile
        elif onlyHaveRead2:
            cmd = headCMD + " %s.2.fq" % inputFile
        else:
            print("Not processing: Both have no reads")
            cmd = ""
        result = None
        count = 0
        if debug:
            print("Cmd: %s" % cmd)
            print(self.info)
        else:
            print("Assembling %s..." % self.ko)
            print("executing: %s" % cmd)
            while result is None:
                count = count + 1
                if (count > 6):
                    print("Exceeded 5 tries, not going to try again")
                    result = True
                try:
                    #doesnt really capture newbler's error messages
                    subprocess.run(cmd, shell=True, check=True, timeout=timeoutlimit) #timeout 3600 seconds ie. 1 hour. so wait for max 2 hours before we timeout
                except subprocess.CalledProcessError as err:
                    #error is empty
                    print("newbler error:\n", err.output)
                    self.__cleanup(self.ko, self.root)
                except subprocess.TimeoutExpired as err:
                    print("Assembly of %s took more than 2 hours. Aborting" % ko)
                    self.__cleanup(self.ko, self.root)
                else:
                    if self.__checkNewblerIsDone():
                        result = True
                    else:
                        pass
            print("Done Assembling")

    def mdrCentricAssembly(self, debug=False):
        '''
        Running assembler: NEWBLER second time to generate gene centric assemblies
        '''

        self.info['fq']['filePath'] = "%s/%s/%s" % (self.root, self.ko, self.ko)
        self.info['fq']['status']   = self.__check(self.info['fq']['filePath'], 1)

        if self.info['fq']['status']:
            cmd = self.assm + " -cpu " + self.cpu + " -force -m -urt -rip -o %s/%s %s" % (self.root, self.ko, self.info['fq']['filePath'])
        else:
            print("Not processing: Both have no reads")
            cmd = ""

        result = None
        if debug:
            print("Cmd: %s" % cmd)
            print(self.info)
        else:
            print("Assembling %s..." % self.ko)
            print("executing: %s" % cmd)
            while result is None:
                try:
                    #doesnt really capture newbler's error messages
                    subprocess.run(cmd, shell=True, check=True)
                except subprocess.CalledProcessError as err:
                    #error is empty
                    print("newbler error:\n", err.output)
                    self.__cleanup(self.ko, self.root)
                else:
                    if self.__checkNewblerIsDone():
                        result = True
                    else:
                        pass
            print("Done Assembling")
