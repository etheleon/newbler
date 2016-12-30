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

    def __init__ (self, root, ko, cpu):
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

    def geneCentricAssembly(self, debug=False):
        '''
        Running assembler: NEWBLER first time to generate gene centric assemblies
        '''
        inputFile = "%s/%s/input/%s" % (self.root, self.ko, self.ko)

        def check(fq, num):
            fqExists = os.path.isfile(fq)
            if fqExists:
                num_lines = sum(1 for line in open(fq))
                return fqExists and num_lines > 0
            else:
                print("fastQ read%s does not exists for %s"%(num, self.ko))
                return False

        for i in ("1","2"):
            self.info['fq'+i] = {
                'filePath' : "%s.%s.fq"%(inputFile, i)
            }
            self.info['fq'+i]['status'] = check(self.info['fq'+i]['filePath'], i)

        haveBothReads = (self.info['fq1']['status'] and self.info['fq2']['status'])
        onlyHaveRead1 = (self.info['fq1']['status'] and  not self.info['fq2']['status'])
        onlyHaveRead2 = (not self.info['fq1']['status'] and self.info['fq2']['status'])

        assm = "/home/uesu/Downloads/newbler/opt/454/apps/mapper/bin/runAssembly"
        headCMD = assm + " -cpu " + self.cpu + " -force -m -urt -rip -o %s/%s" % (self.root, self.ko)
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
