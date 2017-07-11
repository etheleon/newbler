# NEWBLER Manipulation for Gene Centric Analyses

The scripts and classes are tested using Anaconda installation of python 3.5,
Requires NEWBLER2.9 installed on system.

## Project Folder structure

```
/your/project/root/directory
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
```

## Docker

```
docker pull etheleon/python3
```

Run interactively

```
cd $projectFolder

docker run --rm -u 507 -v $PWD:/w -w /w -it etheleon/python3:0.1 /bin/bash
```



## Description

Much of the gene centric pipeline is only are scripts. pASS is just one part of this.

Classes

2. Newbler
    * Newbler._geneCentricAssembly_ - runs Newbler2.9 assembler on a KO by KO basis, generating contigs for use later for pAss
    * Newbler._mdrCentricAssembly_  - runs Newbler2.9 assembler on a KO by KO basis, but only for READs found in the MDR region, generating contigs for use later for pAss

1. Alignment
    Alignment._doPile_ -         Generates a short read pileup for each of the contigs to be used later for assembly.

pileup 

