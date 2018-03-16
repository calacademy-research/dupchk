# dupchk
dupchk is a tool to check a sequence file in FASTA or FASTQ format for duplicate sequences
  
### Contents
README.md : this README  
dupchk.go : Go source code  
  
### How to compile
\# We have provided an executable compiled for use on 64-bit Linux systems in each release.  
\# If you need to compile the code for use on other architectures, you will need Go tools.  
\# If you do not have them, download and install Go tools as described here https://golang.org/doc/install  
$ git clone https://github.com/calacademy-research/dupchk.git  
$ cd dupchk  
$ go build dupchk.go  
  
### Usage
$ dupchk  
```
Usage: dupchk [[-n <NtopDups>] [fa]] <filename.fq>
       Checks file for duplicate reads by looking at the first 21 bases and the last 21 bases
       of each read. Those that match this fingerprint are considered dups for dupchk.
       Optional:
           -n <NtopDups> sets how many of the most frequent dups are shown (default: 5).
           -fa writes to stdout a read in fasta format for each of the top dups shown.
```

### Citing

#### Authorship
Code author: James B. Henderson  
README.md authors: <a href="https://orcid.org/0000-0002-0210-7261" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">Zachary R. Hanna</a>, James B. Henderson  

#### Version 1.0.0
[![DOI](https://zenodo.org/badge/67077965.svg)](https://zenodo.org/badge/latestdoi/67077965)
