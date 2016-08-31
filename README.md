# dupchk

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
Code author: James B. Henderson, jhenderson@calacademy.org  
README.md authors: Zachary R. Hanna, James B. Henderson  
