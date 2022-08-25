#!/bin/sh
mkdir /Users/marlena/Documents/FileZilla/rct/$1

scp mnorwood@bayes.biostat.washington.edu:/home/students/mnorwood/rct/$1/\{summary.csv,DESCRIPTION.txt\} /Users/marlena/Documents/FileZilla/rct/$1

