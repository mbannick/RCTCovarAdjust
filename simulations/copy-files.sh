#!/bin/sh
mkdir /Users/marlena/Documents/FileZilla/rct/$1

scp mnorwood@bayes0.biostat.washington.edu:/home/students/mnorwood/rct/$1/DESCRIPTION.txt /Users/marlena/Documents/FileZilla/rct/$1

scp mnorwood@bayes0.biostat.washington.edu:/home/students/mnorwood/rct/$1/summary.csv /Users/marlena/Documents/FileZilla/rct/$1
