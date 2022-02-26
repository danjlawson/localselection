#!/bin/sh
# fs test.cp -N -v -phasefiles test.chr1.phase test.chr2.phase -recombfiles test.chr1.recomb test.chr2.recomb -idfile test.ids -s2args:-b\ -k\ 10 -s1indfrac 0.1 -indsperproc 1 -go
fs test.cp -N -v -phasefiles test.chr1.phase test.chr2.phase -recombfiles test.chr1.recomb test.chr2.recomb -popidfile test.ids -s7args:-b -s6indfrac 0.1 -indsperproc 1 -hpc 1 -go
cat test/commandfiles/commandfile6.txt | parallel --bar
fs test.cp -go
cat test/commandfiles/commandfile7.txt | parallel --bar

