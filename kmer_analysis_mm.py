#!/usr/bin/python

import re, os, sys, getopt
import itertools
from itertools import combinations

def calc_score(motif, a, c, g, t, glen, obs):
    p = 1
    bpp = {'A': a, 'C': c, 'G': g, 'T': t, 'N': 1, 'M': 0, 'R': 0, 'W': 0, 'S': 0, 'K': 0, 'Y': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0}
    mlen = len(motif)
    bps = list(motif)
    for i in range(mlen):
        bp = bps[i]
        pbp = bpp[bp]
        p = p * pbp
    exp = glen * p
    score = float((obs - exp))/(obs + exp)
    return score, exp

def main(argv):
    inputfile = ''
    outputfile = ''
    motif = ''
    pos = 0
    n_ = 0
    n = dict()
    bplist = list()
    acgt = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'M': 0, 'R': 0, 'W': 0, 'S': 0, 'K': 0, 'Y': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0}
    genome_len = 0
    nchr = 0
    spec = ''
    mm = 0

    try:
        opts, args = getopt.getopt(argv,"hi:o:s:n:m:",["ifile=","ofile=","spec=","n=","mismatchs="])
    except getopt.GetoptError:
        print ('kmer_analysis.py -i <inputfile> -o <outputfile> -s <species name> -n <motif length> -m <mismatches>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('kmer_analysis.py -i <inputfile> -o <outputfile> -s <species name> -n <motif length> -m <mismatches>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-n"):
            n_ = int(arg)
        elif opt in ("-s"):
            spec = arg
        elif opt in ("-m"):
            m = int(arg)
            if m/n_ > 0.5:
                print ('Mismatches cannot be longer than half of the input k-mer length!')
                sys.exit()

    fi = open(inputfile,'r')
    fo = open(outputfile,'w')

    for line in fi:
        lin = line.rstrip('\n')
        line = str(lin)
        x = re.search(">",line)
        if x :
            pos = 0; motif = ""; bplist = list(); nchr = nchr + 1
        else:
            w = list(line)
            for b in w:
                acgt[b.upper()] = acgt[b.upper()] + 1
                pos = pos+1
                genome_len = genome_len + 1

                bplist.append(b)
                if pos >= n_ + 1:
                    bplist.pop(0)

                if pos >= n_:
                    motif = ''.join(bplist).upper()
                    z = re.search("[NMRWSYKBDHV]",motif)
                    if not z:
                        # normal motif
                        if motif in n:
                            n[motif] = n[motif] + 1
                        else:
                            n[motif] = 1
                        # cycle through mismatched forms
                        ms = list()
                        for ii in range(1, m+1):
                            c = combinations(range(len(motif)),ii)
                            for jj in c:
                                mot = motif
                                for kk in jj:
                                    mot = mot[:kk]+'N'+mot[kk+1:]
                                if mot in n:
                                    n[mot] = n[mot] + 1
                                else:
                                    n[mot] = 1

    tacgt = acgt['A'] + acgt['C'] + acgt['G'] + acgt['T']
    ap = float(acgt['A'])/tacgt
    cp = float(acgt['C'])/tacgt
    gp = float(acgt['G'])/tacgt
    tp = float(acgt['T'])/tacgt

    fo.write("#Species\tNo. chr.\tGenome length\tA%\tC%\tG%\tT%\n")
    fo.write("#"+spec+"\t"+str(nchr)+"\t"+str(genome_len)+"\t"+str(ap)+"\t"+str(cp)+"\t"+str(gp)+"\t"+str(tp)+"\n")
    fo.write("#Motif\tObserved\tExpected\tScore\n")

    for motif in sorted(n):
        score, exp = calc_score(motif, ap, cp, gp, tp, genome_len, n[motif])
        fo.write(motif+"\t"+str(n[motif])+"\t"+str(exp)+"\t"+str(score)+"\n")

    fo.close()
    fi.close()

if __name__ == "__main__":
    main(sys.argv[1:])
