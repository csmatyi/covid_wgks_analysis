#!/usr/bin/python

import re, os, sys, getopt
import itertools
from itertools import combinations

def calc_score_2(motif, n, glen):
    mlen = len(motif)
    mlen1 = mlen - 1
    mlen2 = mlen - 2
    m1a = motif[0:mlen1]
    m1b = motif[1:mlen]
    m2 = motif[1:mlen1]
    q = 0
    r = 0
    s = 1
    if m1a in n.keys():
        q = float(n[m1a])
    if m1b in n.keys():
        r = float(n[m1b])
    if m2 in n.keys():
        s = float(n[m2])
    #print(motif+"     "+m1a+" "+str(q)+" "+m1b+" "+str(r)+" "+m2+" "+str(s))
    exp = q * r / s
    score = float((n[motif] - exp)/(n[motif] + exp))
    return score, exp

def mm_motif_list(motif, m):
    ms = list()
    ms.append(motif)
    for mm in range(1, m + 1):
        c = combinations(range(len(motif)),mm)
        for jj in c:
            mot = motif
            for kk in jj:
                mot = mot[:kk]+'N'+mot[kk+1:]
                ms.append(mot)
    return list(set(ms))

def main(argv):
    inputfile = ''
    outputfile = ''
    motif = ''
    pos = 0
    k_ = 0
    ns = 0
    n = dict()
    n1 = dict()
    n2 = dict()
    bplist = list()
    acgt = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'M': 0, 'R': 0, 'W': 0, 'S': 0, 'K': 0, 'Y': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0}
    genome_len = 0
    nchr = 0
    spec = ''
    m = 0

    try:
        opts, args = getopt.getopt(argv,"hi:o:s:k:m:",["ifile=","ofile=","spec=","k=","mismatches="])
    except getopt.GetoptError:
        print ('kmer_analysis_km.py -i <inputfile> -o <outputfile> -s <species name> -k <kmer length> -m <mismatches>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('kmer_analysis_km.py -i <inputfile> -o <outputfile> -s <species name> -k <kmer length> -m <mismatches>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-k"):
            k_ = int(arg)
        elif opt in ("-m"):
            m = int(arg)
            if m/k_ > 0.5:
                print ('Mismatches cannot be longer than half of the input k-mer length!')
                sys.exit()
        elif opt in ("-s"):
            spec = arg

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
                pos = pos + 1
                genome_len = genome_len + 1

                bplist.append(b)
                if pos >= k_ + 1:
                    bplist.pop(0)

                if pos >= k_:
                    motif = ''.join(bplist).upper()
                    motifs = mm_motif_list(motif, m) #
                    for mot in motifs:
                        x1 = k_ - 1
                        x2 = k_ - 2
                        motif1 = mot[0:x1]
                        motif2 = mot[0:x2]
                        n[mot] = n[mot] + 1 if mot in n else 1
                        n[motif1] = n[motif1] + 1 if motif1 in n else 1
                        n[motif2] = n[motif2] + 1 if motif2 in n else 1

    tacgt = acgt['A'] + acgt['C'] + acgt['G'] + acgt['T']
    ns = acgt['B'] + acgt['D'] + acgt['H'] + acgt['V'] + acgt['M'] + acgt['R'] + acgt['W'] + acgt['S'] + acgt['Y'] + acgt['K'] + acgt['N']
    np = float(ns)/genome_len
    ap = float(acgt['A'])/tacgt
    cp = float(acgt['C'])/tacgt
    gp = float(acgt['G'])/tacgt
    tp = float(acgt['T'])/tacgt

    fo.write("#Species\tNo. chr.\tGenome length\tA%\tC%\tG%\tT%\tN%\n")
    fo.write("#"+spec+"\t"+str(nchr)+"\t"+str(genome_len)+"\t"+str(ap)+"\t"+str(cp)+"\t"+str(gp)+"\t"+str(tp)+"\t"+str(np)+"\n")
    fo.write("#Motif\tObserved\tExpected\tScore\n")

    for motif in sorted(n):
        if len(motif) == k_:
            score, exp = calc_score_2(motif, n, genome_len)
            fo.write(motif+"\t"+str(n[motif])+"\t"+str(exp)+"\t"+str(score)+"\n")

    fo.close()
    fi.close()

if __name__ == "__main__":
    main(sys.argv[1:])
