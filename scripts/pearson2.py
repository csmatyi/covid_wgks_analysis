#!/usr/bin/python

import re, os, sys, getopt
import scipy
from scipy import stats

def file2dict(infile):
    fi = open(infile,'r')
    c = 0
    d = {}
    for line in fi:
        c = c + 1
        if c > 3:
            lin = line.rstrip('\n')
            motif, exp, occ, score = lin.split('\t')
            d[motif] = score
    fi.close()
    return d

def intersect(sc1, sc2):
    sc1b = []
    sc2b = []
    for k1 in sc1:
        if k1 in sc2.keys():
            sc1b.append(float(sc1[k1]))
            sc2b.append(float(sc2[k1]))
    return sc1b, sc2b

def scoredist(sc1, sc2):
    d = 0
    c = 0
    for k1 in sc1:
        if k1 in sc2.keys():
            c = c + 1
            d = d + abs(float(sc1[k1]) - float(sc2[k1]))
    dist = float(d/c)
    return dist

def main(argv):
    inputfile1 = ''
    inputfile2 = ''
    outputfile = ''

    try:
        opts, args = getopt.getopt(argv,"hi:j:o:",["ifile=","jfile=","ofile="])
    except getopt.GetoptError:
        print('person.py -i <inputfile 1> ij <inputfile 2> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('pearson.py -i <inputfile 1> -j <inoutfile 2> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile1 = arg
        elif opt in ("-j", "--jfile"):
            inputfile2 = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    scores1 = file2dict(inputfile1)
    scores2 = file2dict(inputfile2)
    (scores1b, scores2b) = intersect(scores1, scores2)
    (r, p) = scipy.stats.pearsonr(scores1b, scores2b)
    (w, pw) = scipy.stats.wilcoxon(scores1b, scores2b)
    (x, px) = scipy.stats.chisquare(scores1b, scores2b)
    d = scoredist(scores1, scores2)
    fo = open(outputfile,'w')
    fo.write(inputfile1+"\t"+inputfile2+"\t"+str(r)+"\t"+str(p)+"\t"+str(w)+"\t"+str(pw)+"\t"+str(x)+"\t"+str(px)+"\t"+str(d)+"\n")
    fo.close()

if __name__ == "__main__":
    main(sys.argv[1:])
