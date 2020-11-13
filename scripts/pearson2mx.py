#!/usr/bin/python

import re, os, sys, getopt, math

def getsg(sgfile):
    sg = dict()
    fsg = open(sgfile,'r')
    for line in fsg:
        lin = line.rstrip('\n')
        (s, g) = lin.split()
        #sg[g] = s
        sg[s] = g
    return sg

def getmx(ifile, sg):
    corr = dict()
    species = []
    fi = open(ifile,'r')
    for line in fi:
        lin = line.rstrip('\n')
        (sp1, sp2, cc, pcc, wcx, pwcx, x2, px2, d) = lin.split('\t')
        spa = sg[sp1]
        spb = sg[sp2]
        corr[(spa, spb)] = cc
        corr[(spb, spa)] = cc
        species.append(spa)

    species.sort()
    spec = list(set(species))
    spec.sort()

    return corr, spec

def main(argv):
    sgfile = ''
    inputfile = ''
    mxfile = ''

    try:
        opts, args = getopt.getopt(argv,"hs:i:m:",["sfile=","ifile=","mfile="])
    except getopt.GetoptError:
        print ('pearson2mx.py -s <speciesfile> -i <inputfile> -m <mxfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('pearson2mx.py -s <speciesfile> -i <inputfile> -m <mxfile>')
            sys.exit()
        elif opt in ("-s", "--sfile"):
            sgfile = arg
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-m", "--mfile"):
            mxfile = arg

    sg = getsg(sgfile)
    rmx, species = getmx(inputfile, sg)

    fo = open(mxfile,'w')
    # header
    fo.write("ID")
    for s in species:
        fo.write("\t"+s)
    fo.write('\n')
    # row by row
    for s2 in species:
        fo.write(s2)
        for s3 in species:
            r = str(rmx[(s2,s3)])
            fo.write("\t"+r)
        fo.write('\n')
    fo.close()

if __name__ == "__main__":
    main(sys.argv[1:])

