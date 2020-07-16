# covid_wgks_analysis
Repository for R scripts used for analysis of covid19 WGS

This repository contains several scripts, WGKS databases and test data to generate WGKS databases and also match WGKS to databases. In the following, several steps are described in the following:

I. Generation of WGKS
In order to generate a WGKS from a covid WGS, you can use either one of two Python k-mer analysis scripts. One, kmer_analysis.py anakyzes all k-mers, whereas kmer_analysis_mm.py analyzes all k-mers with up to [k/2] mutations. These scripts were written in Python3.

You can run the scripts in this manner:

A. k-mers without mismatches:
python3 kmer_analysis.py -i <input WGS> -o <output file name> -n <length of k-mer> -s <species name>
  
  Note: -n is an integer value, the species name must be a single string, i.e. "Homo_sapiens", and the input file must be in fasta format
  Example output:
  
  Sample output:
  
  head -6 ../../k6/MT012098.1.fasta.k6
#Species        No. chr.        Genome length   A%      C%      G%      T%
#SARS-CoV-2_human_IND_29_2020   1       29854   0.2986199504254036      0.18389495544985596     0.19628860454210492     0.32119648958263547
#Motif  Observed        Expected        Score
AAAAAA  2       21.169736025945483      -0.8273610024939085
AAAAAC  16      13.036663015416792      0.10205501172809861
AAAAAG  22      13.9152723625396        0.22510556389077938

Where the first 3 rows are annotation:
* line 1: Column names for row #2.
* line 2: name of species/WGS, no. chromosomes, genomelength, ACGT%
* line 3: column names for the rest of the file: motif (k-mer), observed occurrence, expected occurrence based on background ACGT%, and lastly, the score value.

B. k-mers with mismathes:
python3 kmer_analysis.py -i <input WGS> -o <output file name> -n <length of k-mer> -s <species name> -m <no. mismatches>

i.e. python3 kmer_analysis.py -i KY938558.1.fasta -o KY938558.1.fasta.k6m -n 6 -s Bat_coronavirus_strain_16BO133

Here -m is an integer value designating the number of mismctahes allowed in the k-mer, at most [k/2],so if k=5, at most 2 mismatches can be set.
The rest of the parameters are the same as in the first version, and the output is also the same, except that certain k-mers also appear with mismatches, i.e. AACNAT

II. Generationof WGKS database and matching WGKS for species identification

C. Generation of k-mer database:
A k-mer database is similar to a WGKS (either with or without mismatches), where the rows stanbd for k-mers, and the columns stand for WGS, i.e.:

        MN908947.3      MN938384.1      MN975262.1      MN985325.1
AAAAAA  0.163827258123403       -0.82714650627275       -0.0868566031942658     -0.406475179822576
AAAAAC  0.0951280852975193      0.102935323309245       0.0975457456821831      0.0992932867068953
AAAAAG  0.218384432437582       0.225512066892021       0.220620500632462       0.222298589267177
AAAAAN  0.0770444931525735      -0.116756375445301      0.0026720826109998      -0.0620983965378556
AAAAAT  -0.181473652401695      -0.174140264379581      -0.17925123440612       -0.177542760797033

This is the first six rows and first five columns of the file covid19_92_matrix.txt. For example, the k-mer AAAAAC has a score value of 0.102935... in WGS MN938384.1.

In order to generate a database, you need a number of WGKS (either with or without mismatches) in a directory. A test directory with 92 hexamer signarures has been provided in the file covid19.tar.gz. You need to gunzip and untar this dorectory into your working directory:

tar xvzf covid19.tar.gz

which will give you the directory with 92 .k6m files in it.

In order to generate the covid19 human database, you need to run the generateDB.R script sin the following way in R:

Rscript generateDB.R <directory with WGKS> <output file name>
  
  i.e. Rscript generateDB.R covid19 covid19_92_matrix.txt
  
  Run time for 92 WGKS is about 15 minutres.

D. Matching a WGKS against agiven WGKS database:
In order to match your own WGKS, generated by one of the two Python scripts in steps A or B, run the following command:

Rscript matchQuery.R <DB matrix file> <your WGKS file>
  
  i.e.
  
  Rscript matchQuery.R covid19_92_matrix.txt NC_045512.2.fasta.k6m
Reading WGKS DB...
Adding query to database...
Calculating correlation matrix...
[1] "P-value is 0.000761095776760526"
[1] "Mean DB PCC is 0.999214189705153"
[1] "Mean Query PCC is 0.99936955830903"

Here we ran a human covid sample, NC_045512.2.fasta.km6, which is the WGKS of the Wuhan-Hu-1 strain, with hexmares and 1 mismatch.

The output says that the p-value is 7.6e-4, and that the mean PCC values don't really differ.

Another example:

Rscript matchQuery.R batCovid19_48_matrix.txt KY938558.1.fasta.k6m
Reading WGKS DB...
Adding query to database...
Calculating correlation matrix...
[1] "P-value is 0.824113639614854"
[1] "Mean DB PCC is 0.842932155420701"
[1] "Mean Query PCC is 0.840870342390432"

Here KY938558.1.fasta.k6m is the hexamer WGKS for Bat_coronavirus_strain_16BO133, matched against the bat covid19 WGKS database batCovid19_48_matrix.txt. The p-value is 0.843, which means there is no statistical significance between the WGKS and the database, meaning that we are matching a bat WGKS with a bat database.
