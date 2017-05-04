#!/usr/bin/python



"""This program takes as input a tab-delimted file (below named pairwise.txt) that is output from Singletons_to_pairwise.py, 
and has the following structure:

GeneID1 species1 GeneID2 species2

This program calculates the dN, dS, and dN/dS values for each gene pair using codon-aligned CDS sequences. To do this, there must 
be a file of all proteins seqs (below names PEP.fasta) and CDS seqs (below named CDS.fasta). The files produces by this program are:

alignments.csv (contains each pair of genes, and their aligned protein and CDS sequences)
final_dNdS.csv (contains each pair of genes, and their dN, dS, and dN/dS values)

Usage: .py pairwise.txt PEP.fasta CDS.fasta"""


import time
import sys
import numpy as np
import pandas as pd
import re
from Bio.Phylo.PAML import codeml
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import os

start_time = time.time()
######################################################################################
###########################    START FUNCTIONS   #####################################
######################################################################################


def gapsFromPeptide( peptide_seq, nucleotide_seq ):
    """ Transfers gaps from aligned peptide seq into codon partitioned nucleotide seq (codon alignment)
          - peptide_seq is an aligned peptide sequence with gaps that need to be transferred to nucleotide seq
          - nucleotide_seq is an un-aligned dna sequence whose codons translate to peptide seq"""
    def chunks(l, n):
        """ Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield l[i:i+n]
    codons = [codon for codon in chunks(nucleotide_seq,3)]  #splits nucleotides into codons (triplets)
    gappedCodons = []
    codonCount = 0
    for aa in peptide_seq:  #adds '---' gaps to nucleotide seq corresponding to peptide
        if aa!='-':
            gappedCodons.append(codons[codonCount])
            codonCount += 1
        else:
            gappedCodons.append('---')
    return(''.join(gappedCodons))

def run_codeml_on_everything(CDS_a, CDS_b):
    """Takes a dictionary of gene NAMES and prints the aligned CDS sequences in Phylip format"""
    assert len(CDS_a) == len(CDS_b)
    with open("alignment.phylip", "w") as handle:
        handle.write("%s%i\n%s\n%s\n%s\n%s" % ("2 ", len(CDS_a), "query", CDS_a, "subject", CDS_b))
    with open("treefile.txt", "w") as handle:
        handle.write("(query, subject)")
    cml = codeml.Codeml()
    cml.tree = "treefile.txt"
    cml.alignment = "alignment.phylip"
    cml.working_dir = "./work"
    cml.out_file = "codeml_outfile.txt"
    #cml.Small_Diff = .5e-6 #
    cml.set_options(RateAncestor = 1, #
            Small_Diff = .5e-6, #
            #aaRatefile = 0, #
            #icode = 0, #
            model = 0, #
            noisy = 0, #
            runmode = -2, #
            #Mgene = 0, #
            #method = 0 #
            #Malpha = 0, #
            seqtype = 1, #
            #omega = 0, #
            #getSE = 0, #
            verbose = True, #
            CodonFreq = 2, #
            #cleandata = 0, #
            #fix_blength = 0, #
            #NSsites = 0, #
            #fix_omega = 0, #
            #clock = 0, #
            kappa = 2, #
            ndata = 1, #
            #aaDist = 0, #
            fix_alpha = 1, #
            #alpha = 0, #
            #fix_kappa = 0, #
            ncatG = 1) #
    results = cml.run()
    path = results['pairwise']['query']['subject']
    return path['dN'], path['dS'], path['omega']

def prot_SeqIO_to_pandas(IDs):
    return prot_seqs[IDs].seq[:-1]

def CDS_SeqIO_to_pandas(IDs):
    return CDS_seqs[IDs].seq[:-3]

def pairwise_align(first, second):
    matrix = matlist.blosum62
    answer = pairwise2.align.globaldx(first, second, matrix)
    return answer[0][0], answer[0][1]


######################################################################################
###########################    END OF FUNCTIONS    ###################################
######################################################################################

start_time = time.time()

"""This part of the script filters out reciprocal hits and self-hits"""

prot_seqs = SeqIO.index(sys.argv[2], "fasta")
CDS_seqs = SeqIO.index(sys.argv[3], "fasta")

end_time = time.time()
total_time = end_time - start_time
print("Time to index fasta files:")
print(total_time)
start_time = time.time()

name = ["qseqid", "qspecies", "sseqid", "sspecies"]
df = pd.read_table(sys.argv[1], names = name) 


df = df[(df.qseqid != df.sseqid)] #eliminate self-hits
idx = (df['qseqid'] > df['sseqid'])
df.loc[idx,['qseqid', 'qspecies', 'sseqid', 'sspecies']] = df.loc[idx,['sseqid', 'sspecies', 'qseqid', 'qspecies']].values
df = df[df.duplicated()]

df.to_csv("doesitwork.csv")

end_time = time.time()
total_time = end_time - start_time
print("Time to filter blast file:")
print(total_time)
start_time = time.time()
print("The number of gene pairs identified is:")
print(df.shape[0])








"""THIS PART OF THE SCRIPT PREPARES ALIGNS PROTEIN SEQS AND CDS SEQS"""


df['qprot'] = df.apply(lambda row: prot_SeqIO_to_pandas(row['qseqid']), axis = 1)
df['sprot'] = df.apply(lambda row: prot_SeqIO_to_pandas(row['sseqid']), axis = 1)
df['qprot'] = df['qprot'].astype(str)
df['sprot'] = df['sprot'].astype(str)

end_time = time.time()
total_time = end_time - start_time
print("Time to pull out prot seqs:")
print(total_time)
start_time = time.time()



df[['qprotalign','sprotalign']] = df.apply(lambda row : pd.Series(pairwise_align(row['qprot'], row['sprot'])), axis = 1)
df.drop('qprot', axis = 1, inplace = True)
df.drop('sprot', axis = 1, inplace = True)

end_time = time.time()
total_time = end_time - start_time
print("Time to align prot seqs:")
print(total_time)
start_time = time.time()

df['qCDS'] = df.apply(lambda row: CDS_SeqIO_to_pandas(row['qseqid']), axis = 1)
df['sCDS'] = df.apply(lambda row: CDS_SeqIO_to_pandas(row['sseqid']), axis = 1)
df['qCDS'] = df["qCDS"].astype(str)
df['sCDS'] = df["sCDS"].astype(str)

#df['qprotlen'] = df['qprot'].str.len()
#df['qCDSlen'] = df['qCDS'].str.len()
#df['sprotlen'] = df['sprot'].str.len()
#df['sCDSlen'] = df['sCDS'].str.len()

end_time = time.time()
total_time = end_time - start_time
print("Time to pull out CDS seqs:")
print(total_time)
start_time = time.time()


df['qCDSalign'] = df.apply(lambda row: gapsFromPeptide(row['qprotalign'], row['qCDS']), axis = 1) #repeatedly fails. Dont know what I changed. 
df['sCDSalign'] = df.apply(lambda row: gapsFromPeptide(row['sprotalign'], row['sCDS']), axis = 1)
df.drop('qCDS', axis = 1, inplace = True)
df.drop('sCDS', axis = 1, inplace = True)


df.to_csv("alignments.csv")

df.drop('qprotalign', axis = 1, inplace=True)
df.drop('sprotalign', axis = 1, inplace = True)

end_time = time.time()
total_time = end_time - start_time
print("Time to align CDS:")
print(total_time)
start_time = time.time()




"""THIS PART OF THE SCRIPT 
1. TAKES THE ALIGNED CDS FILES AND FORMATS THEM FOR CODEML 
2. RUNS CODEML
3. AND RETURNS RESULTS INTO SAME DATAFRAME """

with open("treefile.txt", "w") as handle:
    handle.writelines("(query, subject)")

df[['dN', 'dS', 'dN/dS']] = df.apply(lambda row: pd.Series(run_codeml_on_everything(row['qCDSalign'], row['sCDSalign'])), axis = 1)




df.drop('qCDSalign', axis = 1, inplace = True)
df.drop('sCDSalign', axis = 1, inplace = True)





end_time = time.time()
total_time = end_time - start_time
print("Time to run Codeml:")
print(total_time)

#print(df)
df.to_csv("final_dNdS.csv")  
#print(df.shape[0])

