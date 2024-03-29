#Calculates a spectrum of long MOV branches from the UShER tree

import argparse
from collections import OrderedDict
from Bio import SeqIO
import numpy as np

#Creates an empty mutational spectrum dictionary for RNA datasets, i.e. not combining symmetric mutations
def getRNADict():
    mutation = OrderedDict()
    mutation["ACAA"] = 0
    mutation["ACAC"] = 0
    mutation["ACAG"] = 0
    mutation["ACAT"] = 0
    mutation["CCAA"] = 0
    mutation["CCAC"] = 0
    mutation["CCAG"] = 0
    mutation["CCAT"] = 0
    mutation["GCAA"] = 0
    mutation["GCAC"] = 0
    mutation["GCAG"] = 0
    mutation["GCAT"] = 0
    mutation["TCAA"] = 0
    mutation["TCAC"] = 0
    mutation["TCAG"] = 0
    mutation["TCAT"] = 0
    mutation["ACGA"] = 0
    mutation["ACGC"] = 0
    mutation["ACGG"] = 0
    mutation["ACGT"] = 0
    mutation["CCGA"] = 0
    mutation["CCGC"] = 0
    mutation["CCGG"] = 0
    mutation["CCGT"] = 0
    mutation["GCGA"] = 0
    mutation["GCGC"] = 0
    mutation["GCGG"] = 0
    mutation["GCGT"] = 0
    mutation["TCGA"] = 0
    mutation["TCGC"] = 0
    mutation["TCGG"] = 0
    mutation["TCGT"] = 0
    mutation["ACTA"] = 0
    mutation["ACTC"] = 0
    mutation["ACTG"] = 0
    mutation["ACTT"] = 0
    mutation["CCTA"] = 0
    mutation["CCTC"] = 0
    mutation["CCTG"] = 0
    mutation["CCTT"] = 0
    mutation["GCTA"] = 0
    mutation["GCTC"] = 0
    mutation["GCTG"] = 0
    mutation["GCTT"] = 0
    mutation["TCTA"] = 0
    mutation["TCTC"] = 0
    mutation["TCTG"] = 0
    mutation["TCTT"] = 0
    mutation["ATAA"] = 0
    mutation["ATAC"] = 0
    mutation["ATAG"] = 0
    mutation["ATAT"] = 0
    mutation["CTAA"] = 0
    mutation["CTAC"] = 0
    mutation["CTAG"] = 0
    mutation["CTAT"] = 0
    mutation["GTAA"] = 0
    mutation["GTAC"] = 0
    mutation["GTAG"] = 0
    mutation["GTAT"] = 0
    mutation["TTAA"] = 0
    mutation["TTAC"] = 0
    mutation["TTAG"] = 0
    mutation["TTAT"] = 0
    mutation["ATCA"] = 0
    mutation["ATCC"] = 0
    mutation["ATCG"] = 0
    mutation["ATCT"] = 0
    mutation["CTCA"] = 0
    mutation["CTCC"] = 0
    mutation["CTCG"] = 0
    mutation["CTCT"] = 0
    mutation["GTCA"] = 0
    mutation["GTCC"] = 0
    mutation["GTCG"] = 0
    mutation["GTCT"] = 0
    mutation["TTCA"] = 0
    mutation["TTCC"] = 0
    mutation["TTCG"] = 0
    mutation["TTCT"] = 0
    mutation["ATGA"] = 0
    mutation["ATGC"] = 0
    mutation["ATGG"] = 0
    mutation["ATGT"] = 0
    mutation["CTGA"] = 0
    mutation["CTGC"] = 0
    mutation["CTGG"] = 0
    mutation["CTGT"] = 0
    mutation["GTGA"] = 0
    mutation["GTGC"] = 0
    mutation["GTGG"] = 0
    mutation["GTGT"] = 0
    mutation["TTGA"] = 0
    mutation["TTGC"] = 0
    mutation["TTGG"] = 0
    mutation["TTGT"] = 0
    mutation["AGTA"] = 0
    mutation["AGTC"] = 0
    mutation["AGTG"] = 0
    mutation["AGTT"] = 0
    mutation["CGTA"] = 0
    mutation["CGTC"] = 0
    mutation["CGTG"] = 0
    mutation["CGTT"] = 0
    mutation["GGTA"] = 0
    mutation["GGTC"] = 0
    mutation["GGTG"] = 0
    mutation["GGTT"] = 0
    mutation["TGTA"] = 0
    mutation["TGTC"] = 0
    mutation["TGTG"] = 0
    mutation["TGTT"] = 0
    mutation["AGCA"] = 0
    mutation["AGCC"] = 0
    mutation["AGCG"] = 0
    mutation["AGCT"] = 0
    mutation["CGCA"] = 0
    mutation["CGCC"] = 0
    mutation["CGCG"] = 0
    mutation["CGCT"] = 0
    mutation["GGCA"] = 0
    mutation["GGCC"] = 0
    mutation["GGCG"] = 0
    mutation["GGCT"] = 0
    mutation["TGCA"] = 0
    mutation["TGCC"] = 0
    mutation["TGCG"] = 0
    mutation["TGCT"] = 0
    mutation["AGAA"] = 0
    mutation["AGAC"] = 0
    mutation["AGAG"] = 0
    mutation["AGAT"] = 0
    mutation["CGAA"] = 0
    mutation["CGAC"] = 0
    mutation["CGAG"] = 0
    mutation["CGAT"] = 0
    mutation["GGAA"] = 0
    mutation["GGAC"] = 0
    mutation["GGAG"] = 0
    mutation["GGAT"] = 0
    mutation["TGAA"] = 0
    mutation["TGAC"] = 0
    mutation["TGAG"] = 0
    mutation["TGAT"] = 0
    mutation["AATA"] = 0
    mutation["AATC"] = 0
    mutation["AATG"] = 0
    mutation["AATT"] = 0
    mutation["CATA"] = 0
    mutation["CATC"] = 0
    mutation["CATG"] = 0
    mutation["CATT"] = 0
    mutation["GATA"] = 0
    mutation["GATC"] = 0
    mutation["GATG"] = 0
    mutation["GATT"] = 0
    mutation["TATA"] = 0
    mutation["TATC"] = 0
    mutation["TATG"] = 0
    mutation["TATT"] = 0
    mutation["AAGA"] = 0
    mutation["AAGC"] = 0
    mutation["AAGG"] = 0
    mutation["AAGT"] = 0
    mutation["CAGA"] = 0
    mutation["CAGC"] = 0
    mutation["CAGG"] = 0
    mutation["CAGT"] = 0
    mutation["GAGA"] = 0
    mutation["GAGC"] = 0
    mutation["GAGG"] = 0
    mutation["GAGT"] = 0
    mutation["TAGA"] = 0
    mutation["TAGC"] = 0
    mutation["TAGG"] = 0
    mutation["TAGT"] = 0
    mutation["AACA"] = 0
    mutation["AACC"] = 0
    mutation["AACG"] = 0
    mutation["AACT"] = 0
    mutation["CACA"] = 0
    mutation["CACC"] = 0
    mutation["CACG"] = 0
    mutation["CACT"] = 0
    mutation["GACA"] = 0
    mutation["GACC"] = 0
    mutation["GACG"] = 0
    mutation["GACT"] = 0
    mutation["TACA"] = 0
    mutation["TACC"] = 0
    mutation["TACG"] = 0
    mutation["TACT"] = 0

    return(mutation)

#Updates contexts leading to a given mutation if there are previous surrounding mutations
def updateContexts(ref, l, m):
    p = int(m[1:-1])

    #Upstream and downstream base
    u = ref[p - 2]
    d = ref[p]

    #Update contexts if needed
    for node in l.strip().split("\t")[1].split(m)[0].split(" ")[:-1]:
        for eM in node.split(":")[1].split(","):
            if int(eM[1:-1]) == (p - 1):
                u = eM[-1]
            elif int(eM[1:-1]) == (p + 1):
                d = eM[-1]

    return(u, d)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", help = "sample-paths file")
    parser.add_argument("-r", help = "Wuhan-Hu1 reference used to get contexts")
    parser.add_argument("-o", help = "Name of output file")
    args = parser.parse_args()

    #Empty spectrum
    spectrum = getRNADict()

    #Import and extract reference to an array
    for record in SeqIO.parse(args.r, "fasta"):
        ref = np.array(record.seq)
    
    #Analysed nodes
    nodes = set()

    #Iterate through the sample-paths, identify MOV branches and extract their contextual mutations
    with open(args.s) as f:
        for l in f:
            for n in l.strip().split("\t")[1].split(" "):
                #Check if the branch has at least 10 mutations
                if (n.count(",") >= 9) and (n.split(":")[0] not in nodes):
                    for eM in n.split(":")[1].split(","):
                        #Update surrounding context if necessary
                        u, d = updateContexts(ref, l, eM)
                        spectrum[u + eM[0] + eM[-1] + d] += 1
                        nodes.add(n.split(":")[0])
    
    #Write spectrum
    out = open(args.o, "w")
    out.write("Substitution,Number_of_mutations\n")

    for m in spectrum:
        out.write(m[0] + "[" + m[1] + ">" + m[2] + "]" + m[3] + "," + str(spectrum[m]) + "\n")
    
    out.close()