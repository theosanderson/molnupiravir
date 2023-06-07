#Extracts mutations from the second trial samples
#Identifies mutations as those present at >=5% at positions with >=100 reads

import argparse
from collections import OrderedDict

#Creates an empty mutational spectrum dictionary for RNA datasets, i.e. not combining symmetric mutations
#Copied from reconstruct_spectrum.py in MutTui
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

#Generates reference dictionary from a variant file
#Sites as keys, bases as values
def getRef(vF):
    ref = dict()

    with open(vF) as f:
        for l in f:
            ref[int(l.strip().split("\t")[0])] = l.strip().split("\t")[1]

    return(ref)

parser = argparse.ArgumentParser()
parser.add_argument("-v", nargs = "+", help = "Variant files")
args = parser.parse_args()

#Mutations identified in each patient
pM = dict()

#Nucleotides and positions in files
nucPos = {"A":2, "C":3, "G":4, "T":5}
nucleotides = ["A", "C", "G", "T"]

#Mutation lists
mM = list()
nM = list()
pM = list()

#Spectra
mSpectrum = getRNADict()
nSpectrum = getRNADict()
pSpectrum = getRNADict()

#Output files
outM = open("MOV_mutations.csv", "w")
outM.write("Patient,Position,Reference,Variant,Mutation\n")
outN = open("naive_mutations.csv", "w")
outN.write("Patient,Position,Reference,Variant,Mutation\n")
outP = open("paxlovid_mutations.csv", "w")
outP.write("Patient,Position,Reference,Variant,Mutation\n")

outMSpectrum = open("MOV_spectrum.csv", "w")
outMSpectrum.write("Substitution,Number_of_mutations\n")
outNSpectrum = open("naive_spectrum.csv", "w")
outNSpectrum.write("Substitution,Number_of_mutations\n")
outPSpectrum = open("paxlovid_spectrum.csv", "w")
outPSpectrum.write("Substitution,Number_of_mutations\n")

#Iterate through the variant files, extract variants
for vF in args.v:
    #Drug treatment and patient
    fn = vF.split("/")[-1]
    d = fn.split("_")[0]
    patient = fn.split("_")[0] + "_" + fn.split("_")[1]

    #Extract reference
    ref = getRef(vF)
    
    with open(vF) as f:
        for l in f:
            #Sum reads at position
            ls = l.strip().split("\t")
            nr = int(ls[2]) + int(ls[3]) + int(ls[4]) + int(ls[5])

            #Only analyse sites with >=100 reads and where the consensus is a nucleotide
            if (nr >= 100) and (ls[1] != "N"):
                #Calculate proportion of non-reference bases
                for eN in nucPos:
                    if eN != ls[1]:
                        #Include variants at >=5%
                        if int(ls[nucPos[eN]])/nr >= 0.05:
                            #Extract context
                            s = int(ls[0])
                            if ((s - 1) in ref) and ((s + 1) in ref):
                                context = ref[s-1], ref[s+1]

                                if (context[0] in nucleotides) and (context[1] in nucleotides):
                                    #Add mutation to respective list
                                    if d == "MOLNUPIRAVIR":
                                        mM.append(patient + "__" + ls[0] + "__" + ls[1] + "__" + eN + "__" + context[0] + ls[1] + eN + context[1])
                                    elif d == "NAIVE":
                                        nM.append(patient + "__" + ls[0] + "__" + ls[1] + "__" + eN + "__" + context[0] + ls[1] + eN + context[1])
                                    elif d == "PAXLOVID":
                                        pM.append(patient + "__" + ls[0] + "__" + ls[1] + "__" + eN + "__" + context[0] + ls[1] + eN + context[1])

#Write mutations and add to spectra
for eM in set(mM):
    outM.write(",".join(eM.split("__")) + "\n")
    mSpectrum[eM.split("__")[-1]] += 1
for eM in set(nM):
    outN.write(",".join(eM.split("__")) + "\n")
    nSpectrum[eM.split("__")[-1]] += 1
for eM in set(pM):
    outP.write(",".join(eM.split("__")) + "\n")
    pSpectrum[eM.split("__")[-1]] += 1

#Write spectra
for eM in mSpectrum:
    outMSpectrum.write(eM[0] + "[" + eM[1] + ">" + eM[2] + "]" + eM[3] + "," + str(mSpectrum[eM]) + "\n")
    outNSpectrum.write(eM[0] + "[" + eM[1] + ">" + eM[2] + "]" + eM[3] + "," + str(nSpectrum[eM]) + "\n")
    outPSpectrum.write(eM[0] + "[" + eM[1] + ">" + eM[2] + "]" + eM[3] + "," + str(pSpectrum[eM]) + "\n")

outM.close()
outN.close()
outP.close()
outMSpectrum.close()
outNSpectrum.close()
outPSpectrum.close()