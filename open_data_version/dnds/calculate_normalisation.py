from Bio import SeqIO, Entrez
from Bio.Data import CodonTable

Entrez.email = "your.email@example.com"  # Replace with your email address

def all_possible_codons(codon):
    bases = ['A', 'T', 'C', 'G']
    for i in range(3):
        for base in bases:
            if base != codon[i]:
                yield codon[:i] + base + codon[i+1:]

def count_synonymous_nonsynonymous(codon_table, original, mutated):
    original_aa = codon_table[str(original)]
    mutated_aa = codon_table[str(mutated)]
    return (original_aa == mutated_aa, original_aa != mutated_aa)

def dnds_ratio(record):

    initial = dict(CodonTable.unambiguous_dna_by_id[1].forward_table)
    for s in CodonTable.unambiguous_dna_by_id[1].stop_codons:
      initial[s] = "*"
    
    
    synonymous_count = 0
    nonsynonymous_count = 0
    

    for feature in record.features:
        if feature.type == "CDS":
            if "ORF1a polyprotein" in feature.qualifiers.get('product',[]):
              print("ignoring orf1A duplication")
              continue
            #print(feature)
            cds = feature.extract(record.seq)
            for i in range(0, len(cds), 3):
                codon = cds[i:i+3]
                if len(codon) == 3:
                    for mutated_codon in all_possible_codons(codon):
                        syn, non_syn = count_synonymous_nonsynonymous(initial, codon, mutated_codon)
                        synonymous_count += syn
                        nonsynonymous_count += non_syn

    print(nonsynonymous_count,synonymous_count)
    return nonsynonymous_count , synonymous_count

def download_genbank(accession):
    with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
        return SeqIO.read(handle, "genbank")

accession = "NC_045512.2"
record = download_genbank(accession)
nonsynonymous_count , synonymous_count = dnds_ratio(record)
total_nucs_counted = (nonsynonymous_count+ synonymous_count)/3
total_nucs = len(record.seq)
diff = total_nucs - total_nucs_counted
print(f"We calculated that {diff} nucleotides were not included in the calculation - we will assume all muts at these sites are synonymous")
synonymous_count += diff*3 
ratio = nonsynonymous_count/synonymous_count

print("Ratio of mutations that give non-syn. changes to those that give syn. changes:", ratio)
