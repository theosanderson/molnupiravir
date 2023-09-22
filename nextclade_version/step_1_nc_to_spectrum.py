import argparse
import gzip

def open_file(file_path, mode='r'):
    """Opens a file, compressed or not, based on its extension."""
    if file_path.endswith('.gz'):
        return gzip.open(file_path, mode + 't')
    else:
        return open(file_path, mode)

def parse_fasta(ref_file):
    """Parses a FASTA file and returns the sequence as a string."""
    sequence = "-" # to make it one-indexed
    with open_file(ref_file, 'r') as file:
        for line in file:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence


def extract_contextual_mutations(nextclade_file, ref_file, output_file):
    reference_sequence = parse_fasta(ref_file)
    

    
    with open_file(nextclade_file, 'r') as file:
        header = file.readline().strip().split('\t')
        if 'privateNucMutations.unlabeledSubstitutions' not in header or "substitutions" not in header:
            raise ValueError("The required columns are not present in the tsv file.")
        private_mutations_index = header.index('privateNucMutations.unlabeledSubstitutions')
        substitutions_index = header.index('substitutions')
        seqname_index = header.index('seqName')
        all_possible_mutation_contexts = []
        all_bases = ['A', 'C', 'G', 'T']
        for base1 in all_bases:
            for base2 in all_bases:
                for base3 in all_bases:
                    for base4 in all_bases:
                        if base3!=base4:
                            all_possible_mutation_contexts.append(base1+'['+base3+'>'+base4+']'+base2)
        all_simple_mutations = []
        for base1 in all_bases:
            for base2 in all_bases:
                if base1!=base2:
                    all_simple_mutations.append(base1+'>'+base2)
        print("seqName\t" + "\t".join(all_possible_mutation_contexts) + "\t" + "\t".join(all_simple_mutations) )
        
        for line in file:
            line_data = line.strip().split('\t')
            seq_name = line_data[seqname_index]
            substitutions = line_data[substitutions_index].strip().split(',')
            substitution_dict = {}
            for substitution in substitutions:
                if substitution.strip() == '':
                    continue
                position = int(substitution[1:-1]) - 1
                final_nucleotide = substitution[-1]
                substitution_dict[position] = final_nucleotide
                
            
            private_mutations = line_data[private_mutations_index].strip().split(',')
            basic_mutation_counts = {}
            mutation_counts = {}
            for mutation in private_mutations:
                mutation = mutation.strip()
                if mutation == '':
                    continue
                try:
                    position = int(mutation[1:-1]) - 1  # Convert to 0-indexed position
                    if 0 < position < len(reference_sequence) - 1:  # Check if position is valid
                        context = (reference_sequence[position - 1] if not (position-1) in substitution_dict else substitution_dict[(position-1)]) + '[' + mutation[0] + '>' + mutation[-1] + ']' + (reference_sequence[position + 1] if not (position+1) in substitution_dict else substitution_dict[(position+1)])
                        simple_mutation = mutation[0] + '>' + mutation[-1]
                        mutation_counts[context] = mutation_counts.get(context, 0) + 1
                        basic_mutation_counts[simple_mutation] = basic_mutation_counts.get(simple_mutation, 0) + 1
                except Exception as e:
                    print("AA",line_data[private_mutations_index],'AA')
                    print("Error processing mutation: " + private_mutations)
                    raise e
        
            print(seq_name + "\t" + "\t".join([str(mutation_counts.get(context, 0)) for context in all_possible_mutation_contexts]) + "\t" + "\t".join([str(basic_mutation_counts.get(context, 0)) for context in all_simple_mutations]) )

def main():
    parser = argparse.ArgumentParser(description="Extract and count 192 contextual mutations from nextclade.tsv using a reference FASTA.")
    parser.add_argument("--nextclade_file", help="Path to the nextclade.tsv file. Supports gzipped files.")
    parser.add_argument("--ref_file", help="Path to the reference FASTA file. Supports gzipped files.")
    parser.add_argument("--output_file", help="Path to the output file where mutation counts will be saved. Supports gzipped files.")
    
    args = parser.parse_args()
    
    extract_contextual_mutations(args.nextclade_file, args.ref_file, args.output_file)
    print(f"Mutation counts written to {args.output_file}.")

if __name__ == "__main__":
    main()

    """ python nc_to_spectrum.py --ref_file ~/ref.fa.fasta --nextclade_file ~/nc.tsv.gz --output_file out.tsv.gz"""