import argparse
import gzip

def open_file(file_path, mode='r'):
    """Opens a file, compressed or not, based on its extension."""
    if file_path.endswith('.gz'):
        return gzip.open(file_path, mode + 't')
    else:
        return open(file_path, mode)


def analyze_output(file_path, metadata):
    """Analyze the output file to find cases where there are more than 15 mutations 
    and >30% of them are GtoA or >30% of them are CtoT."""
    results = []
    epi_isl_output_file = open('epi_isl_output.txt', 'w')
    with open_file(file_path, 'r') as file:
        header = file.readline().strip().split('\t')
        g_to_a_index = header.index("G>A")
        c_to_t_index = header.index("C>T")
        t_to_c_index = header.index("T>C")
        a_to_g_index = header.index("A>G")
        
        for line in file:
            data = line.strip().split('\t')
            seq_name = data[0]
            basic_mutations = [int(count) for count in data[-12:]]
            total_mutations = sum(basic_mutations)
            g_to_a_count = int(data[g_to_a_index])
            c_to_t_count = int(data[c_to_t_index])
            t_to_c_count = int(data[t_to_c_index])
            a_to_g_count = int(data[a_to_g_index])
            transition_count = g_to_a_count + c_to_t_count + t_to_c_count + a_to_g_count
            
            if  ( total_mutations > 20 and (g_to_a_count / total_mutations > 0.25) and (c_to_t_count / total_mutations > 0.25) and (transition_count / total_mutations > 0.9)) or ( total_mutations > 30 and (g_to_a_count / total_mutations > 0.2) and (c_to_t_count / total_mutations > 0.2) and (transition_count / total_mutations > 0.85)) or ( total_mutations > 13 and (g_to_a_count / total_mutations > 0.3) and (c_to_t_count / total_mutations > 0.3) and (transition_count / total_mutations > 0.95) )or ( total_mutations > 10 and (g_to_a_count / total_mutations > 0.4) and (c_to_t_count / total_mutations > 0.4) and (transition_count / total_mutations > 0.99) ):
                    results.append(seq_name)
                    row = metadata.loc[seq_name]
                    # include date and country
                    print(seq_name + '\t' + str(total_mutations) + '\t' + str(g_to_a_count) + '\t' + str(c_to_t_count) + '\t' + row['date'] + '\t' + row['country'] + '\t' + row['gisaid_epi_isl'])
                    epi_isl_output_file.write(row['gisaid_epi_isl'] + '\n')

    
    return results
import pandas as pd
def main():
    parser = argparse.ArgumentParser(description="Analyze the output file for specific mutation patterns.")
    parser.add_argument("--file_path", help="Path to the output file. Supports gzipped files.", required=True)
    parser.add_argument("--metadata", help="Path to the metadata file. Supports gzipped files.", required=True)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep='\t', usecols=['strain', 'date', 'country', 'gisaid_epi_isl'])
    metadata = metadata.set_index('strain')
    
    cases_of_interest = analyze_output(args.file_path, metadata)
 
if __name__ == "__main__":
    main()
