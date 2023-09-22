# Alternative reproduction with Nextclade

Our primary analysis is based on a mutation-annotated tree built with all data from GISAID and the INSDC. Due to GISAID policies we are not allowed to redistribute this file online. As an alternative approach to allow anyone with access to GISAID "Genomic Epidemiology" packages to reproduce our core results, we provide a version that uses Nextclade here. Nextclade calls of private mutations for each sequence allow the identification of molnupiravir-associated sequences.

(In another folder you can see a version of our analysis using just open data sequences from the UCSC public tree.)

We would also recommend a Nextclade-based approach for those looking to prospectively identify molnupiravir-based sequences in the future.

## Step 0:
Run [Nextclade](https://clades.nextstrain.org/) CLI on the sequences file, and output a nextclade.tsv file.

## Step 1:
Run step_1_nc_to_spectrum.py on the nextclade.tsv to identify the types of mutations present in the `privateNucMutations.unlabeledSubstitutions` of each sequence. Provide the Hu1 spectrum to allow calculation of contexts.

## Step 2:
Run step_2_spectrum_to_molnupiravir.py on the output of step 1 to identify the sequences with molnupiravir-associated mutations. Provide the GISAID metadata file to include date data.

## Step 3:
Run step_3_molnupiravir_to_tree.py on the output of step 2 to count results by country and age.

# Example outputs
We include the `epi_isl_output.txt` and `aggregated.tsv` files for those who don't want to run the scripts. The `epi_isl_output.txt` file contains EPI_ISL ids for potentially molnupiravir-associated sequences identified through this method.