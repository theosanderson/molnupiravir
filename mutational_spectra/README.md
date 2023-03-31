# Description of files, scripts and commands in mutational spectra analyses
The location of files is shown in each section

### Calculation of long branch spectrum (Figure 7)
Files and scripts are in the long_phylogenetic_branches directory

Initially the SBS spectrum of the long phylogenetic branches was calculated:
python3 calculate_long_branch_spectrum.py -s sample_paths.txt -r NC_045512.2.fasta -o long_branch_spectrum.csv

Outputs a spectrum of mutations on branches with at least 10 mutations, this was then rescaled using MutTui (available at [https://github.com/chrisruis/MutTui](https://github.com/chrisruis/MutTui)):
python3 MutTui/MutTui/rescale_spectrum.py -s long_branch_spectrum.csv -r NC_045512.2 -o long_branch_spectrum_rescaled.csv --rna

The rescaled spectrum can be plotted with MutTui:
MutTui plot -s long_branch_spectrum_rescaled.csv -o long_branch_spectrum_rescaled.pdf --rna --plot_proportion
