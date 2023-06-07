### Calculating spectrum for molnupiravir-treated patients

Initially extracted mutations and combined into SBS spectra using:
```
python3 extract_MOV_second_trial.py -v refs_results/*.tsv
```

Then rescaled the spectrum using MutTui (https://github.com/chrisruis/MutTui):
```
python3 rescale_spectrum.py -s MOV_spectrum.csv -r wuhan-hu-1.fasta -o MOV_spectrum_rescaled.csv --rna
```

Spectrum can be plotted with MutTui:
```
MutTui plot -s MOV_spectrum_rescaled.csv -o MOV_spectrum_rescaled.pdf --rna --plot_proportion
```
