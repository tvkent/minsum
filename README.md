# Coalescent simulations and haplotype scoring
Scripts from Sas et al.

### Scripts
* min_sum.py: python script which calculates various measures of the minimum sum of hapcut results. Paper uses minsum--other measures have autocorrelation problems
* minsum.sh: bash script to run python script above
* hapcut.sh: example script to get hapcut results
* simsum.py: python script to get minsum metric from simulation results
* hap_plots.R: R plotting scripts
* minsum_plots.R: R plotting scripts

ms commands were: ms 349 100000 -s 244 -r (5,20,100) 2000 -I 2 348 1 -n 2 0.05 -ej 0.04 2 1.
