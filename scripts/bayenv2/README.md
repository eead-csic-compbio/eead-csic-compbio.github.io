
## bayenv2 scripts 

These scripts were first published as part of the 
[barley-agroclimatic-association](https://eead-csic-compbio.github.io/barley-agroclimatic-association/) repository.  

They support the parallel computation of bayenv2 Bayes factors and XtX provided that a 
[bayenv2 binary](https://bitbucket.org/tguenther/bayenv2_public/src) 
is placed in the same folder. Note that a copy of perl module
[ForkManager](https://metacpan.org/pod/Parallel::ForkManager) is also included.

Example of use:

    perl calc_bf_parallel.pl 

    Usage: perl calc_bf_parallel.pl <number of processors/cores> <bayenv2 command>

    <number of processors/cores> : while 1 is accepted, at least 2 should be requested
    <bayenv2 command> : is the explicit command that you would run in your terminal

    Example: bayenv2 -t -i input.tsv -p 20 -e env.tsv -n 87 -m matrix.txt -k 100000 -r 12345 -c -o test.env

Note: make sure you split env.tsv if there are more than 100 variables.
