Regression analysis
========================

Since version 0.15, we have also the ability to perform a regression analysis at
the drug level. 

Information on the analysis and code can be found in the reference:
:mod:`gdsctools.regression.`, in a notebook (see https://github.com/CancerRxGene/gdsctools/tree/master/notebooks).


We first show how to perform the analysis. We then show how to use a snakemake
pipeline to perform the analysis on a cluster to get also nice HTML reports.


Documentation coming soon



Snakemake pipeline
--------------------


in pipelines/, there is a regression.rules and a config.yaml file to be copied
in a local directory.

Then, type::

    snakemake -s regression.rules -j 4

Or a cluster, you may add the following information (for instance on a slurm
system)::


    snakemake -s regression.rules -j 40 --cluster "sbatch --qos normal"

where -j 40 means uses 40 cores. When until it is finished. You should have an
index.html


See an example in _static/regression


