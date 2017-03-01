There are two snakemake pipelines.


One called anova_pipelines.rules, which does not need any config file and is
provided as an example.

The second one is called regression_pipelines.rules and rely on a 
configuration file (config.yaml)


For instance all parameters related to the regression may be placed in ::

regression:
   - param1
   - param2
   - param3
