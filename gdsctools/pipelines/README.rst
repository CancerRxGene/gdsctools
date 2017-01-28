the configutation file is called config.yaml (not regression.yaml) since
it may be used not only by regression.rules but any other pipelines as long
as we keep it structured. 

For instance all parameters related to the regression may be placed in ::

regression:
   - param1
   - param2
   - param3
