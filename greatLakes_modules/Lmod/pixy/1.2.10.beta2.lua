-- -*- lua -*-
help([[
Name: pixy
Version: 1.2.10.beta2
Description: A command line tool for calculating the population genetic summary statistics pi (average per site heterozygosity) and dxy (average number of nucleotide differences between populations per site) from a VCF file.

** Before load this module, run 'purge modules' **

input files: VCF that includes invariant sites, population or samples file]])


-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/conda-envs/pixy_e"

-- Prepend the bin directory to the PATH environment variable
prepend_path("PATH", pathJoin(base, "bin"))

-- Set the PYTHONPATH environment variable
prepend_path("PYTHONPATH", pathJoin(base, "lib/python3.8/site-packages"))

-- Set the conda prefix
setenv("CONDA_PREFIX", base)
