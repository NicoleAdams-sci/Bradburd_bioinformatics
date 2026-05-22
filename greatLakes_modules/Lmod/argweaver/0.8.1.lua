-- -*- lua -*-
help([[
Name: ARGweaver
Version: 0.8.1
Description: Sampling and manipulating genome-wide ancestral recombination graphs (ARGs).
             Implements the SMC/ARG-sampling method of Rasmussen & Siepel (2013).

Key programs:
  arg-sample        -- infer ARGs from sequence data (main inference)
  arg-sim           -- simulate ARGs and sequences
  arg-summarize     -- compute summary statistics from sampled ARGs
  arg-extract-tmrca -- extract TMRCA across sampled ARGs

Citation: Rasmussen & Siepel (2013) arXiv:1306.5110
]])

whatis("Name: ARGweaver")
whatis("Version: 0.8.1")
whatis("Description: Genome-wide ancestral recombination graph sampling")


local base = "/nfs/turbo/lsa-bradburd/shared/programs/argweaver"

prepend_path("PATH", pathJoin(base, "bin"))
prepend_path("PYTHONPATH", pathJoin(base, "lib64/python2.7/site-packages"))
