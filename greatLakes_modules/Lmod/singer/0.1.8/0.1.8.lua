-- -*- lua -*-
help([[
Name: SINGER
Version: 0.1.8
Description: Sampling and INference of GEnealogies with Recombination. A Bayesian method to do posterior sampling of Ancestral Recombination Graph under Sequentially Markovian Coalescent
]])

-- Load required modules
load("slim-postprocess/1.0.0")

-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/programs/"

-- Define the bin directory path
local bindir = pathJoin(base, "singer_0.1.8")

-- Set the PATH environment variable
prepend_path("PATH", bindir)

