-- -*- lua -*-
help([[
Name: feems
Version: 1.0.0
Description: Fast Estimation of Effective Migration Surfaces (feems) is a python package implementing a statistical method for inferring and visualizing gene-flow in spatial population genetic data.

Initial One-Time Setup: 
module load mamba
mamba init
source ~/.bashrc
]])

-- Load necessary modules
load("gcc")
load("geos")
load("mamba")


-- Set the FEEMS environment variable
local feems_dir = "/nfs/turbo/lsa-bradburd/shared/programs/feems"
setenv("FEEMS", feems_dir)

-- Activate the feems Conda environment
-- Using the `execute{cmd = ...}` method to source the shell setup files and activate the environment within this module.
execute{cmd = 'source ~/.bashrc && mamba activate ' .. feems_dir, modeA = {"load"}}

-- Optional: Add bin directory of the Conda environment to the PATH, if it's not automatically added
prepend_path("PATH", pathJoin(feems_dir, "bin"))

