-- -*- lua -*-
help([[
Name: Relate
Version: 1.2.4
Description: Estimates genome-wide genealogies in the form of trees that adapt to changes in local ancestry caused by recombination. (ARG inference)
]])

-- Load required modules
load("Rtidyverse")

-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/programs/relate_v1.2.4/"

-- Define the bin directory path
local bindir = pathJoin(base, "bin")

-- Set the PATH environment variable
prepend_path("PATH", bindir)
