-- -*- lua -*-
help([[
Name: SLiM
Version: 4.3
Description: SLiM is an evolutionary simulation framework that combines a powerful engine for population genetic simulations with the capability of modeling arbitrarily complex evolutionary scenarios.
]])

-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/programs/SLiM/"

-- Set the PATH environment variable
prepend_path("PATH", base)
