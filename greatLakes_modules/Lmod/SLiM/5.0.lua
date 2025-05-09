-- -*- lua -*-
help([[
Name: SLiM
Version: 5.0
Description: SLiM is an evolutionary simulation framework that combines a powerful engine for population genetic simulations with the capability of modeling arbitrarily complex evolutionary scenarios.]])

whatis("Name: SLiM")
whatis("Version: 5.0")

-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/programs/SLiM_v5/build/"

-- Set the PATH environment variable
prepend_path("PATH", base)
