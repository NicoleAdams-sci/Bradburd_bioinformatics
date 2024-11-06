-- -*- lua -*-
help([[
Name: angsd
Version: 0.940
Description: Analysis of next generation Sequencing Data (ANGSD)
]])

whatis("Name: ANGSD")
whatis("Version: 0.940")

-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/programs/angsd0.940"

-- Set the PATH environment variable
prepend_path("PATH", base)
