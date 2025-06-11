-- -*- lua -*-
help([[
Name: admixture
Version:  1.3.0
Description: maximum likelihood estimation of individual ancestries from multilocus SNP genotype datasets
]])

whatis("Name: admixture")
whatis("Version: 1.3.0")

-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/programs/admixture_1.3.0"

-- Set the PATH environment variable
prepend_path("PATH", base)
