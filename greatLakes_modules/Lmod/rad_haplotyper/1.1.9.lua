-- -*- lua -*-
help([[
Name: rad_haplotyper
Version:  1.1.9
Description: Produce SNP haplotypes from RAD-seq data with fixed-size RAD loci
]])

whatis("Name: rad_haplotyper")
whatis("Version: 1.1.9")

-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/programs/rad_haplotyper_e/bin/"

-- Set the PATH environment variable
prepend_path("PATH", base)
