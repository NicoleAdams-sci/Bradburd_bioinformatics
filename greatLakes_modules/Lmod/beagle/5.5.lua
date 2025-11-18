-- -*- lua -*-
help([[
Name: Beagle
Version: 5.5
Description: Beagle is a software package for phasing genotypes and imputing ungenotyped markers
Usage: beagle_5.5 [options]
Memory: Set BEAGLE_MEMORY environment variable to change memory allocation (default: 4000m)
]])

whatis("Name: Beagle")
whatis("Version: 5.5")
whatis("Description: Software for phasing genotypes and imputing ungenotyped markers")
whatis("URL: https://faculty.washington.edu/browning/beagle/beagle.html")

-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/programs/beagle"

-- Add the directory containing the wrapper script to PATH
prepend_path("PATH", base)

-- Optionally set a default memory allocation
-- setenv("BEAGLE_MEMORY", "4000m")
