-- -*- lua -*-
help([[
Name: fastp
Version:  0.24.0
Description: Ultrafast all-in-one preprocessing and quality control for short-read FastQ data
]])

whatis("Name: fastp")
whatis("Version: 0.24.0")

-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/programs/fastp"

-- Set the PATH environment variable
prepend_path("PATH", base)
