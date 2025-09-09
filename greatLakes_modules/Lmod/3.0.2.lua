-- -*- lua -*-
help([[
Name: SCAT
Version: 3.0.2
Description: Smoothed and Continuous Assignment Tests
]])
-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/programs/scat"
-- Define the bin directory path
local bindir = pathJoin(base, "src")
-- Set the PATH environment variable
prepend_path("PATH", bindir)
