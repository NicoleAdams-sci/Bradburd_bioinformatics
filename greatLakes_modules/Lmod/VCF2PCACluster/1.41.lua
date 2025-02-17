-- -*- lua -*-
help([[
Name: VCF2PCACluster
Version: 1.41
Description: VCF2PCACluster: Software to perform PCA and clustering analysis for population VCF file.
]])

-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/programs/VCF2PCACluster"
-- Define the bin directory path
local bindir = pathJoin(base, "bin")
-- Set the PATH environment variable
prepend_path("PATH", bindir)
