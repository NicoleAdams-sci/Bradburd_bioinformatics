-- -*- lua -*-

help([[
MalariaGEN Data 15.8.0

Python package for accessing and analyzing MalariaGEN vector genomic
surveillance data, including the Af1 resource.

Usage:
  module load malariagen_data/15.8.0

Protected datasets require per-user Google Application Default
Credentials and appropriate MalariaGEN authorization.

Documentation:
  https://malariagen.github.io/vector-data/
  https://malariagen.github.io/vector-data/af1/cloud.html
]])

whatis("Name: malariagen_data")
whatis("Version: 15.8.0")
whatis("Description: Access and analyze MalariaGEN vector genomic data")
whatis("URL: https://malariagen.github.io/vector-data/")

local release = "/nfs/turbo/lsa-bradburd/shared/programs/malariagen_data/2026-02"
local venv    = pathJoin(release, "venv")

conflict("malariagen_data")
depends_on("python/3.12.1")

prepend_path("PATH", pathJoin(venv, "bin"))

setenv("VIRTUAL_ENV", venv)
setenv("MALARIAGEN_DATA_HOME", release)
setenv("PYTHONNOUSERSITE", "1")

