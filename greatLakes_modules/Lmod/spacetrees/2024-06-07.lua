-- -*- lua -*-
help([[
SpaceTrees

Ppatial phylogenetic inference using tree sequences,
tsconvert, tskit, Snakemake, and Relate.

Usage:
  module use /nfs/turbo/lsa-bradburd/shared/Lmod
  module load spacetrees/2024-06-07

This module loads Relate and the shared SpaceTrees Python virtual
environment.

Important:
  This module is intended for running SpaceTrees, not for installing/updating
  Python packages. For pip installs, use the build environment documented
  separately, because relate/1.2.4 loads gcc/15.1.0.
]])

whatis("Name: SpaceTrees")
whatis("Version: 2024-06-07")
whatis("Description: Spatial phylogenetic inference workflows using tree sequences and Relate")
whatis("URL: https://github.com/osmond-lab/spacetrees")

-- Avoid loading multiple SpaceTrees versions at once
conflict("spacetrees")

-- Load runtime dependency
depends_on("relate/1.2.4")

local base = "/nfs/turbo/lsa-bradburd/shared/programs/spacetrees/2024-06-07"
local repo = pathJoin(base, "spacetrees")
local env  = pathJoin(base, "venv")

-- put venv/bin first so python, pip, snakemake, etc. come from the SpaceTrees env.
prepend_path("PATH", pathJoin(env, "bin"))

-- Useful environment variables
setenv("SPACETREES_HOME", base)
setenv("SPACETREES_REPO", repo)
setenv("SPACETREES_VENV", env)
setenv("VIRTUAL_ENV", env)

prepend_path("PYTHONPATH", repo)

-- For plotting
local jupyter = pathJoin(base, "jupyter", "share", "jupyter")
prepend_path("JUPYTER_PATH", jupyter)

-- Avoid issues from inherited Python settings
unsetenv("PYTHONHOME")