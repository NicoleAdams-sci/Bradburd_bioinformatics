-- -*- lua -*-
help([[
Name: feems
Version: 1.0.0
Description: Fast Estimation of Effective Migration Surfaces (feems) is a Python package implementing a statistical method for inferring and visualizing gene-flow in spatial population genetic data.
]])

-- Load necessary modules
load("gcc")
load("geos")
load("mamba")

-- Set the FEEMS environment variable
local feems_dir = "/nfs/turbo/lsa-bradburd/shared/programs/feems"
setenv("FEEMS", feems_dir)

-- Set the Conda environment paths
local conda_bin = pathJoin(feems_dir, "bin")
local conda_lib = pathJoin(feems_dir, "lib")
local conda_python = pathJoin(feems_dir, "lib/python3.8/site-packages")  -- Update Python version if necessary

-- Add the bin, lib, and python paths of the Conda environment
prepend_path("PATH", conda_bin)
prepend_path("LD_LIBRARY_PATH", conda_lib)
prepend_path("PYTHONPATH", conda_python)

-- Ensure necessary Conda environment variables are set
setenv("CONDA_PREFIX", feems_dir)
setenv("CONDA_DEFAULT_ENV", feems_dir)
setenv("CONDA_PYTHON_EXE", pathJoin(conda_bin, "python"))
setenv("CONDA_EXE", pathJoin(conda_bin, "mamba"))

-- Set whatis and description
whatis("Name: feems")
whatis("Version: 1.0.0")
whatis("Description: Fast Estimation of Effective Migration Surfaces (feems)")
