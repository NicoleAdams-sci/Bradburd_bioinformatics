-- -*- lua -*-
help([[
Name: locator
Version: 1.2
Description: 
]])

whatis("Name: Locator")
whatis("Version: 1.2")

-- Load necessary modules


-- Define the base installation directory
local envPath = "/nfs/turbo/lsa-bradburd/shared/programs/locator_env"

-- Set environment paths
prepend_path("PATH", pathJoin(envPath, "bin"))
prepend_path("MANPATH", pathJoin(envPath, "share/man"))
prepend_path("LD_LIBRARY_PATH", pathJoin(envPath, "lib"))
prepend_path("LIBRARY_PATH", pathJoin(envPath, "lib"))
prepend_path("PYTHONPATH", pathJoin(envPath, "lib/python3.11/site-packages"))

-- Set environment variables
setenv("CONDA_PREFIX", envPath)
setenv("CONDA_DEFAULT_ENV", "locator")
