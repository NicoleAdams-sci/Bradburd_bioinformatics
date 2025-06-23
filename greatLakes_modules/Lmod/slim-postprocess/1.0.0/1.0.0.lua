-- -*- lua -*-
help([[
Post-processing environment for SLiM simulation outputs.
Python: tskit, pyslim, msprime
R: gaia (requires R module)
]])

whatis("Name: slim-postprocess")
whatis("Version: 1.0.0")

local base = "/nfs/turbo/lsa-bradburd/shared/programs/slim-postprocess_v1.0.0"

-- Load R module dependency
load("R/4.4.0")  -- or appropriate version

-- Set up Python environment
prepend_path("PATH", pathJoin(base, "env/bin"))
setenv("SLIM_POSTPROCESS_ROOT", base)

-- Optional: Set R library path if you install gaia to module directory
-- prepend_path("R_LIBS_USER", pathJoin(base, "R_libs"))

-- Display version info when loading
if (mode() == "load") then
    LmodMessage("Loading slim-postprocess v1.0.0")
    LmodMessage("  Python packages: tskit v0.6.4, pyslim v1.0.4, msprime v1.3.4, fastgaia-0.1.1")
    LmodMessage("  R packages: gaia")
    LmodMessage("  R version: 4.4.0")
end
