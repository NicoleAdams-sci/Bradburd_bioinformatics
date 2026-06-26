-- -*- lua -*-
help([[
ROHan

ROHan estimates local heterozygosity and runs of homozygosity from BAM files

Usage: rohan
]])

whatis("Name: ROHan")
whatis("Version: 20230216")
whatis("Description: Estimate local heterozygosity and runs of homozygosity from BAM files")
whatis("URL: https://github.com/grenaud/rohan")

local base = "/nfs/turbo/lsa-bradburd/shared/programs/rohan"
local env  = pathJoin(base, "rohan-env")
local src  = pathJoin(base, "rohan-src")
local bin  = pathJoin(src, "bin")

prepend_path("PATH", bin)
prepend_path("PATH", pathJoin(env, "bin"))

-- Needed if ROHan dynamically links to libraries from the Mamba environment
prepend_path("LD_LIBRARY_PATH", pathJoin(env, "lib"))

setenv("ROHAN_HOME", base)
setenv("ROHAN_SRC", src)
setenv("ROHAN_ENV", env)
setenv("ROHAN_BIN", bin)
setenv("ROHAN_EXE", pathJoin(bin, "rohan"))
