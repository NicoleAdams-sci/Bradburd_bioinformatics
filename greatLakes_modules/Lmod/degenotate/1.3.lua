-- -*- lua -*-
help([[
Degenotate

Annotate degeneracy of sites in coding regions of a genome.

Usage:
  degenotate
  degenotate.py
]])

whatis("Name: Degenotate")
whatis("Version: 1.3")
whatis("Description: Annotate degeneracy of sites in coding regions of a genome")
whatis("URL: https://github.com/harvardinformatics/degenotate")

local base = "/nfs/turbo/lsa-bradburd/shared/programs/degenotate"
local env  = pathJoin(base, "env")

prepend_path("PATH", pathJoin(env, "bin"))

setenv("DEGENOTATE_HOME", base)
