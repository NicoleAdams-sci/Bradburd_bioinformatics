-- -*- lua -*-
help([[
tsbrowse v0.0.3: A command-line interface for visualizing tree sequences.
This module uses a standalone virtual environment.

Included Packages:
- tskit: 1.0.0  (Strict validation enabled)
- zarr:  2.18.4 (Pinned for compatibility)]])

local base = "/nfs/turbo/lsa-bradburd/shared/programs/tsbrowse_v0.0.3"

-- We NO LONGER load slim-postprocess here. 
-- This environment is completely self-contained.

prepend_path("PATH", pathJoin(base, "env/bin"))

-- Still a good idea to prevent home-directory interference
setenv("PYTHONNOUSERSITE", "1")

if (mode() == "load") then
    LmodMessage("tsbrowse v0.0.3 with tskit v1.0.0")
    
end
