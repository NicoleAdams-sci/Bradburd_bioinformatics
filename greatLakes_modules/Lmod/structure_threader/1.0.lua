-- -*- lua -*-
help([[
Name: structure_threader
Version: 1.0
Description: structure_threader: A program to parallelize and automate the runs of Structure, fastStructure, MavericK and ALStructure software.

To run fastStructure: 
structure_threader run -K <K> -R <R> -i <input.fam> --ind <indiv_list.txt> -o <out.st> -t <t> -fs /nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/fastStructure

where K = No. of groups; R = No. of replicates; t = No. of threads 

input files: Plink .fam file and a text file of sample IDs

]])


whatis("Name: structure_threader")
whatis("Version: 1.0")
whatis("To run fastStructure: structure_threader run -K <K> -R <R> -i <input.fam> --ind <indiv_list.txt> -o <out.st> -t <t> -fs /nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/fastStructure")


-- Define the base installation directory
local base = "/nfs/turbo/lsa-bradburd/shared/programs/structure_threader"

-- Define the bin directory path
local bindir = pathJoin(base, "bin")

-- Define the directory containing the structure_threader package
local pypath = pathJoin(base, "lib/python3.10/site-packages")

-- Set the PATH environment variable
prepend_path("PATH", bindir)

-- Set the PYTHONPATH environment variable
prepend_path("PYTHONPATH", pypath)
