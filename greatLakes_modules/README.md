# Custom modules available on Great Lakes

**Use the following to initiate all custom modules**
`module use /nfs/turbo/lsa-bradburd/shared/Lmod/` 

## ANGSD
version: v0.940
to use: `module load angsd/0.940`

## Structure Threader (fastStructure)
version:1.0
to use: `module load structure_threader/1.0`
**Need to specify fastStructure location**
fastStructure: "/nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/fastStructure"
structure: "/nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/structure"
structure threader: "/nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/structure_threader"
example run: `structure_threader run -K <K> -R <R> -i <input.fam> --ind <indiv_list.txt> -o <out.st> -t <t> -fs /nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/fastStructure`

## SLiM
version: 4.3
to use: `module load SLiM/4.3`

## VCF2PCACluster
version: 1.40
to use: `module load VCF2PCACluster/1.40`