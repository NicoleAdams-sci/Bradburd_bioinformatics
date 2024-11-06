# Custom modules available on Great Lakes

**Use the following to initiate all custom modules**
`module use /nfs/turbo/lsa-bradburd/shared/Lmod/` 

## ANGSD
version: v0.940<br>
to use: `module load angsd/0.940`

## SLiM
version: 4.3<br>
to use: `module load SLiM/4.3`

## Structure Threader (fastStructure)
version:1.0<br>
to use: `module load structure_threader/1.0`<br>
<br>
**Need to specify fastStructure location**<br>
fastStructure: "/nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/fastStructure"<br>
structure: "/nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/structure"<br>
structure threader: "/nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/structure_threader"<br>
<br>
example run: `structure_threader run -K <K> -R <R> -i <input.fam> --ind <indiv_list.txt> -o <out.st> -t <t> -fs /nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/fastStructure`


## VCF2PCACluster
version: 1.40<br>
to use: `module load VCF2PCACluster/1.40`