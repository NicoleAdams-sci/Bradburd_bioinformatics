# Custom modules available on Great Lakes

**Need to use the following to initiate all custom modules**<br>
`module use /nfs/turbo/lsa-bradburd/shared/Lmod/` 

## ANGSD
version: 0.940<br>
to use: `module load angsd/0.940`<br>
<br>
[ANGSD website](https://www.popgen.dk/angsd/index.php/ANGSD)<br>

## feems
version: 1.0.0<br>
to use: `module load feems`<br>
<br>
**Initial one-time setup**
```
module load mamba
mamba init
source ~/.bashrc
module load feems
```

[feems website](https://github.com/NovembreLab/feems)<br>

## pixy
version: 1.2.10.beta2<br>
to use: `module purge; module load pixy`<br>
<br>
[pixy website](https://pixy.readthedocs.io/en/latest/index.html)<br>

## SLiM
version: 4.3<br>
to use: `module load SLiM/4.3`<br>
<br>
[SLiM website](https://messerlab.org/slim/)<br>

## Structure Threader (fastStructure)
version: 1.0<br>
to use: `module load structure_threader/1.0`<br>
<br>
**Need to specify program location**<br>
fastStructure: "/nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/fastStructure"<br>
structure: "/nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/structure"<br>
structure threader: "/nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/structure_threader"<br>
<br>
example fastStructure run: `structure_threader run -K <K> -R <R> -i <input.fam> --ind <indiv_list.txt> -o <out.st> -t <t> -fs /nfs/turbo/lsa-bradburd/shared/programs/structure_threader/bin/fastStructure`<br>
<br>
[Structure threader website](https://structure-threader.readthedocs.io/en/latest/usage/)<br>

## VCF2PCACluster
version: 1.40<br>
to use: `module load VCF2PCACluster/1.40`<br>
<br>
[VCF2PCACluster website](https://github.com/hewm2008/VCF2PCACluster)<br>