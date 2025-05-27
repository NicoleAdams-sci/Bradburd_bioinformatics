# Custom modules available on Great Lakes

**Need to use the following to initiate all custom modules**<br>
`module use /nfs/turbo/lsa-bradburd/shared/Lmod/` 

## ANGSD
version: 0.940<br>
to use: `module load angsd/0.940`<br>
<br>
[ANGSD website](https://www.popgen.dk/angsd/index.php/ANGSD)<br>

## fastp
version: 0.24.0<br>
to use: `module load fastp`<br>
<br>
[fastp website](https://github.com/OpenGene/fastp)<br>

## feems
version: 1.0.0<br>
to use: `module load feems`<br>
<br>
[feems website](https://github.com/NovembreLab/feems)<br>

## locator
version: 1.2<br>
to use: `module load locator`<br>
<br>
[locator website](https://github.com/kr-colab/locator)<br>

## pixy
version: 1.2.10.beta2<br>
to use: `module purge; module load pixy`<br>
<br>
[pixy website](https://pixy.readthedocs.io/en/latest/index.html)<br>

## rad_haplotyper
version: 1.1.9<br>
to use: `module purge; module load rad_haplotyper`<br>
<br>
[rad_haplotyper website](https://github.com/chollenbeck/rad_haplotyper)<br>

## SLiM
Available versions: 4.3 (default), 5.0<br>
### SLiM 4.3
to use: `module load SLiM` or `module load SLiM/4.3/4.3`<br>
### SLiM 5.0
to use: `module load SLiM/5.0/5.0`<br>
<br>
[SLiM website](https://messerlab.org/slim/)<br>

## SLiM post processing
purpose: SLiM simulation output analysis and tree sequence processing<br>
version: 1.0.0<br>
to use: `module load slim-postprocess/1.0.0`<br>
contains: tskit v0.6.4, pyslim v1.0.4, msprime v1.3.4 (Python) + gaia (R)<br>
<br>
Links: [tskit](https://tskit.dev/) | [pyslim](https://pyslim.readthedocs.io/) | [msprime](https://msprime.readthedocs.io/) | [gaia](https://github.com/blueraleigh/gaia)<br>

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
version: 1.41<br>
to use: `module load Bioinformatics VCF2PCACluster`<br>
<br>
[VCF2PCACluster website](https://github.com/hewm2008/VCF2PCACluster)<br>