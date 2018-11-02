# Juicebox utilities

This repo includes a set of scripts that we use in the course of interacting with the [Juicebox tool](http://aidenlab.org/juicebox/). Juicebox is helpful for visualizing and manually editing assemblies based on Hi-C data.

We provide these scripts in the hope that they are helpful to others in interacting with Hi-C data. 

## Dependencies
These examples require:

* Linux OS
* Python 2.7, 3.5, or 3.6
* [3d-dna](https://github.com/theaidenlab/3d-dna)
* [matlock](https://github.com/phasegenomics/matlock)

## Workflows
Here are some of the common tasks that one might want to use these scripts for.

### Getting a `.assembly` file from an AGP file
```
python agp2assembly.py in.agp out.assembly
```
### Getting an AGP file from a FASTA file
```
python makeAgpFromFasta.py in.fasta out.agp
```
### Getting a `.hic` file from a BAM file and a `.assembly` file. 
```
matlock bam2 juicer in.bam out.links.txt  # this step sometimes crashes on mem
sort -k2,2 -k6,6 out.links.txt > out.sorted.links.txt
bash 3d-dna/visualize/run-assembly-visualizer.sh -p false in.assembly out.sorted.links.txt # creates a .hic file
```

