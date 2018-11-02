# Juicebox utilities

This repo includes a set of scripts that we use in the course of interacting with the [Juicebox tool](http://aidenlab.org/juicebox/). Juicebox is helpful for visualizing and manually editing assemblies based on Hi-C data.

We provide these scripts in the hope that they are helpful to others in interacting with Hi-C data. 

## Included Tools
The following scripts are provided in the `juicebox_scripts/` directory:

### agp2assembly.py
From an AGP file, create a .assembly file for use as a Map Assembly or Modified Assembly in Juicebox.

Usage:
```
agp2assembly.py usage:	agp2assembly.py <input_agp_file> <output_assembly_file>
```

### juicebox_assembly_converter.py
Given an original .fasta file and a .assembly file which you have created by modifying and exporting from Juicebox,
generate new .fasta, .agp, and .bed files describing the modified assembly. Also generates a report describing any
contig breaks which were introduced via Juicebox.

Usage:
```
usage: juicebox_assembly_converter.py [-h] -a ASSEMBLY -f FASTA [-p PREFIX] [-c] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -a ASSEMBLY, --assembly ASSEMBLY
                        juicebox assembly file
  -f FASTA, --fasta FASTA
                        the fasta file
  -p PREFIX, --prefix PREFIX
                        the prefix to use for writing outputs. Default: the
                        assembly file, minus the file extension
  -c, --contig_mode     ignore scaffold specification and just output contigs.
                        useful when only trying to obtain a fasta reflecting
                        juicebox breaks. Default: False
  -v, --verbose         print summary of processing steps to stdout, otherwise
                        silent. Default: True
```

### makeAgpFromFasta.py
Given a FASTA file, make an AGP file that represents it.

Usage:
```
makeAgpFromFasta.py usage:	makeAgpFromFasta.py <fasta_file> <agp_out_file>
```

## Example Workflows
Here are some of the common tasks that one might want to use these scripts for.

### Generating new SCAFFOLD `.fasta`, `.agp`, and `.bed` files plus making a break report after modifying an assembly in Juicebox
This will produce assembly outputs that match the scaffolding represented in Juicebox.
```
python juicebox_assembly_converter.py -a <modified_assembly> -f <original_premodification_fasta>
```


### Generating new CONTIG `.fasta`, `.agp`, and `.bed` files plus making a break report after modifying an assembly in Juicebox
This will produce assembly outputs that ONLY reflect the contigs present in the assembly. They WILL NOT match the scaffolding represented in Juicebox. This is most useful when you have used Juicebox to break possible misjoined contigs,
and you just want to get the broken contigs back with no scaffolding.
```
python juicebox_assembly_converter.py -a <modified_assembly> -f <original_premodification_fasta> -c
```

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

## Dependencies
The examples above require:

* Linux OS
* Python 2.7, 3.5, or 3.6
* [3d-dna](https://github.com/theaidenlab/3d-dna)
* [matlock](https://github.com/phasegenomics/matlock)
