# PeptoVar

PeptoVar - **Pept**ides **o**f **Var**iations

## Overview

PeptoVar is a program for annotation of genomic variants in protein coding genes and generation of variant peptides.

- Annotation of variants:
   - comprehensive information of potential compound effects of combining alternate alleles across multiple variant loci

- Variant peptides:
   - variant peptides of a population
   - personalized peptidome for a sample (phased and unphased data can be used)
   - unique peptides for a pair of samples


## The program was developed under the research project supported by RFBR according to the project number 17-04-02186.

## Requirements

* Linux or MacOS
* python >= 3.5
* pysam module >= 0.11.2.2

## Installation

* install pysam package using pip:

   pip3 install pysam

* [download](https://github.com/open-projects/PeptoVar/zipball/master) latest stable PeptoVar build from this page
* unzip the archive
* add resulting folder to your ``PATH`` variable
  * or add symbolic link for ``PeptoVar`` script to your ``bin`` folder
  * or use PeptoVar directly by specifying full path to the executable script

## Examples of usage

#### Annotation of variations
The example of annotation of genomic variations in the protein coding regions:

    PeptoVar.py -var nonsyn -gff ./testdata/test.gff -vcf ./testdata/test.vcf.gz
    
or
    
    PeptoVar.py -var nonsyn -indir ./testdata


#### Peptides for combinations of variations in the set
This example illustrates usage virtual sample with all variations:

    PeptoVar.py -peptlen 9 10 11 -gff ./testdata/test.gff -vcf ./testdata/test.vcf.gz
    
or
    
    PeptoVar.py -peptlen 9 10 11 -indir ./testdata
    


#### Peptides for a sample
This example illustrates usage for sample SAMPLE01:

    PeptoVar.py -samples SAMPLE01 -peptlen 9 -gff ./testdata/test.gff -vcf ./testdata/test.vcf.gz
    
or
    
    PeptoVar.py -samples SAMPLE01 -peptlen 9 -indir ./testdata


#### Unique peptides for a pair of samples
This example illustrates usage for the pair of samples: SAMPLE01 and SAMPLE02:

    PeptoVar.py -samples SAMPLE01 SAMPLE02 -peptlen 9 -gff ./testdata/test.gff -vcf ./testdata/test.vcf.gz
    
or
    
    PeptoVar.py -samples SAMPLE01 SAMPLE02 -peptlen 9 -indir ./testdata


## Documentation

Detailed PeptoVar description can be found in the [manual](https://github.com/open-projects/PeptoVar/blob/master/UserManual.pdf)

If you haven't found the answer to your question in the docs, or have any suggestions concerning new features, feel free to create an issue here, on GitHub, or write an email to dmitry.malko at gmail.com:
<br />![my mail](https://user-images.githubusercontent.com/5543031/28415000-8bea641e-6d56-11e7-85ca-4287500a4192.png)

## License
Copyright (c) 2017, D. Malko
All Rights Reserved

PeptoVar is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).


