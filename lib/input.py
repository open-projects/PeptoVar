#!/usr/bin/env python3

# Copyright (C) 2017 D. Malko
# This file is part of PeptoVar (Peptides on Variations): the program for personalization of protein coding genes and population-wide peptidome generation.
#
# PeptoVar is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PeptoVar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PeptoVar.  If not, see <http://www.gnu.org/licenses/>.

import re
import glob
import os.path
from lib.output import ColorPrint

class InFile:
    def __init__(self, dir_path):
        if os.path.isdir(dir_path):
            self.bulk = self._globbing(dir_path)
        else:
            self.bulk = []
    
    def _globbing(self, path):
        input_set = []
        re.sub(path, '', r'\/$')
        for vcf_name in glob.glob(path + "/*.vcf.gz"):
            re_name = re.match(r'(.*\/(.+))\.vcf.gz$', vcf_name)
            if re_name:
                full_name = re_name.group(1)
                name = re_name.group(2)
                gff_name = None
                gff_name1 = full_name + ".gff"
                gff_name2 = full_name + ".gff3"
                fasta_name = None
                fasta_name1 = full_name + ".fasta"
                fasta_name2 = full_name + ".fa"
                if os.path.isfile(gff_name1):
                    gff_name = gff_name1
                elif os.path.isfile(gff_name2):
                    gff_name = gff_name2
                if gff_name:
                    if os.path.isfile(fasta_name1):
                        fasta_name = fasta_name1
                    elif os.path.isfile(fasta_name2):
                        fasta_name = fasta_name2
                    input_set.append({'name': name, 'vcf': vcf_name, 'gff': gff_name, 'fasta': fasta_name, 'seq': None})
                else:
                    cprint = ColorPrint()
                    cprint.printWarning("No {}.gff or {}.gff3 file ...skipped!".format(name))
        return input_set
# end of InFile