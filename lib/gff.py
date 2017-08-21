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
import os.path
from operator import itemgetter
from lib.fasta import Fasta

class gffItem:
    _pattern_id = re.compile("ID=([^;]+)")
    _pattern_parent_id = re.compile("Parent=([^;]+)")
    
    def __init__(self, fields):
        self.id = self._getID(fields[8])
        self.parent_id = self._getParentID(fields[8])
        self.seq = None
        self.chrom = fields[0]
        self.source = fields[1]
        self.type = fields[2]
        self.beg = int(fields[3])
        self.end = int(fields[4])
        self.score = fields[5]
        self.strand = fields[6]
        self.phase = int(fields[7]) if fields[7] != '.' else fields[7]
        self.attributes = fields[8]
        
    def _getID(self, string):
        id_search = re.search(self._pattern_id, string)
        if id_search:
            return id_search.group(1)
        return None
    
    def _getParentID(self, string):
        parent_search = re.search(self._pattern_parent_id, string)
        if parent_search:
            parents_id = parent_search.group(1)
            return parents_id.split(',')
        return None
# end of gffItem

class Gff:
    def __init__(self, gff_file, tmp_dir = '.'):
        self._transcripts = {}
        self._tmp = tmp_dir
        self._fasta = Fasta(tmp_dir)
        self._rawseq = {}
        
        pos = 0
        gff = None
        if os.path.isfile(gff_file):
            gff = open(gff_file, 'r')
        elif os.path.isfile(gff_file + '3'):
            gff = open(gff_file + '3', 'r')
        else:
            raise ValueError("no GFF file")
        
        for line in gff:
            pos += len(line)
            if re.match("##gff-version", line):
                for line in gff:
                    pos += len(line)
                    if re.match("##FASTA", line):
                        gff.close()
                        self._sort_exons()
                        self._fasta.appendFasta(gff_file, pos)
                        return
                    elif not re.match("#", line):
                        line = line.rstrip()
                        fields = line.split()
                        if len(fields) == 9 and fields[2] == 'CDS':
                            item = gffItem(fields)
                            if item.parent_id:
                                for pid in item.parent_id:
                                    if pid not in self._transcripts:
                                        self._transcripts[pid] = []
                                    self._transcripts[pid].append(item)
                            else:
                                raise ValueError("no 'Parent' attribute for CDS in locus {} ({}..{})".format(item.chrom, item.beg, item.end))
        gff.close()
        self._sort_exons()
    
    def _sort_exons(self):
        for exons in self._transcripts.values():
            exons.sort(key=lambda x: x.beg)
            if exons[-1].strand == '-':
                exons[-1].end -= exons[-1].phase
            else:
                exons[0].beg += exons[0].phase
    
    def attachSeq(self, file = None):
        if file:
            if re.search(r"\.seq$", file):
                self._fasta.appendSeq(file)
            elif re.search(r"\.fa(sta)?$", file):
                self._fasta.appendFasta(file)
        
        for exons in self._transcripts.values():
            for exon in exons:
                exon.seq = self._fasta.getSeq(exon.chrom, exon.beg, exon.end)
    
    def getTranscriptsID(self):
        return sorted(list(self._transcripts.keys()))
    
    def getTranscriptExons(self, trn_id):
        if trn_id in self._transcripts:
            return self._transcripts[trn_id]
        return []
# end of Gff