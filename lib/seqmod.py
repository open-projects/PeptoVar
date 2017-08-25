#!/usr/bin/env python3

# Copyright (C) 2017 Dmitry Malko
# This file is part of PeptoVar (Peptides of Variations): the program for personalized and population-wide peptidome generation.
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
from lib.seqtools import translate, complement

# The package is used for debugging only

class Seq:
    def __init__(self, strand, beg, end, seq):
        self._seq = list(seq)
        self._ins = list() # insertion at the begin of exon
        self.strand = strand
        self.beg = beg
        self.end = end
    
    def setSeq(self, lbeg, lend, seq):
        llen = lend - lbeg + 1
        if llen > 0:
            seq = list(seq)
            if len(seq) < llen:
                seq.extend('-' for i in range (llen - len(seq)))
            elif len(seq) > llen:
                tail = seq[llen - 1:]
                seq = seq[:llen - 1]
                seq.append(tail)
            if len(seq) != llen:
                print("ERROR: wrong allele_seq length {}..{}".format(lbeg, lend))
                exit()
        elif llen == 0:
            if lend < 0:
                self._ins = list(seq)
            else:
                self._seq[lend] = list(self._seq[lend] + seq)
            return None
        else:
            print("ERROR: seq.beg > seq.end {}..{}".format(lbeg, lend))
            exit()
        
        for i in range(len(seq)):
            if type(self._seq[lbeg + i]) is list: # there are insertions in the sequence
                if self._seq[lbeg + i][0] != '-': # some VCF files have wrong overlapped alleles
                    self._seq[lbeg + i][0] = seq[i]
            else:
                if self._seq[lbeg + i] != '-': # some VCF files have wrong overlapped alleles
                    self._seq[lbeg + i] = seq[i]
        return None
    
    def getSeq(self, beg = 0, end = None):
        seq = ''
        if end is None:
            end = len(self._seq)
            seq = ''.join(self._ins)
        else:
            end += 1
        for nucl in self._seq[beg:end]:
            if type(nucl) is list:
                seq += ''.join(nucl)
            elif nucl != '-':
                seq += nucl
        return seq
# end of Seq

class Sequence:
    def __init__(self, chrom, trn_id, samples):
        self.chrom = chrom
        self.trn_id = trn_id
        self.strand = None
        self._smplseq = {}
        self._refseq = []
        self._smplprot = {}
        self._smplpept = {}
        self.refprot = ''
        self.samples = samples
        self._strangers = {}
    
    def append(self, strand, beg, end, seq):
        if not self.strand:
            self.strand = strand
        elif self.strand != strand:
            print("ERROR: Wrong strand!")
        for sample in self.samples:
            if sample.id not in self._smplseq:
                self._smplseq[sample.id] = []
            self._smplseq[sample.id].append(Seq(strand, beg, end, seq))
        self._refseq.append(Seq(strand, beg, end, seq))
    
    def modify(self, snp):
        for sample in self.samples:
            for allele in snp.getSampleAlleles(sample):
                allele_seq = complement(allele.seq) if snp.strand == '-' else allele.seq
                if allele.isReference():
                    if allele_seq != '-':
                        for exseq in self._refseq:
                            if exseq.beg <= snp.beg and snp.end <= exseq.end:
                                lbeg = snp.beg - exseq.beg
                                lend = snp.end - exseq.beg
                                refseq = exseq.getSeq(lbeg, lend)
                                if refseq.upper() != allele_seq:
                                    print("Oy Gevalt! Wrong reference allale!")
                                break
                else:
                    for exseq in self._smplseq[sample.id]:
                        if exseq.beg <= snp.beg and snp.end <= exseq.end:
                            lbeg = snp.beg - exseq.beg
                            lend = snp.end - exseq.beg
                            exseq.setSeq(lbeg, lend, allele_seq)
                            break
    
    def translate(self):
        self._refseq = sorted(self._refseq, key = lambda x: x.beg)
        seq = "".join(x.getSeq() for x in self._refseq)
        if self.strand == '-':
            self.refprot = (translate(complement(seq[::-1])))[0]
        else:
            self.refprot = (translate(seq))[0]
        
        for sample in self.samples:
            self._smplseq[sample.id] = sorted(self._smplseq[sample.id], key = lambda x: x.beg)
            
            seq = "".join(x.getSeq() for x in self._smplseq[sample.id])
            if self.strand == '-':
                protein = (translate(complement(seq[::-1])))[0]
                protein = re.sub(r"\*.*", "", protein)
                self._smplprot[sample.id] = protein
            else:
                protein = (translate(seq))[0]
                protein = re.sub(r"\*.*", "", protein)
                self._smplprot[sample.id] = protein
        return self._smplprot
    
    def peptides(self, length_set):
        for sample_id in self._smplprot:
            for length in length_set:
                if length: # length == 0 - full length protein
                    for i in range(len(self._smplprot[sample_id])):
                        pept = self._smplprot[sample_id][i : i + length]
                        if len(pept) == length:
                            if sample_id not in self._smplpept:
                                self._smplpept[sample_id] = set()
                            self._smplpept[sample_id].add(pept)
                        else:
                            break
    
    def refpeptides(self, length):
        peptset = set()
        for i in range(len(self.refprot)):
            pept = self.refprot[i : i + length - 1]
            if len(pept) == length:
                peptset.add(pept)
            else:
                break
        return peptset
    
    def add2compare(self, sample_peptides):
        for item in sample_peptides:
            sample_id = "_".join((item[2], str(item[3]), str(item[4])))
            peptide = item[9]
            if sample_id not in self._strangers:
                self._strangers[sample_id] = []
            self._strangers[sample_id].append(peptide)
    
    def compare(self):
        misspept = []
        for sample_id, peptides in self._strangers.items():
            for peptide in peptides:
                if sample_id not in self._smplpept or peptide not in self._smplpept[sample_id]:
                    misspept.append("graphpept\t{}\t{}\t{}\t{}\n".format(self.chrom, self.trn_id, sample_id, peptide))
        for sample_id, peptides in self._smplpept.items():
            for peptide in peptides:
                if sample_id not in self._strangers or peptide not in self._strangers[sample_id]:
                    misspept.append("dummypept\t{}\t{}\t{}\t{}\n".format(self.chrom, self.trn_id, sample_id, peptide))
        return misspept
# end of Sequence
