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

codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

def is_stopcodon(codon):
    codon = codon.upper()
    if codon in codontable and codontable[codon] == '*':
        return True
    return False

def translate(seq):
    seq_len = len(seq)
    protein = ''
    tail = None
    for i in range(0, seq_len, 3):
        if seq_len < i + 3:
            tail = seq[i:]
            break
        codon = seq[i:i+3].upper()
        if codon in codontable:
            protein += codontable[codon]
        else:
            protein += 'X'
    return (protein, tail)
# end of translate_dna

trantab = str.maketrans('ACGTMRYKVHDBacgtmrykvhdb', 'TGCAKYRMBDHVtgcakyrmbdhv')
def complement(seq):
    return seq.translate(trantab)
# end of complement

class PeptComparator:
    def __init__(self):
        self.refpept = {}
        self.altpept = {}
    
    def addRefPept(self, pseq):
        length = len(pseq)
        if length not in self.refpept:
            self.refpept[length] = {}
        if pseq not in self.refpept[length]:
            self.refpept[length][pseq] = 0
        self.refpept[length][pseq] += 1
    
    def addAltPept(self, pseq):
        length = len(pseq)
        if length not in self.altpept:
            self.altpept[length] = {}
        if pseq not in self.altpept[length]:
            self.altpept[length][pseq] = 0
        self.altpept[length][pseq] += 1
    
    def getRefNotInAlt(self, length = None):
        pept = {}
        if length and (length not in self.refpept or length not in self.altpept):
            return pept
        for peptlen in self.refpept:
            if length and peptlen != length:
                continue
            for pseq in self.refpept[peptlen]:
                if pseq not in self.altpept[peptlen]:
                    pept[pseq] = self.refpept[peptlen][pseq]
                else:
                    if self.altpept[peptlen][pseq] < self.refpept[peptlen][pseq]:
                        pept[pseq] = self.refpept[peptlen][pseq] - self.altpept[peptlen][pseq]
            return pept
    
    def getAltNotInRef(self, length = None):
        pept = {}
        if length and (length not in self.refpept or length not in self.altpept):
            return pept
        for peptlen in self.refpept:
            if length and peptlen != length:
                continue
            for pseq in self.altpept[peptlen]:
                if pseq not in self.refpept[peptlen]:
                    pept[pseq] = self.altpept[peptlen][pseq]
                else:
                    if self.altpept[peptlen][pseq] > self.refpept[peptlen][pseq]:
                        pept[pseq] = self.altpept[peptlen][pseq] - self.refpept[peptlen][pseq]
        return pept
# end of PeptComparator
