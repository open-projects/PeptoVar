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
from lib.seqtools import complement
from lib.sample import SampleContainer, Sample

class Allele:
    def __init__(self, snp_id, seq, strand, freq, ref = None):
        self.snp_id = snp_id
        self.seq = seq.upper()
        self.length = len(seq)
        self.freq = freq
        self.strand = strand
        self._is_ref = True if ref else False
        self.id = self._get_id()
        self._is_nonsyn = False
        self._indel = 0
        self._fshift = False
        self._samples = SampleContainer()
        self._phased = set() # samples for which the allele is phased
    
    def _get_id(self):
        return ":".join((self.snp_id.lower(), self.seq.upper())) + ("=" + str(self.freq) if self.freq else "") + ("(ref)" if self._is_ref else "(alt)")
    
    def rebuildFormat(self): # convert VCF format to dbSNP format
        self.seq = self.seq[1:]
        if len(self.seq) == 0:
            self.seq = '-'
        self.length -= 1
        self.id = self._get_id()
    
    def setFreq(self, freq):
        self.freq = freq
        self.id = self._get_id()
    
    def addSample(self, sample, phased = False):
        if phased:
            self._phased.add(sample.name)
        return self._samples.appendSample(sample)
    
    def getSamples(self):
        return self._samples.getSample()
    
    def getSample(self, sample):
        return self._samples.getSample(sample)
    
    def isPhased(self, sample_name):
        if sample_name in self._phased:
            return True
        return False
    
    def isReference(self):
        return self._is_ref
    
    def setNonSyn(self):
        self._is_nonsyn = True
    
    def isNonSyn(self):
        return self._is_nonsyn
    
    def setInDel(self, indel):
        self._indel = indel
    
    def isInDel(self):
        return self._indel
    
    def setFrameShift(self):
        self._fshift = True
    
    def isFrameShift(self):
        return self._fshift
        
# end of Allele

class SNP:
    def __init__(self, snp, samples, optimization):
        if hasattr(snp, 'info'): # normalization
            if 'MATCHED_REV' in snp.info:
                if snp.info['MATCHED_REV'] == True:
                    new_alt = []
                    for allele in snp.alleles:
                        new_alt.append(complement(allele)[::-1])
                    snp.alleles = tuple(new_alt)
        
        self.id = snp.id
        self.beg = snp.start + 1
        self.end = snp.stop
        self.strand = '+'
        self.ref = Allele(snp.id, snp.ref, self.strand, 0, 'reference')
        self.alleles = [self.ref] # all alleles in the vasnp.id == 'rs529174822;rs563866452'snp.id == 'rs529174822;rs563866452'riation (self.alleles[0] - the reference allele forever!!!)
        self._sample_allele_index = set() # self.allele indexes for samples (they will be taken for graph construction)
        self._fshift = False
        self._indel = False
        self._complex = True if re.search("[,; ]", snp.id) else False
        self.err = []
        if not optimization:
            self.ref.setNonSyn()
        
        ref_af = -1
        if hasattr(snp, 'info'):
            if 'AF' in snp.info: # according to VCF file description allele frequency (AF) specified ONLY for alternative alleles
                if len(snp.alts) == len(snp.info['AF']):
                    ref_af = 1
                    alt_alleles = []
                    for i, allele in enumerate(snp.alts):
                        if allele == snp.ref: continue
                        af = round(float(snp.info['AF'][i]), 4)
                        ref_af -= af # AF for reference allele should be calculated
                        new_allele = Allele(snp.id, allele, self.strand, af)
                        if not optimization:
                            new_allele.setNonSyn()
                        indel = len(allele) - len(self.ref.seq) # indel>0 - insertion, indel<0 - deletion
                        if indel:
                            self._indel = True
                        frameshift = indel % 3
                        if frameshift:
                            self._fshift = True
                            new_allele.setFrameShift()
                            self.ref.setFrameShift()
                        new_allele.setInDel(indel)
                        alt_alleles.append(new_allele)
                    if ref_af < 0:
                        self.err.append("wrong allele frequency for reference allele")
                    else:
                        self.ref.setFreq(round(ref_af, 4))
                        self.alleles.extend(alt_alleles)
                else:
                    self.err.append("wrong 'AF' data")
        if ref_af < 0:
            # set AF = 0 for all alleles
            for allele in snp.alts:
                # self.ref has been initialized by '0'
                new_allele = Allele(snp.id, allele, self.strand, 0)
                if not optimization:
                    new_allele.setNonSyn()
                indel = len(allele) - len(self.ref.seq) # indel>0 - insertion, indel<0 - deletion
                if indel:
                    self._indel = True
                frameshift = indel % 3
                if frameshift:
                    self._fshift = True
                    new_allele.setFrameShift()
                new_allele.setInDel(indel)
                self.alleles.append(new_allele)
        
        for sample in samples:
            if sample == 'virtual': # the virtual sample (all alleles are unphased)
                for i, allele in enumerate(self.alleles):
                    allele.addSample(Sample(sample, 1, 1))
                    self._sample_allele_index.add(i)
            else:
                if hasattr(snp, 'samples'):
                    try:
                        sample_data = snp.samples[sample]
                    except:
                        self.err.append("no data for sample {}".format(sample))
                        continue
                    else:
                        if hasattr(sample_data, 'alleles') and len(sample_data.alleles) > 0:
                            if len(snp.samples[sample].alleles) > 1:
                                allele1_idx = self._findAlleleIndex(sample_data.alleles[0])
                                allele2_idx = self._findAlleleIndex(sample_data.alleles[1])
                                
                                if allele1_idx is None and allele2_idx is None:
                                    allele1_idx = self._findAlleleIndex(snp.ref)
                                    allele2_idx = allele1_idx
                                    self.err.append("no alleles for sample {} (replaced by reference)".format(sample))
                                elif allele2_idx is None:
                                    allele2_idx = allele1_idx
                                    self.err.append("no allele2 for sample {} (replaced by allele1)".format(sample))
                                elif allele1_idx is None:
                                    allele1_idx = allele2_idx
                                    self.err.append("no allele1 for sample {} (replaced by allele2)".format(sample))
                                    
                                if hasattr(sample_data, 'phased') and sample_data.phased:
                                    self.alleles[allele1_idx].addSample(Sample(sample, 1, 0), True) # True is phased
                                    self.alleles[allele2_idx].addSample(Sample(sample, 0, 1), True) # True is phased
                                else:
                                    self.alleles[allele1_idx].addSample(Sample(sample, 1, 1))
                                    self.alleles[allele2_idx].addSample(Sample(sample, 1, 1))
                                self._sample_allele_index.update([allele1_idx, allele2_idx])
                            else:
                                allele1_idx = self._findAlleleIndex(sample_data.alleles[0])
                                if allele1_idx is None:
                                    allele1_idx = self._findAlleleIndex(snp.ref)
                                    self.err.append("no allele1 for sample {} (replaced by reference)".format(sample))
                                self.alleles[allele1_idx].addSample(Sample(sample, 1, 1))
                                self._sample_allele_index.update([allele1_idx])
                        else:
                            raise ValueError("no alleles for sample {}".format(sample))
                else:
                    raise ValueError("no samples in VCF file")
        
        if self._indel and not self._complex:
            for allele in self.alleles:
                allele.rebuildFormat()
            self.beg += 1
    
    def _findAlleleIndex(self, seq): # find the allele with a given sequence
        for i, allele in enumerate(self.alleles):
            if allele.seq == seq:
                return i
        return None
    
    def setComplement(self):
        self.strand = '-' if self.strand == '+' else '+'
        for allele in self.alleles:
            allele.seq = complement(allele.seq)
            allele.strand = self.strand
    
    def filterAF(self, min_af):
        if min_af > 0:
            self._sample_allele_index = set(i for i in self._sample_allele_index if self.alleles[i].freq >= min_af)
            if len(self._sample_allele_index) > 1:
                return True
            elif len(self._sample_allele_index) == 1 and 0 not in self._sample_allele_index: # reference allele has index = 0
                return True
            else:
                return False
        return True
    
    def getSampleAlleles(self, sample = None):
        alleles = []
        for i in self._sample_allele_index:
            if sample:
                if self.alleles[i].getSample(sample):
                    alleles.append(self.alleles[i])
            else:
                alleles.append(self.alleles[i])
        return alleles
    
    def isInDel(self):
        return self._indel
    
    def isComplex(self):
        return self._complex
# end of SNP
