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

import operator

class Sample:
    def __init__(self, sample_name = None, allele_1 = 0, allele_2 = 0):
        self.id = None
        self.name = sample_name
        self.allele_1 = allele_1
        self.allele_2 = allele_2
        self._setID()
    
    def setAlleleNum(self, num):
        if num == 1:
            self.allele_1 = 1
        elif num == 2:
            self.allele_2 = 1
        else:
            raise ValueError("wrong sample allele number")
        self._setID()
    
    def _setID(self):
        self.id = '_'.join((self.name, str(self.allele_1), str(self.allele_2)))
    
    def __or__(self, other):
        if self.name != other.name:
            raise ValueError("wrong sample addition")
        return Sample(self.name, self.allele_1 | other.allele_1, self.allele_2 | other.allele_2)
    
    def __ior__(self, other):
        if self.name != other.name:
            raise ValueError("wrong sample addition")
        return Sample(self.name, self.allele_1 | other.allele_1, self.allele_2 | other.allele_2)
    
    def __eq__(self, other): # self == other
        if self.name == other.name and self.allele_1 == other.allele_1 and self.allele_2 == other.allele_2:
            return True
        return False
    
    def __le__(self, other): # self <= other
        if self.name == other.name and (self.allele_1 == other.allele_1 and self.allele_2 == other.allele_2 or other.allele_1 and other.allele_2 and (self.allele_1 == other.allele_1 or self.allele_2 == other.allele_2)):
            return True
        return False
    
    def __ge__(self, other): # self >= other
        if self.name == other.name and (self.allele_1 == other.allele_1 and self.allele_2 == other.allele_2 or self.allele_1 and self.allele_2 and (self.allele_1 == other.allele_1 or self.allele_2 == other.allele_2)):
            return True
        return False
    
    def __sub__(self, other): # self - other
        if self.name != other.name:
            raise ValueError("wrong sample addition")
        if other.allele_1 and other.allele_2:
            return Sample(self.name, 0, 0)
        if self.allele_1 >= other.allele_1 and self.allele_2 >= other.allele_2:
            return Sample(self.name, self.allele_1 - other.allele_1, self.allele_2 - other.allele_2)
        return Sample(self.name, self.allele_1, self.allele_2)
        
# end of Sample

class SampleContainer:
    def __init__(self):
        self._samples = {}
    
    def removeSample(self, sample):
        if sample.name in self._samples:
            self._samples[sample.name] -= sample
    
    def removeSamples(self, samples):
        for sample in samples:
            self.removeSample(sample)
    
    def appendSample(self, sample):
        if sample.name in self._samples:
            self._samples[sample.name] |= sample
        else:
            self._samples[sample.name] = sample
        return self._samples[sample.name]
    
    def appendSamples(self, samples):
        for sample in samples:
            self.appendSample(sample)
    
    def getSample(self, sample = None):
        if sample:
            if sample.name in self._samples and self._samples[sample.name] >= sample:
                return sample
        else:
            return list(self._samples.values())
        return None
# end of SampleContainer
