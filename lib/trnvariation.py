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
from operator import itemgetter

DEBUG = 0

class TrnAllele:
    def __init__(self, allele, prefix, suffix, mtx, trn):
        self.id = allele.id
        self.trn_id = prefix.id + ':' + allele.id + ':' + (suffix.id if suffix.id != '-' else '') # may be prefix.seq is better?
        self.allele = allele
        self.prefix = prefix
        self.suffix = suffix
        self.context = prefix.id + ':' + (suffix.id if suffix.id != '-' else '') # may be prefix.seq is better?
        self.mtx = mtx
        self.trn = trn
        self.sample_names = set(sample.name for sample in allele.getSamples())
# end of TrnAllele

class TrnVariation:
    def __init__(self, snp, beg, end, trn_strand):
        self.id = snp.id
        self.ref_allele_id = None
        self.snp = snp
        self.beg_vertebra = beg
        self.end_vertebra = end
        self.strand = trn_strand
        self._alleles = {}
        self._alleles_in_context = {}
        self._nonsyn_alleles = {}
        self._allele_usage = {}
        self._allele_rank = {}
        self._sample_names = set()
    
    def appendTrnAllele(self, trn_allele):
        if not trn_allele.context in self._alleles_in_context:
            self._alleles_in_context[trn_allele.context] = {}
        if not trn_allele.trn in self._alleles_in_context[trn_allele.context]:
            self._alleles_in_context[trn_allele.context][trn_allele.trn] = []
        self._alleles_in_context[trn_allele.context][trn_allele.trn].append(trn_allele)
        if trn_allele.allele.isReference():
            self.ref_allele_id = trn_allele.id
        if not trn_allele.id in self._alleles:
            self._alleles[trn_allele.id] = []
        self._alleles[trn_allele.id].append(trn_allele)
        self._sample_names = self._sample_names.union(trn_allele.sample_names)
        return None
        
    def getTrnAlleles(self):
        return self._alleles
    
    def getNonSynAlleles(self):
        return self._nonsyn_alleles
    
    def joinSynonymAlleles(self):
        if DEBUG and self.id == 'rs1130929':
            bp=1
        
        # may be not necessary
        context_list = list(self._alleles_in_context.keys())
        context_list.sort(key=len)
        context2remove = set()
        for i in range(len(context_list)):
            if i + 1 < len(context_list):
                for j in range(i + 1, len(context_list)):
                    if re.match(context_list[i], context_list[j]):
                        for trn in self._alleles_in_context[context_list[i]]:
                            if trn not in self._alleles_in_context[context_list[j]]:
                                self._alleles_in_context[context_list[j]][trn] = []
                            self._alleles_in_context[context_list[j]][trn].extend(self._alleles_in_context[context_list[i]][trn])
                        context2remove.add(context_list[i])
        for context in context2remove:
            if len(self._alleles_in_context[context]) < 2:
                self._alleles_in_context.pop(context)
        # may be not necessary
        
        n = 0
        for context in self._alleles_in_context:
            for trn in self._alleles_in_context[context]:
                n += 1
                l = len(self._alleles_in_context[context][trn])
                if l < 1:
                    raise ValueError("wrong allele context: {}".format(self.id))
                l = 1 / l
                for trn_allele in self._alleles_in_context[context][trn]:
                    if trn_allele.id in self._allele_usage:
                        self._allele_usage[trn_allele.id] += l
                    else:
                        self._allele_usage[trn_allele.id] = l
        if n:
            allele_ranks = []
            for allele_id in self._allele_usage:
                self._allele_usage[allele_id] /= n # allele usage calculation
                allele_ranks.append({'id': allele_id, 'usage': self._allele_usage[allele_id], 'ref': 1 if allele_id == self.ref_allele_id else 0})
            allele_ranks = sorted(allele_ranks, key=itemgetter('usage', 'ref', 'id'))
            for i in range(len(allele_ranks)):
                self._allele_rank[allele_ranks[i]['id']] = i
            
            for context in self._alleles_in_context:
                for trn in self._alleles_in_context[context]:
                    max_rank = -1
                    best_allele_id = None
                    
                    sample_rank = {}
                    sample_allele_id = {}
                    prefix_seq = {}
                    
                    for sample_name in self._sample_names:
                        if sample_name not in sample_rank:
                            sample_rank[sample_name] = -1
                            sample_allele_id[sample_name] = None
                            prefix_seq[sample_name] = None
                        for trn_allele in self._alleles_in_context[context][trn]:
                            if sample_name in trn_allele.sample_names and sample_rank[sample_name] < self._allele_rank[trn_allele.id]:
                                sample_rank[sample_name] = self._allele_rank[trn_allele.id]
                                sample_allele_id[sample_name] = trn_allele.id
                                prefix_seq[sample_name] = trn_allele.prefix.seq
                            if max_rank < self._allele_rank[trn_allele.id]:
                                max_rank = self._allele_rank[trn_allele.id]
                                best_allele_id = trn_allele.id
                        
                        if sample_allele_id[sample_name]:
                            for node in self.beg_vertebra.getNodes():
                                for allele in node.getAlleles():
                                    if allele.isPhased(sample_name) or allele.id == sample_allele_id[sample_name]:
                                        node.appendPrefixSeq(prefix_seq[sample_name])
                                        break
                    
                    if not best_allele_id:
                        raise ValueError("no best allele: {}".format(self.id))
                    if len(self._alleles_in_context[context]) > 1:
                        self._nonsyn_alleles[best_allele_id] = self._alleles[best_allele_id]
                        self._nonsyn_alleles[best_allele_id][0].allele.setNonSyn()
        return self._nonsyn_alleles
    
    def isNonSynonymous(self, allele_id = None):
        if allele_id:
            if allele_id in self._nonsyn_alleles:
                return True
            return False
        if len(self._nonsyn_alleles) > 1:
            return True
        return False
# end of TrnVariation
