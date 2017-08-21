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

import weakref
from lib.seqtools import translate
from lib.sample import SampleContainer

class Node:
    def __init__(self, pos, overpos, nucl, allele = None):
        self.pos = pos
        self.overpos = overpos
        self.nucl = nucl
        self._is_ref = False if allele else True # TRUE if it belongs to the reference sequence
        self._next = [] # the next nodes
        self._vertebra = None # reference to the master vertebra
        self._allele_id = allele.id if allele else 'ref'
        self._alleles = [allele] if allele else [] # alleles of variations
        self._is_fshift = True if allele and allele.isInDel() else False # True if the node belongs to a frameshift
        self._fshift_path = {}
        self._is_used = False # True if the node has been used
        self._prefix_seq = set()
        self._samples = SampleContainer()
        self._in_graph = False
        if allele:
            self.appendSamples(allele.getSamples())
    
    def setInGraph(self):
        self._in_graph = True
    
    def isInGraph(self):
        return self._in_graph
    
    def setVertebra(self, vertebra):
        self._vertebra = weakref.proxy(vertebra)
    
    def getVertebra(self):
        return self._vertebra
    
    def appendSamples(self, samples):
        self._samples.appendSamples(samples)
    
    def removeSample(self, sample):
        self._samples.removeSample(sample)
    
    def appendPrefixSeq(self, prefix_seq):
        self._prefix_seq.add(prefix_seq.lower())
        return self._prefix_seq
    
    def appendAllele(self, allele):
        if not self._is_ref: # an alternative (not reference) allele has only one allele_id
            raise ValueError("wrong allele ID for node in pos: {}".format(self.pos))
        self._alleles.append(allele)
        self.appendSamples(allele.getSamples())
        if allele.isFrameShift():
            self._is_fshift = True
        if self._allele_id == 'ref':
            self._allele_id = allele.id
        else:
            self._allele_id += ',' + allele.id
        return len(self._alleles)
    
    def getAlleleID(self):
        return self._allele_id
    
    def getAlleles(self, sample = None):
        if sample:
            alleles = []
            for allele in self._alleles:
                if allele.getSample(sample):
                    alleles.append(allele)
            return alleles
        return self._alleles
    
    def isReference(self):
        return self._is_ref
    
    def getSample(self, sample = None):
        return self._samples.getSample(sample)
    
    def getPrefixSeq(self, prefix_seq):
        if prefix_seq.lower() in self._prefix_seq:
            return self._prefix_seq
        return None
    
    def getNext(self, sample = None, pefix_seq = None):
        if sample:
            next_nodes = []
            for next_node in self._next:
                if next_node.getSample(sample):
                    if pefix_seq != None and len(self._next) > 1: # check prefix if branches exist
                        if next_node.getPrefixSeq(pefix_seq):
                            next_nodes.append(next_node)
                    else:
                        next_nodes.append(next_node)
            return next_nodes
        return self._next
    
    def setNext(self, next_nodes):
        self._next = next_nodes
        return self._next
    
    def setUsed(self):
        self._is_used = True
        
    def cleanUsed(self):
        self._is_used = False
    
    def isUsed(self):
        return self._is_used
    
    def appendFShiftPath(self, fshift_path):
        if fshift_path.sample_id not in self._fshift_path:
            self._fshift_path[fshift_path.sample_id] = set()
        if fshift_path.id not in self._fshift_path[fshift_path.sample_id]:
            self._fshift_path[fshift_path.sample_id].add(fshift_path.id)
            return True
        return False
    
    def getFShiftPathSet(self, sample):
        if sample.id in self._fshift_path:
            return self._fshift_path[sample.id]
        return set()
    
    def isFrameShift(self):
        return self._is_fshift
    
    def attachPrefix(self, prefix):
        if self._vertebra:
            self._vertebra.appendPrefix(prefix)
            return True
        return False
    
    def getPos(self):
        return {'pos': self.pos, 'overpos': self.overpos, 'id': self._allele_id}
# end of Node

class Vertebra:
    def __init__(self, pos, nucl, samples):
        self.pos = pos
        self._ref_node = Node(pos, 0, nucl)
        self._nodes = [] # nodes in the graph
        self._next = None # the next item in the backbone
        self._before = None # the previous item in the backbone
        self._prefixes = uContainer() # prefixes (CodonPrefix objects) to build the first codons for variation alleles
        self._samples = SampleContainer()
        self._alt_allele_samples = SampleContainer()
        self._all_samples = samples
    
    def appendSample(self, sample):
        self._samples.appendSample(sample)
        
    def getSample(self, sample):
        return self._samples.getSample(sample)
    
    def appendPrefix(self, prefix):
        return self._prefixes.append(prefix)
    
    def getPrefixes(self, sample):
        prefixes = self._prefixes.get(sample)
        if not len(prefixes):
            raise ValueError("no prefix in position: {}".format(self.pos))
        return prefixes
    
    def getSuffixes(self, sample):
        suffixes = []
        nodes1 = self._nodes
        while len(nodes1):
            new_nodes1 = []
            empty_suffix1 = False
            for node1 in nodes1:
                if not node1.getSample(sample):
                    continue
                if node1.nucl == '-':
                    next_nodes1 = node1.getNext(sample)
                    if len(next_nodes1):
                        new_nodes1.extend(next_nodes1)
                    else:
                        if not empty_suffix1:
                            suffixes.append(CodonSuffix(sample, [Node(0, 0, 'N'), Node(0, 0, 'N')]))
                            empty_suffix1 = True
                    continue
                nodes2 = node1.getNext(sample)
                while len(nodes2):
                    new_nodes2 = []
                    empty_suffix2 = False
                    for node2 in nodes2:
                        if node2.nucl == '-':
                            next_nodes2 = node2.getNext(sample)
                            if len(next_nodes2):
                                new_nodes2.extend(next_nodes2)
                            else:
                                if not empty_suffix2:
                                    suffixes.append(CodonSuffix(sample, [node1, Node(0, 0, 'N')]))
                                    empty_suffix2 = True
                            continue
                        suffixes.append(CodonSuffix(sample, [node1, node2]))
                    nodes2 = new_nodes2
            nodes1 = new_nodes1
        if not len(suffixes):
            suffixes.append(CodonSuffix(sample, [Node(0, 0, 'N'), Node(0, 0, 'N')]))
        return suffixes
    
    def addNodeToGraph(self, node = None): # add alleles to the graph
        if node: # add the variation to the graph
            if not node.isInGraph():
                if node.isReference():
                    if self._next:
                        node.setNext(self._next.getNodes())
                elif self._all_samples[0].name != 'virtual':
                    for node_sample in node.getSample():
                        if self._alt_allele_samples.getSample(node_sample):
                            raise ValueError("ambiguous variation")
                        else:
                            self._alt_allele_samples.appendSample(node_sample)
                    
                node.setInGraph()
                node.setVertebra(self)
                self._nodes.append(node)
        else: # add the reference allele to the graph if no alternatives
            node = self._ref_node
            if not node.isInGraph() and not len(self._nodes): # if reference not in the graph
                if self._next:
                    node.setNext(self._next.getNodes())
                node.setInGraph()
                node.setVertebra(self)
                node.appendSamples(self._all_samples)
                self._nodes.append(node)
        return None
    
    def getRefNode(self):
        return self._ref_node
    
    def getNodes(self):
        return self._nodes
    
    def setNext(self, item):
        self._next = item
    
    def setBefore(self, item):
        self._before = item
    
    def getNext(self):
        return self._next
    
    def getBefore(self):
        return self._before
# end of Vertebra

class Fibro(Vertebra):
    def __init__(self, pos, samples):
        super().__init__(pos + 0.5, '-', samples)
# end of Fibro    

class Codon:
    def __init__(self, nodes):
        self.nodes = nodes
        self.codon = nodes[0].nucl + nodes[1].nucl + nodes[2].nuc
        self.aa = translate(self.codon)[0]
# end of Codon

class CodonPrefix:
    def __init__(self, sample, nodes = []):
        self.id = "-"
        self.seq = ""
        self.sample = sample
        self._fsh_path_set = set()
        self.nodes = nodes
        if len(nodes):
            self.id = "-".join(str(id(node)) for node in nodes)
            self.seq = "".join(node.nucl.lower() for node in nodes)
            self._setFShiftPathSet(self.nodes[0])
        
    def _setFShiftPathSet(self, node):
        self._fsh_path_set = node.getFShiftPathSet(self.sample)
        
    def getAllelesID(self):
        ids = []
        for node in self.nodes:
            for allele in node.getAlleles():
                if not allele.isNonSyn():
                    allele_id = '[' + allele.id + ']'
                else:
                    allele_id = allele.id
                if allele_id not in ids:
                    ids.append(allele_id)
        return ids
        
    def getFShiftPathSet(self):
        return self._fsh_path_set
# end of CodonPrefix

class CodonEmptyPrefix(CodonPrefix):
    def __init__(self, sample, node):
        super().__init__(sample)
        super()._setFShiftPathSet(node)
# end of CodonEmptyPrefix

class CodonSuffix (CodonPrefix):
    def __init__(self, sample, nodes = []):
        super().__init__(sample, nodes)
    
    def getTrimmed(self, length = 0):
        return CodonSuffix(self.sample, self.nodes[:length])
# end of CodonSuffix

class CodonEmptySuffix(CodonSuffix):
    def __init__(self, sample):
        super().__init__(sample, [Node(0, 0, 'N'), Node(0, 0, 'N')])
# end of CodonEmptyPrefix

class uContainer:
    def __init__(self):
        self._samples = {}
        self._ids = set()
    
    def append(self, item):
        container_id = item.sample.id + item.id
        if container_id not in self._ids:
            if item.sample.id not in self._samples:
                self._samples[item.sample.id] = []
            self._samples[item.sample.id].append(item)
            self._ids.add(container_id)
            return True
        return False
    
    def get(self, sample):
        if sample.id in self._samples:
            return self._samples[sample.id]
        return []
# end of uContainer

class FShiftPath:
    def __init__(self, sample_id, path_id = '-'):
        self.id = path_id
        self.sample_id = sample_id
    
    def appendAlleleID(self, allele_id_set):
        set_id = ",".join(allele_id_set)
        if self.id != "-":
            if not self.id.endswith(set_id):
                self.id += "," + set_id
                return True
            return False
        else:
            self.id = set_id
            return True
    
    def clonePath(self):
        return FShiftPath(self.sample_id, self.id)
    
    def checkConsistency(self):
        checker = set()
        for allele_id in self.id.split(","):
            if allele_id in checker:
                return allele_id
            else:
                checker.add(allele_id)
        return None
# end of FShiftPath
