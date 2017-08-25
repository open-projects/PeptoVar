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

from lib.exon import Exon
from lib.sample import Sample
from lib.graph import Node, Fibro, FShiftPath, CodonEmptyPrefix, CodonPrefix, CodonEmptySuffix
from lib.trnvariation import TrnAllele, TrnVariation
from lib.peptide import Peptide, PeptideDriver
from lib.seqtools import translate, complement, is_stopcodon

DEBUG = 0

class Transcript:
    def __init__(self, mrna_id , chrom = None, strand = None):
        self.id = mrna_id
        self.chrom = chrom
        self.strand = strand
        self._exon_number = 0
        self._exons = [] # to handle variations at the end of CDS the last exon must include the STOP CODON !!!!
        self._backbone_head = None # the head vertebra of the backbone
        self._graph_start = None # start nodes of the graph
        self._used_nodes = [] # array of used graph nodes
        self._variations = [] # array of variation linked to the backbones
        self._samples = [] # all samples in the run
        self._codon_start_nodes = [] # array of the first nodes in codons (for each sample separately)
    
    def cleanUsed(self):
        for used_node in self._used_nodes:
            used_node.cleanUsed() # graph clean up
        nodes = self._used_nodes
        self._used_nodes = []
        return nodes
    
    def setSamples(self, sample_names):
        for name in sample_names:
            if name == 'virtual':
                self._samples.append(Sample(name, 1, 1))
            else:
                self._samples.append(Sample(name, 1, 0))
                self._samples.append(Sample(name, 0, 1))
    
    def getSamples(self):
        return self._samples
    
    def appendExon(self, chrom, strand, beg, end, seq):
        beg = int(beg)
        end = int(end)
        if not self.chrom:
            self.chrom = chrom
        if not self.strand:
            self.strand = strand
        
        if not seq:
            raise ValueError("mRNA failure (no exon sequence) - {}".format(' '.join((self.id, self.chrom, self.strand, str(beg), str(end)))))
        
        if beg < 1 or end < 1:
            raise ValueError("mRNA failure (wrong exon location) - {}".format(' '.join((self.id, self.chrom, self.strand, str(beg), str(end)))))
        if end < beg:
            raise ValueError("mRNA failure (wrong exon orientation) - {}".format(' '.join((self.id, self.chrom, self.strand, str(beg), str(end)))))
        if len(seq) != end - beg + 1:
            raise ValueError("mRNA failure (wrong exon sequence length) - {}".format(' '.join((self.id, self.chrom, self.strand, str(beg), str(end)))))
        if len(self._exons):
            if beg <= self._exons[-1].end:
                raise ValueError("mRNA failure (wrong exon priority) - {}".format(' '.join((self.id, self.chrom, self.strand, str(beg), str(end)))))
        
        new_exon = Exon(self, beg, end, complement(seq) if strand == '-' else seq)
        
        
        fibro = Fibro(end, self._samples)
        vertebra1 = new_exon.getFirstVertebra()
        vertebra2 = new_exon.getLastVertebra()
        if self.strand == '-':
            self._backbone_head = fibro
            if not len(self._exons):
                vertebra1.setNext(Fibro(beg - 1, self._samples))
            else:
                last_vert = self._exons[-1].getLastVertebra()
                last_fibro = last_vert.getBefore()
                last_fibro.setBefore(vertebra1)
                vertebra1.setNext(last_fibro)
            vertebra2.setBefore(fibro)
            fibro.setNext(vertebra2)
        else:
            if not len(self._exons):
                self._backbone_head = Fibro(beg - 1, self._samples)
                self._backbone_head.setNext(vertebra1)
                vertebra1.setBefore(self._backbone_head)
            else:
                last_vert = self._exons[-1].getLastVertebra()
                last_fibro = last_vert.getNext()
                last_fibro.setNext(vertebra1)
                vertebra1.setBefore(last_fibro)
            vertebra2.setNext(fibro)
            fibro.setBefore(vertebra2)
        self._exons.append(new_exon)
        
        # check data consistency #
        vert = new_exon.getLastVertebra() if self.strand == '-' else  self._exons[0].getFirstVertebra()
        vert_num = 0
        while vert:
            vert_num += 1
            vert = vert.getNext()
            if vert:
                vert = vert.getNext()
        mtx_len = 0
        for exon in self._exons:
            mtx_len += exon.end - exon.beg + 1
        if vert_num != mtx_len:
            raise ValueError("faulty node sequence - {}".format(self.id))
        
        self._exon_number += 1
        return self._exon_number
    
    def getExons(self):
        return self._exons
    
    def getVariations(self):
        return self._variations
    
    def appendVariation(self, snp, beg_backbone_item, end_backbone_item):
        self._variations.append(TrnVariation(snp, beg_backbone_item, end_backbone_item, self.strand))
        return len(self._variations)
    
    def makeGraph(self):
        if not self._backbone_head:
            return None
        elif self._graph_start:
            return self._graph_start
        
        backbone_item = self._backbone_head
        self._graph_start = Node(0, 0, '')
        self._graph_start.setNext(backbone_item.getNodes())
        while backbone_item:
            backbone_item.addNodeToGraph()
            backbone_item = backbone_item.getNext()
        
        for sample in self._samples:
            self._findFrameShifts(sample)
            nodes = self._attachPrefixes(sample)
            self._codon_start_nodes.append((sample, nodes))
        return self._codon_start_nodes
    
    def _findFrameShifts(self, sample):
        def pathFinder(path):
            node = path['node']
            frame = path['frame']
            fshift_path = path['fsh_path']
            next_path = []
            if len(node.getNext(sample)):
                if node.nucl != '-':
                    frame += 1
                if frame > 2: frame = 0
                
                if node.isFrameShift():
                    fshift_path = fshift_path.clonePath()
                    if DEBUG:
                        print(fshift_path.id)
                    fsh_allele_id = []
                    for allele in node.getAlleles(): # don't use node.getAlleleID(): the node could belong to two and more alleles!
                        if allele.isFrameShift():
                            fsh_allele_id.append(allele.id)
                    fshift_path.appendAlleleID(fsh_allele_id)
                for next_node in node.getNext(sample):
                    if frame == 0:
                        if next_node.appendFShiftPath(fshift_path): # check if the path has not been used yet (STOP if it's used)
                            next_path.append({'node': next_node, 'frame': frame, 'fsh_path': fshift_path})
                    else:
                        next_path.append({'node': next_node, 'frame': frame, 'fsh_path': fshift_path})
            return next_path
        # end of pathFinder()
        
        tree = []
        fshift_startpath = FShiftPath(sample.id)
        for node in self._graph_start.getNext(sample):
            node.appendFShiftPath(fshift_startpath)
            tree.append({'node': node, 'frame': 0, 'fsh_path': fshift_startpath})
        while len(tree):
            new_tree = []
            for path in tree:
                new_path = pathFinder(path)
                new_tree.extend(new_path)
            tree = new_tree
    
    def _attachPrefixes(self, sample):
        def prefixDriver(node1): # node1 is the first nucleotide in the codon
            next_nodes = []
            if node1.isUsed() or not node1.getSample(sample):
                return next_nodes # STOP - the node have been used with other path or the node is not in the sample path
            
            node1.attachPrefix(CodonEmptyPrefix(sample, node1))
            vertebra1 = node1.getVertebra()
            if vertebra1:
                vertebra1.appendSample(sample)
            
            if node1.nucl == '-':
                next_nodes = node1.getNext(sample)
                return next_nodes
            
            nodes2 = node1.getNext(sample) # the second nucleotide
            while len(nodes2):
                new_nodes2 = []
                for node2 in nodes2:
                    node2.attachPrefix(CodonPrefix(sample, [node1]))
                    vertebra2 = node2.getVertebra()
                    if vertebra2:
                        vertebra2.appendSample(sample)
                    if node2.nucl == '-':
                        new_nodes2.extend(node2.getNext(sample))
                        continue
                    nodes3 = node2.getNext(sample)
                    while len(nodes3):
                        new_nodes3 = []
                        for node3 in nodes3: # the third nucleotide
                            node3.attachPrefix(CodonPrefix(sample, [node1, node2]))
                            vertebra3 = node3.getVertebra()
                            if vertebra3:
                                vertebra3.appendSample(sample)
                            if node3.nucl == '-':
                                new_nodes3.extend(node3.getNext(sample))
                                continue
                            
                            if DEBUG:
                                alleles1 = ",".join(allele.id for allele in node1.getAlleles())
                                alleles2 = ",".join(allele.id for allele in node2.getAlleles())
                                alleles3 = ",".join(allele.id for allele in node3.getAlleles())
                                print("{} pos:{}({})/{}/{} codon:{} alleles:{}/{}/{}".format(sample.id, node1.pos, node1.overpos, node2.pos, node3.pos, "".join((node1.nucl, node2.nucl , node3.nucl)), alleles1 if alleles1 else '-', alleles2 if alleles2 else '-', alleles3 if alleles3 else '-'))
                                if node1.pos == 152312604 or node2.pos == 152312604 or node3.pos == 152312604:
                                    bp=1
                            
                            if is_stopcodon(node1.nucl + node2.nucl + node3.nucl):
                                continue # there is no reason to attach prefixes to untranslated codons
                            next_nodes.extend(node3.getNext(sample))
                        nodes3 = new_nodes3
                nodes2 = new_nodes2
            node1.setUsed() # mark the node as `used`
            self._used_nodes.append(node1) # save used nodes to clean the marks before the next run
            
            return next_nodes
        
        tree = self._graph_start.getNext(sample)
        while len(tree):
            new_tree = []
            for node in tree:
                new_nodes = prefixDriver(node)
                new_tree.extend(new_nodes)
            tree = new_tree
        return self.cleanUsed() # array of the first nodes in codons can be used to make peptides 
    
    def translateVariations(self):
        for var in self._variations:
            beg_vertebra = var.beg_vertebra
            end_vertebra = var.end_vertebra.getNext()
            
            samples = []
            alleles = []
            alleles_id = set()
            for sample in self._samples:
                if beg_vertebra.getSample(sample) and len(beg_vertebra.getPrefixes(sample)): # take only translated alleles
                    for allele in var.snp.getSampleAlleles():
                        if allele.getSample(sample):
                            samples.append(sample)
                            if allele.id not in alleles_id:
                                alleles.append(allele)
                                alleles_id.add(allele.id)
                            
            
            for allele in alleles:
                if allele.seq == '-':
                    allele_seq = ''
                else:
                    allele_seq = allele.seq if self.strand == '+' else allele.seq[::-1]
                for sample in samples:
                    if end_vertebra:
                        siffixes = end_vertebra.getSuffixes(sample) # the last exon must have STOP CODON to get suffix from the last codon
                    else:
                        siffixes = [CodonEmptySuffix(sample)]
                    for prefix in beg_vertebra.getPrefixes(sample): # get allele translations for each pair prefix/suffix
                        for suffix in siffixes:
                            suffix_len = 3 - (len(prefix.seq) + len(allele_seq)) % 3
                            if suffix_len == 3: suffix_len = 0
                            old_suffix_id = suffix.id
                            suffix = suffix.getTrimmed(suffix_len)
                            suffix.id = old_suffix_id
                            mtx = prefix.seq + allele_seq + suffix.seq
                            if len(allele_seq) == 0 and len(mtx) == 0:
                                var.appendTrnAllele(TrnAllele(allele, prefix, suffix, '-', '-'))
                            elif len(mtx) < 3:
                                continue
                            else:
                                (trn, tail) = translate(mtx)
                                var.appendTrnAllele(TrnAllele(allele, prefix, suffix, mtx, trn))
            point_for_breakpoint = True
        return self._variations
    
    def getProteins(self, optimization = False):
        return self.getPeptides([0], optimization)
    
    def getPeptides(self, pept_len, optimization = False):
        pept_container = {}
        codon_start_nodes = []
        if pept_len and len(pept_len) > 0:
            if len(pept_len) == 1 and pept_len[0] == 0: # proteins
                for sample in self._samples:
                    codon_start_nodes.append((sample, self._graph_start.getNext(sample, '')))
            elif self._graph_start: # peptides
                codon_start_nodes = self._codon_start_nodes
        
        for sample, start_nodes in codon_start_nodes:
            for node in start_nodes:
                vertebra = node.getVertebra()
                if node.getPrefixSeq('') or not vertebra or vertebra and len(vertebra.getNodes()) == 1:
                    peptDriver = PeptideDriver(self.chrom, self.id, sample, pept_len, optimization)
                    peptDriver.appendNode(node)
                    peptDriver.setPeptContainer(pept_container)
                    driverArray = [peptDriver]
                    while driverArray:
                        new_driverArray = []
                        for driver in driverArray:
                            nodes = driver.getNext()
                            if len(nodes) > 1:
                                for node in nodes:
                                    new_peptDriver = driver.copy() # !!!
                                    if new_peptDriver.appendNode(node):
                                        new_driverArray.append(new_peptDriver)
                            elif len(nodes) == 1:
                                if driver.appendNode(nodes[0]):
                                    new_driverArray.append(driver)
                        driverArray = new_driverArray
        return pept_container
    
    def getBackboneProtein(self):
        mtx = ""
        vertebra = self._backbone_head
        while vertebra:
            nucl = vertebra.getRefNode().nucl
            if nucl != '-':
                mtx += nucl
            vertebra = vertebra.getNext()
        return translate(mtx)
    
    def joinSynonymPathes(self):
        nonsynonym_var = []
        for trn_var in self._variations:
            trn_var.joinSynonymAlleles()
            if trn_var.isNonSynonymous():
                nonsynonym_var.append(trn_var)
        return nonsynonym_var
# end of Transcript
