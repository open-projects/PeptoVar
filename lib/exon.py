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

import weakref
from lib.graph import Node, Vertebra, Fibro

DEBUG = 0

class Exon:
    def __init__(self, transcript, beg, end, seq):
        self.chrom = transcript.chrom
        self.strand = transcript.strand
        self.beg = int(beg)
        self.end = int(end)
        self.seq = seq
        self._transcript = weakref.proxy(transcript)
        self._vertebras = self._makeBackbone()
        self._dirty_references = []
        
        if not len(seq):
            raise ValueError("no exon sequence - {}".format(' '.join((chrom, strand, str(beg), str(end)))))
        elif len(seq) != self.end - self.beg + 1:
            raise ValueError("wrong exon sequence length - {}".format(' '.join((chrom, strand, str(beg), str(end)))))
    
    def _makeBackbone(self):
        vertebras = [] # exon backbone
        samples = self._transcript.getSamples()
        for i, nucl in enumerate(self.seq):
            pos = self.beg + i
            vertebra = Vertebra(pos, nucl, samples)
            if len(vertebras):
                fibro = Fibro(pos - 1, samples)
                if self.strand == '+':
                    vertebras[-1].setNext(fibro)
                    fibro.setNext(vertebra)
                    vertebra.setBefore(fibro)
                    fibro.setBefore(vertebras[-1])
                else:
                    vertebra.setNext(fibro)
                    fibro.setNext(vertebras[-1])
                    vertebras[-1].setBefore(fibro)
                    fibro.setBefore(vertebra)
            vertebras.append(vertebra)
        return vertebras
    
    def _snp_validation(self, snp):
        if snp.beg < self.beg and snp.end < self.beg or self.end < snp.beg and self.end < snp.end:
            if DEBUG == 1:
                print("WARNING: variation is out of the exon - {} ...skipped".format(snp.id))
            return (None, None)
        if snp.beg < self.beg and self.end < snp.end:
            raise ValueError("overlaps with exon border")
        if snp.beg < self.beg and self.beg <= snp.end or snp.beg <= self.end and self.end < snp.end:
            if not (snp.isInDel() and not snp.isComplex() and snp.beg == snp.end + 1):
                raise ValueError("overlaps with exon border")
        
        if self.strand != snp.strand:
            snp.setComplement()
        snp_start = snp.beg - self.beg
        snp_stop = snp.end - self.beg
        if not (snp_start < 0 or len(self._vertebras) < snp_stop + 1):
            snp_ref = self.seq[snp_start : snp_stop + 1]
            if snp.ref.seq != '-' and snp_ref.upper() != snp.ref.seq.upper():
                raise ValueError("wrong reference allele: '{}' VS '{}'".format(snp.ref.seq, snp_ref))
                return (None, None)
        
        if self.strand == '-':
            (snp_start, snp_stop) = (snp_stop, snp_start)
        return (snp_start, snp_stop)
    
    def getFirstVertebra(self):
        if len(self._vertebras):
            return self._vertebras[0]
        return None
    
    def getLastVertebra(self):
        if len(self._vertebras):
            return self._vertebras[-1]
        return None

    def modify(self, snp):
        vertebras = self._vertebras
        (snp_exbeg, snp_exend) = self._snp_validation(snp)
        if snp_exbeg is not None and snp_exend is not None:
            if snp.ref.seq == '-':
                if snp_exbeg <= snp_exend and self.strand == '+' or snp_exbeg >= snp_exend and self.strand == '-':
                    beg_host = vertebras[snp_exbeg].getBefore()
                    end_host = vertebras[snp_exend]
                else:
                    if len(vertebras) < snp_exbeg + 1 and self.strand == '+' or snp_exbeg < 0 and self.strand == '-':
                        beg_host = vertebras[snp_exend].getNext()
                    else:
                        beg_host = vertebras[snp_exbeg].getBefore()
                    if snp_exend < 0 and self.strand == '+' or len(vertebras) < snp_exend + 1 and self.strand == '-':
                        end_host = vertebras[snp_exbeg].getBefore()
                    else:
                        end_host = vertebras[snp_exend].getNext()
            else:
                beg_host = vertebras[snp_exbeg]
                end_host = vertebras[snp_exend]
            
            self._transcript.appendVariation(snp, beg_host, end_host) # save variation backbone to find nonsynonymous SNP
            
            no_reference = True
            for allele in snp.getSampleAlleles():
                allele_seq = list(reversed(allele.seq)) if self.strand == '-' else allele.seq
                if allele.isReference(): # path for the reference allele
                    if allele.seq == '-':
                        ref_node = beg_host.getRefNode()
                        ref_node.appendAllele(allele)
                        beg_host.addNodeToGraph(ref_node)
                    else:
                        for allele_pos, pos in enumerate(range(snp_exbeg, snp_exend - 1, -1) if self.strand == '-' else range(snp_exbeg, snp_exend + 1)):
                            vertebra = vertebras[pos]
                            ref_node = vertebra.getRefNode()
                            ref_node.appendAllele(allele)
                            vertebra.addNodeToGraph(ref_node)
                            if allele_seq[allele_pos] != ref_node.nucl:
                                raise ValueError("reference allele mismatch")
                            if allele_pos > 0:
                                fibro = vertebra.getBefore()
                                ref_node = fibro.getRefNode()
                                ref_node.appendAllele(allele)
                                fibro.addNodeToGraph(ref_node)
                    for sample in self._transcript.getSamples():
                        if not allele.getSample(sample):
                            self._dirty_references.append({'host': beg_host, 'sample': sample})
                        
                    no_reference = False
                else: # path for alternative alleles
                    start_node = None
                    end_node = None
                    
                    for pos, nucl in enumerate(allele_seq):
                        pos = snp.beg + pos
                        pos_overhead = 0
                        if snp.end < pos:
                            pos_overhead = pos - snp.end
                            pos = snp.end
                        new_node = Node(pos, pos_overhead, nucl, allele)
                        if end_node:
                            end_node.setNext([new_node])
                        else:
                            start_node = new_node
                        end_node = new_node
                    
                    beg_host.addNodeToGraph(start_node) # join the allele start to the graph
                    vertebra = end_host.getNext()
                    if vertebra:
                        end_node.setNext(vertebra.getNodes()) # join the allele end to the graph
            if no_reference:
                for sample in self._transcript.getSamples():
                    self._dirty_references.append({'host': beg_host, 'sample': sample})
        return None
    
    def cleanup(self): # remove harmful reference paths (arise with overlapped variations)
        for dirty_ref in self._dirty_references:
            host = dirty_ref['host']
            ref_node = host.getRefNode()
            ref_node.removeSample(dirty_ref['sample'])
            alt_nodes = []
# end of Exon
