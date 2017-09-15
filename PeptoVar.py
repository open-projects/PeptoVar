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

import os
import errno
import re
import argparse
import time
import datetime
import shutil

try:
    from pysam import VariantFile, VariantHeader
except:
    cprint.printFail("\nCan't parse VCF file. Please, install pysam module for python3\n")
    exit()
from lib.gff import Gff
from lib.transcript import Transcript
from lib.output import OutFileContainer, ColorPrint
from lib.input import InFile
from lib.variation import SNP
from lib.seqtools import PeptComparator
from lib.comparator import UniPep

DEBUG = 0

if DEBUG:
    from lib.seqmod import Sequence

def main():
    start_time = time.time()
    total_tr = 0
    cprint = ColorPrint()
    
    input_parser = argparse.ArgumentParser(description='PeptoVar - Peptides on Variations: the program for personalization of protein coding genes and population-wide peptidome generation.')
    input_parser.add_argument('-gff', metavar='file.gff', default=None, help='GFF input file', required=False)
    input_parser.add_argument('-fasta', metavar='file.fasta', default=None, help='FASTA input file', required=False)
    input_parser.add_argument('-vcf', metavar='file.vcf.gz', default=None, help='bgzip-compressed VCF input file (need an index file)', required=False)
    input_parser.add_argument('-tmpdir', metavar='dirpath', default=None, help='TEMP directory', required=False)
    input_parser.add_argument('-samples', metavar='name', nargs='+', default=list(), help='a sample name or a pair of sample names in VCF file; for two samples (donor/recipient) only unique peptides will be represented)', required=False)
    input_parser.add_argument('-tagaf', metavar='TAG_AF', default='AF', help='allele frequency tag in VCF file (for example: EUR_AF, SAS_AF, AMR_AF etc.); use with `-minaf` argument, default=AF', required=False)
    input_parser.add_argument('-minaf', metavar='THRESHOLD', type=float, default=0, help='allele frequency (AF) threshold; alleles with AF < THRESHOLD will be ignored (AF=0 will be set for alleles with no data)', required=False)
    input_parser.add_argument('-var', metavar='all | used', choices=['all', 'used'], help='save translated polymorphisms (all or only used to make peptides)', required=False)
    input_parser.add_argument('-nopt', action='store_false', default=True, help='do not use optimization (may cause high CPU load and memory usage)')
    input_parser.add_argument('-peptlen', metavar='LENGTH', nargs='+', type=int, default=list(), help='lengths of peptides (0 - full-length proteins)', required=False)
    input_parser.add_argument('-outdir', metavar='dirpath', default='./output', help='output directory (will be created if not exists, default=./output)', required=False)
    input_parser.add_argument('-indir', metavar='dirpath', default=None, help='input directory for files *.vcf.gz, *.vcf.gz.tbi, *.gff and *.fasta - if no sequences in GFF file; the files MUST have the same name for each locus (chromosome)', required=False)
    input_parser.add_argument('-trnlist', metavar='transcriptID', nargs='+', default=list(), help='list of transcriptID for processing', required=False)
    input_parser.add_argument('-trnfile', metavar='transcriptID.txt', default=None, help='one column text file with the transcriptID list for processing', required=False)
    if DEBUG:
        input_parser.add_argument('-seq', metavar='file.data', default=None, help='DATA input file', required=False)
    
    args = input_parser.parse_args()
    gff_file = args.gff
    fasta_file = args.fasta
    seq_file = args.seq if DEBUG else None
    vcf_file = args.vcf
    tmp_dir = args.tmpdir
    save_var = args.var
    optimization = args.nopt
    min_af = args.minaf
    tag_af = args.tagaf
    trnlist = args.trnlist
    trnfile = args.trnfile
    pept_len = args.peptlen
    pept_len.sort()
    do_prot = False
    if len(pept_len) and pept_len[0] == 0:
        do_prot = True
        pept_len.pop(0)
    outdir = args.outdir
    indir = args.indir
    samples = args.samples
    
    peptdb = None
    protdb = None
    
    if not min_af:
        cprint.printWarning("\nLow value of -minaf argument can cause high memory usage and increasing computational time!\n")
    
    tmpdir_created = 0
    if tmp_dir:
        if not os.path.exists(tmp_dir):
            try:
                os.makedirs(tmp_dir)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    cprint.printFail("\nCan't create a temporary directory!\n")
                    exit()
            tmpdir_created = 1
        os.environ["TMPDIR"] = tmp_dir
        os.environ["SQLITE_TMPDIR"] = tmp_dir
    
    for plength in pept_len:
        if plength < 0:
            cprint.printFail("\nWrong -peptlen value (must be > 0)!\n")
            exit()
    
    if min_af < 0 or 1 < min_af:
        cprint.printFail("\nWrong -minaf value (must be between 0 and 1)!\n")
        exit()
    
    try:
        outfiles = OutFileContainer(outdir)
    except:
        cprint.printFail("\nCan't create output files!\n")
        exit()
    
    if not(vcf_file or indir):
        cprint.printFail("\nYou have no input data!\n")
        input_parser.print_help()
        exit()
    
    transcript_set = {}
    transcript_num = 0
    if trnfile:
        with open(trnfile, "r") as trnlistfile:
            for line in trnlistfile:
                trnid = line.strip()
                if len(trnid) > 0:
                    transcript_set[trnid] = 0
                    transcript_num += 1
        trnlistfile.close()
    for trnid in trnlist:
        transcript_set[trnid] = 0
        transcript_num += 1
    
    if len(samples) == 0:
        samples.append('virtual')
        cprint.printWarning('MODE: population')
        peptdb = UniPep(samples[0], '-', tmp_dir)
        protdb = UniPep(samples[0], '-', tmp_dir)
    elif len(samples) == 1:
        cprint.printWarning('MODE: sample {}'.format(samples[0]))
        peptdb = UniPep(samples[0], '-', tmp_dir)
        protdb = UniPep(samples[0], '-', tmp_dir)
    elif len(samples) == 2:
        cprint.printWarning('MODE: transplantation')
        peptdb = UniPep(samples[0], samples[1], tmp_dir)
        protdb = UniPep(samples[0], samples[1], tmp_dir)
    else:
        cprint.printFail("\nCan't take more than two samples\n")
        exit()
    
    input_bulk = []
    if vcf_file:
        if gff_file:
            input_bulk.append({'name': 'input', 'vcf': vcf_file, 'gff': gff_file, 'fasta': fasta_file, 'seq': seq_file})
        else:
            vcf_file = re.sub('.*\/', '', vcf_file)
            cprint.printFail("No GFF file for VCF {} ...skipped!".format(vcf_file))
            outfiles.writeWarning([vcf_file, "no GFF file for VCF", "skipped"])
    if indir:
        infiles = InFile(indir)
        input_bulk.extend(infiles.bulk)
    
    if DEBUG:
        errlog = open("ERROR.log", 'w')
    
    all_done = 0
    for fileset in input_bulk:
        print("{} files parsing...\n".format(fileset['name']))
        
        try:
            vcf = VariantFile(fileset['vcf'])
        except:
            cprint.printWarning("{} - wrong VCF file ...skipped".format(fileset['vcf']))
            outfiles.writeWarning([fileset['gff'], "wrong VCF file", "skipped"])
            continue
        
        try:
            vcf.fetch('test', 1, 2) # just to check vcf health
        except:
            cprint.printWarning("{} - no index file (use tabix) ...skipped".format(fileset['vcf']))
            outfiles.writeWarning([fileset['vcf'], "no index file (use tabix)", "skipped"])
            continue
        
        try:
            gff = Gff(fileset['gff'], tmp_dir)
        except ValueError as err:
            cprint.printWarning("{} - {} ...skipped".format(fileset['gff'], err.args[0]))
            outfiles.writeWarning([fileset['gff'], err.args[0], "skipped"])
            continue
        try:
            if fileset['seq']:
                gff.attachSeq(fileset['seq'])
            else:
                gff.attachSeq(fileset['fasta'])
        except ValueError as err:
            cprint.printWarning("Can't attach genome sequences: {} ...skipped".format(err.args[0]))
            outfiles.writeWarning([fileset['gff'], err.args[0], "skipped"])
            continue
        
        miss_samples = set() # just to catch female
        for sample in samples:
            if hasattr(vcf, 'header') and sample != 'virtual' and sample not in vcf.header.samples:
                miss_samples.add(sample)
        
        for trn_id in gff.getTranscriptsID():
            try:
                if len(transcript_set):
                    if trn_id not in transcript_set:
                        continue
                    transcript_set[trn_id] = 1
                    transcript_num -= 1
                
                total_tr += 1
                transcript = Transcript(trn_id)
                transcript.setSamples(samples)
                
                debugseq = None
                n_exons = 0
                for exon in gff.getTranscriptExons(trn_id):
                    try:
                        transcript.appendExon(exon.chrom, exon.strand, exon.beg, exon.end, exon.seq)
                    except ValueError as err:
                        cprint.printWarning("{} - {} ...skipped!".format(trn_id, err.args[0]))
                        outfiles.writeWarning([trn_id, err.args[0], "skipped"])
                        n_exons = 0
                        break
                    else:
                        n_exons += 1
                    
                    if DEBUG:
                        if not debugseq:
                            debugseq = Sequence(transcript.chrom, transcript.id, transcript.getSamples())
                        debugseq.append(exon.strand, exon.beg, exon.end, exon.seq)
                if n_exons == 0:
                    continue
                cprint.printOk("IN PROCESS #{}: {} (locus {})".format(total_tr, transcript.id, transcript.chrom))
                
                print("transript modification...") #print('\x1b[2K\r')
                appended_snp = set()
                for exon in transcript.getExons():
                    for var in vcf.fetch(exon.chrom, exon.beg - 2, exon.end):
                        if var.id in appended_snp:
                            cprint.printWarning("{} duplicated variation ID ...skipped".format(var.id))
                            outfiles.writeWarning([var.id, "duplicated variation ID", "skipped"])
                            continue
                        var_location = "..".join([str(var.start), str(var.stop)])
                        if var_location in appended_snp:
                            cprint.printWarning("{} duplicated variation position ...skipped".format(var.id))
                            outfiles.writeWarning([var.id, "duplicated variation position", "skipped"])
                            continue
                        
                        try:
                            snp = SNP(var, samples, optimization, tag_af)
                        except ValueError as err:
                            cprint.printWarning("{} - {} ...skipped".format(var.id, err.args[0]))
                            outfiles.writeWarning([var.id, err.args[0], "skipped"])
                            if DEBUG:
                                errlog.write("{} - {} ...skipped".format(var.id, err.args[0]))
                            continue
                        
                        if snp.filterAF(min_af):
                            if DEBUG and snp.id == 'rs529174822;rs563866452' or snp.id == 'rs72331392':
                                bp=1
                            
                            try:
                                exon.modify(snp)
                            except ValueError as err:
                                cprint.printWarning("{} - {} ...skipped".format(snp.id, err.args[0]))
                                outfiles.writeWarning([snp.id, err.args[0], "skipped"])
                            else:
                                if DEBUG:
                                    debugseq.modify(snp)
                                appended_snp.add(var.id)
                                appended_snp.add(var_location)
                                for err in snp.err:
                                    cprint.printWarning("{} - {} ...pay attention".format(snp.id, err))
                                    outfiles.writeWarning([snp.id, err, "pay attention"])
                    exon.cleanup()
                
                if DEBUG:
                    backbone_protein = transcript.getBackboneProtein()
                
                print("graph building...")
                transcript.makeGraph()
                trnVariations = transcript.translateVariations()
                
                if optimization:
                    print("graph optimization...")
                    if save_var == 'used':
                        trnVariations = transcript.joinSynonymPathes()
                    else:
                        transcript.joinSynonymPathes()
                
                if save_var: # save translated alleles in file
                    print("variation processing...")
                    for trn_var in trnVariations:
                        if optimization:
                            outfiles.writeVariation(transcript.id, trn_var, save_var)
                        else:
                            outfiles.writeVariation(transcript.id, trn_var)
                
                if do_prot:
                    print("protein processing...")
                    sample_proteins = transcript.getProteins(optimization)
                    for sample_name in sample_proteins:
                        if sample_name in miss_samples and re.search("Y", transcript.chrom):
                            continue # female sample
                        protdb.load(sample_name, sample_proteins[sample_name])
                
                if len(pept_len):
                    print("peptide processing...")
                    sample_peptides = transcript.getPeptides(pept_len, optimization)
                    for sample_name in sample_peptides:
                        if sample_name in miss_samples and re.search("Y", transcript.chrom):
                            continue # female sample
                        peptdb.load(sample_name, sample_peptides[sample_name]) # 'chrom', 'transcript_id', 'sample', 'allele1', 'allele2', 'beg', 'end', 'fshifts_before', 'variations(positions_in_matrix)', 'peptide', 'matrix'
                        if DEBUG:
                            debugseq.add2compare(sample_peptides[sample_name])
                
                if DEBUG:
                    debugseq.translate()
                    debugseq.peptides(pept_len)
                    for misspept in debugseq.compare():
                        errlog.write(misspept)
                
                print("...Ok\n")
                if len(transcript_set) and transcript_num < 1:
                    all_done = 1
                    break
            except ValueError as err:
                cprint.printWarning("transcript {} has a problem ({}) ...skipped".format(trn_id, err.args[0]))
                outfiles.writeWarning([trn_id, err.args[0], "skipped"])
        if all_done: break
    
    print("writing to output...")
    if len(pept_len):
        peptdb.flush()
        peptdb.doIndex()
        for sample in samples:
            if len(samples) > 1:
                output = peptdb.getUnique(sample)
            else:
                output = peptdb.getAll(sample)
            for rec in output:
                if not rec: break
                for row in rec:
                    outfiles.writePeptide(row)
        peptdb.close()
    
    if do_prot:
        protdb.flush()
        protdb.doIndex()
        for sample in samples:
            if len(samples) > 1:
                output = protdb.getUnique(sample)
            else:
                output = protdb.getAll(sample)
            for rec in output:
                if not rec: break
                for row in rec:
                    outfiles.writeProtein(row)
        protdb.close()
    
    notfound = set()
    for trnid, count in transcript_set.items():
        if not count:
            notfound.add(trnid)
    if len(notfound):
        for trnid in notfound:
            outfiles.writeWarning([trnid, "transcript not found", "not found"])
        cprint.printWarning("\nWARNING: {} transcript(s) not found (see warnings.csv file)\n".format(len(notfound)))
    
    stop_time = time.time()
    cprint.printOk("...the job is done (execution time: {})".format(time.strftime("%H:%M:%S", time.gmtime(stop_time - start_time))))
    if tmpdir_created and os.path.exists(tmp_dir):
        try:
            shutil.rmtree(tmp_dir)
        except:
            cprint.printFail("\nCan't remove the temporary directory!\n")
# end of main

if __name__ == '__main__':
    main()
