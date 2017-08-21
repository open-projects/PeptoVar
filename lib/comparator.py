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

import os, errno
import sqlite3
import re
import tempfile

class UniPep:
    def __init__(self, sample_name1, sample_name2, dir_name = ":memory:"):
        self.names = {sample_name1: {}, sample_name2: {}}
        self.names[sample_name1]['insert'] = "INSERT INTO pept1 VALUES (?,?,?,?,?,?,?,?,?,?,?)"
        self.names[sample_name2]['insert'] = "INSERT INTO pept2 VALUES (?,?,?,?,?,?,?,?,?,?,?)"
        self.names[sample_name1]['unique'] = "SELECT t1.chrom, t1.transcript_id, t1.sample, SUM(t1.allele1), SUM(t1.allele2), t1.beg, t1.end, t1.fshifts, t1.snp, t1.peptide, t1.matrix FROM pept1 t1 \
                                                LEFT JOIN pept2 t2 USING(peptide) WHERE t2.peptide IS NULL \
                                                GROUP BY t1.chrom, t1.transcript_id, t1.sample, t1.beg, t1.end, t1.fshifts, t1.snp, t1.peptide, t1.matrix ORDER BY t1.chrom, t1.transcript_id, t1.beg, t1.end"
        self.names[sample_name2]['unique'] = "SELECT t1.chrom, t1.transcript_id, t1.sample, SUM(t1.allele1), SUM(t1.allele2), t1.beg, t1.end, t1.fshifts, t1.snp, t1.peptide, t1.matrix FROM pept2 t1 \
                                                LEFT JOIN pept1 t2 USING(peptide) WHERE t2.peptide IS NULL \
                                                GROUP BY t1.chrom, t1.transcript_id, t1.sample, t1.beg, t1.end, t1.fshifts, t1.snp, t1.peptide, t1.matrix ORDER BY t1.chrom, t1.transcript_id, t1.beg, t1.end"
        self.names[sample_name1]['allpep'] = "SELECT chrom, transcript_id, sample, SUM(allele1), SUM(allele2), beg, end, fshifts, snp, peptide, matrix FROM pept1 \
                                                GROUP BY chrom, transcript_id, sample, beg, end, fshifts, snp, peptide, matrix ORDER BY chrom, transcript_id, beg, end, allele1, allele2"
        self.names[sample_name2]['allpep'] = "SELECT chrom, transcript_id, sample, SUM(allele1), SUM(allele2), beg, end, fshifts, snp, peptide, matrix FROM pept2 \
                                                GROUP BY chrom, transcript_id, sample, beg, end, fshifts, snp, peptide, matrix ORDER BY chrom, transcript_id, beg, end, allele1, allele2"

        self._file, self._filename = tempfile.mkstemp(suffix='.db', dir=dir_name)
        os.close(self._file)
        
        self._con = sqlite3.connect(self._filename)
        self._cur = self._con.cursor()
        self._cur.execute("CREATE TABLE pept1 (chrom TEXT NOT NULL, transcript_id TEXT NOT NULL, sample TEXT NOT NULL, allele1 integer NOT NULL, allele2 integer NOT NULL, beg INTEGER NOT NULL, end INTEGER NOT NULL, fshifts TEXT NOT NULL, snp TEXT NOT NULL, peptide TEXT, matrix TEXT NOT NULL)")
        #self._cur.execute("CREATE INDEX peptide1 ON pept1 (peptide)")
        self._cur.execute("CREATE TABLE pept2 (chrom TEXT NOT NULL, transcript_id TEXT NOT NULL, sample TEXT NOT NULL, allele1 integer NOT NULL, allele2 integer NOT NULL, beg INTEGER NOT NULL, end INTEGER NOT NULL, fshifts TEXT NOT NULL, snp TEXT NOT NULL, peptide TEXT, matrix TEXT NOT NULL)")
        #self._cur.execute("CREATE INDEX peptide2 ON pept2 (peptide)")
    
    def doIndex(self):
        self._cur.execute("CREATE INDEX IF NOT EXISTS peptide1 ON pept1 (peptide)")
        self._cur.execute("CREATE INDEX IF NOT EXISTS peptide2 ON pept2 (peptide)")
        
    
    def load(self, sample_name, data):
        if len(data) and sample_name in self.names:
            self._cur.executemany(self.names[sample_name]['insert'], data)
    
    def flush(self):
        self._con.commit()
    
    def getUnique(self, sample_name):
        self._cur.execute(self.names[sample_name]['unique'])
        return self._cur.fetchall()
    
    def getAll(self, sample_name):
        self._cur.execute(self.names[sample_name]['allpep'])
        return self._cur.fetchall()
    
    def close(self):
        self._con.close()
        try:
            os.remove(self._filename)
        except OSError as e:
            if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
                print("WARNING: can't remove temporary database file")
# end of UniPep
