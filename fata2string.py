#!/usr/bin/env python3

import sys
import re
import argparse

def main():
    input_parser = argparse.ArgumentParser(description='fasta2string - convert fasta files to one string direct access files.')
    input_parser.add_argument('-fasta', metavar='file.fasta', default=None, help='FASTA input file', required=True)
    
    args = input_parser.parse_args()
    fasta_file = args.fasta
    
    path = ''
    re_path = re.match(r"(.*\/)", fasta_file)
    if re_path:
        path = re_path.group(1)
    
    outfile = None
    with open(fasta_file, 'r') as infile:
        for line in infile:
            re_header = re.match(r">(\S+)", line)
            if re_header:
                header = re_header.group(1)
                if outfile:
                    outfile.close()
                outfile = open(path + header + ".seq", "w")
                outfile.write(">" + header + "\n")
                print("file {}.seq created".format(header))
            else:
                line = re.sub(r"[^A-Za-z]", "", line)
                if outfile:
                    outfile.write(line)
    if outfile:
        outfile.close()
    print("...done")

if __name__ == '__main__':
    main()
