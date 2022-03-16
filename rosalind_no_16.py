# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 16:02:29 2021

@author: ashdo
"""
#no 16
from Bio import SeqIO
from Bio import Seq
from Bio import motifs
import urllib.request
import re
import regex

entries = []
for line in open('no16.txt'):
    given_names = line.rstrip('\n')
    entries.append(given_names)
    #print(entries,'\n')
    with open('16out.txt','w') as out:
        for k in entries:
            url = 'http://uniprot.org/uniprot/'+k+'.fasta'
            handle = urllib.request.urlopen(url)
            entry_info = handle.read().decode('utf-8')
            out.write(entry_info)
            #print(entry_info)
    for record in SeqIO.parse('16out.txt','fasta'):
        motif = 'N[^P][ST][^P]'
        names = record.name
        #print(names)
        seqs = str(record.seq)
        find_n = [m.start() + 1 for m in regex.finditer(motif,seqs, flags=0,overlapped=True)]
        motif_finds = str(find_n)
        outs = motif_finds.strip("[").strip("]").strip(",")
        out_motif = outs.replace(',','')
    #print(names, motif_finds)
    print(given_names,'\n',out_motif)