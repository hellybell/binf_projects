#MS/MS VIEWER

import sys
import math

sequence = sys.argv[1]

#amino acid molecular weight dict for theoretical mz values

mw = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
          'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
          'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
          'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }
#print(mw.items())

#Calculate bion values of sequence 
b_dict = {}
bions = 0
for i,aa in enumerate(sequence[:-1], start=1):
    key = 'b'+ str(i)
    bions += mw[aa]
    bmz = bions+1
    value = bmz
    b_dict[key] = value

#print('b-values',b_dict)

#Calculate yion values of sequence
revsequence = sequence[::-1]
y_dict = {} 
yions = 0
for i, aa in enumerate(revsequence[:-1], start=1):
    key = 'y' + str(i)
    yions += mw[aa]
    ymz = yions+19
    value = ymz
    y_dict[key] = value
#print('y', y_dict)

    
#Extract peaks from mzXMLfile    
from matplotlib.offsetbox import AnchoredText 
import sys
import gzip
import xml.etree.ElementTree as ET
from base64 import b64decode
from array import array
import matplotlib.pyplot as plt
import numpy as np

#use input scan number 


file = sys.argv[2]
number = int(sys.argv[3])
xmldoc = gzip.open(file)

ns = '{http://sashimi.sourceforge.net/schema/}'


#extract the correct element of file using ET.iterparse
for event,ele in ET.iterparse(xmldoc):
    if ns == '{http://sashimi.sourceforge.net/schema/}':
        p = ele.tag.find('}')
        if p >= 0:
            ns = ele.tag[:(p+1)]
        if event == 'end' and ele.tag == ns+'scan':
            if int(ele.attrib['num']) == number:
                peaksele = ele.find(ns+'peaks')
##              print(peaksele.text)
            ele.clear()

#Use *special* code to convert peaksele.text into something coherent
            
peaks = array('f',b64decode(peaksele.text))
#print(peaks, '\n')
if sys.byteorder != 'big':
    peaks.byteswap()

#From peaks, calculate mzs and intensity values of the spectra file    
mzs = peaks[::2]
#print(mzs)
ints = peaks[1::2]
max_int = max(ints)
thr = (0.05 * max_int)

#find matches between sequence values and scan values, use data structure to compile matches
matches = []
for blabel,bmz in b_dict.items():
    for i in range(len(ints)):
        mz = mzs[i]
        bintensity = ints[i]
        if abs(bmz-mz) < 0.1:
            if bintensity >= thr:
                matches.append(blabel)
                matches.append(bmz)
                matches.append(bintensity)

for ylabel, ymz in y_dict.items():
    for i in range(len(ints)):
        mz = mzs[i]
        yintensity = ints[i]
        if abs(ymz-mz) < 0.02:
           if yintensity >= thr:
                matches.append(ylabel)
                matches.append(ymz)
                matches.append(yintensity)
                
#print(matches)

#First plot the spectra from the file (aka, scan values)                
x = mzs
y = ints
fig, ax = plt.subplots()
at = AnchoredText(sequence,
                  prop=dict(size=15), frameon=True, loc='upper right')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)
ax.set_title('Peptide fragmentation spectra of Scan No. with matches in Sequence')
ax.set_xlabel('m/z values')                
ax.set_ylabel('intensity')

markerline, stemlines, baseline = plt.stem(x,y,linefmt='grey', markerfmt='D',bottom=1.1,use_line_collection=True)
markerline.set_markerfacecolor('none')

print('The peptide sequence provided has the following matches to the scan No. spectra within threshold criteria')


#Then, annotate the matches from data structure onto spectra, and print out matches
for i in range(0,len(matches),3):
    match_label = matches[i]
    match_mz = matches[i+1]
    match_int = matches[i+2]
    match_row = (match_label,match_mz,match_int)
    
    print(match_row)
    
    plt.text((match_mz-4),(match_int+17),match_label, fontsize=12,)

plt.show()




