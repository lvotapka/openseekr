'''
Created on Mar 31, 2020

Opens an old-style OpenSeekr SeekrCalculation pickle, reading the contents, and 
the writing them to an XML file. A test is also performed to ensure that all 
the components of the class in the pickle file were faithfully written to 

Must be run in Python 2!!

@author: lvotapka
'''

import base
import sys

print "Parse arguments"
if len(sys.argv) < 3:
    print "Usage:\npython convert_pickle.py PICKLE_FILE XML_FILE"
    
    exit()

picklename = sys.argv[1]
xmlname = sys.argv[2]

seekrPickle = base.openSeekrCalcPickle(picklename)
print('xmlname', xmlname)
seekrPickle.save(xmlname)
testSeekr = base.openSeekrCalc(xmlname)

print('Conversion complete.')
