'''
Created on Mar 31, 2020

Opens an old-style OpenSeekr SeekrCalculation pickle, reading the contents, and 
the writing them to an XML file. A test is also performed to ensure that all 
the components of the class in the pickle file were faithfully written to 

Must be run in Python 2!!

@author: lvotapka
'''

import old_base
import sys

print "Parse arguments"
if len(sys.argv) < 3:
    print "Usage:\npython convert_pickle.py PICKLE XML"
    
    exit()

picklename = sys.argv[1]
xmlname = sys.argv[2]

seekrPickle = old_base.openSeekrCalcPickle(picklename)
seekrPickle.save(xmlname)

print('Conversion complete.')