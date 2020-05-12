'''
Created on May 9, 2020

Given an old-style fwd-rev transition statistics pickle, will unpack and write
to XML.

@author: lvotapka
'''

import sys
import cPickle as pickle
import xml.etree.ElementTree as ET
from xml.dom import minidom

import base
import fwd_rev


print "Parse arguments"
if len(sys.argv) < 5:
    print "Usage:\npython convert_pickle.py SEEKR_PICKLE_FILE MILESTONE_INDEX" \
        " TRANSITION_FILE XML_FILE"
    
    exit()

picklename = sys.argv[1]
milestone_index = int(sys.argv[2])
transition_file = sys.argv[3]
xmlname = sys.argv[4]

me = base.openSeekrCalcPickle(picklename)
milestone = me.milestones[milestone_index]

transition_dict, avg_incubation_time, incubation_time_list = \
    fwd_rev.read_data_file_transitions(transition_file, me, milestone)
    
print 'transition_dict:', transition_dict
print 'avg_incubation_time:', avg_incubation_time
#print 'incubation_time_list:', incubation_time_list

root = ET.Element('fwd_transitions')
xmlIndex = ET.SubElement(root, 'milestone_index')
xmlIndex.text = str(milestone.index)
xmlIndex = ET.SubElement(root, 'milestone_siteid')
xmlIndex.text = str(milestone.siteid)
xmlIndex = ET.SubElement(root, 'milestone_radius')
xmlIndex.text = str(milestone.radius)
xmlTime = ET.SubElement(root, 'avg_incubation_time')
xmlTime.text = str(avg_incubation_time)
xmlTransitions = ET.SubElement(root, 'transitions')
for key in transition_dict:
    xmlTransition = ET.SubElement(xmlTransitions, 'transition')
    xmlSrc = ET.SubElement(xmlTransition, 'source')
    src = key.split('_')[0]
    xmlSrc.text = src
    xmlDest = ET.SubElement(xmlTransition, 'destination')
    dest = key.split('_')[1]
    xmlDest.text = dest
    dest_index = int(dest)
    dest_radius = me.milestones[dest_index].radius
    xmlDestRadius = ET.SubElement(xmlTransition, 'destRadius')
    xmlDestRadius.text = str(dest_radius)
    xmlCount = ET.SubElement(xmlTransition, 'count')
    count = transition_dict[key]
    xmlCount.text = str(count)

xmlstr = minidom.parseString(ET.tostring(root)).toprettyxml(
    indent="   ")
with open(xmlname, 'w') as xmlfile:
    xmlfile.write(xmlstr)
