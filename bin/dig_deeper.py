'''
This script extracts a successful forward stage trajectory that ends on a lower
milestone. It extracts the last frame of that trajectory and writes the
structure to the next milestone down's directory.

Created on June 28, 2018

@author: lvotapka
'''

import mdtraj
import seekr
import sys, os
from shutil import copyfile

def read_data_file_transition_down(data_file_name, destination='1', last_frame=True):
  '''Read transition data file, return the first instance of a transition to a
  lower milestone. Return the file index of the forward trajectory to this lower
  milestone.
  Input:
   - data_file_name: string of the filename of a milestone's forward transition
   file
   - destination: the string to search for. In the downward milestone case: '1'.
      Upward would be '3'.
  Output:
   - downward_index: integer of the forward trajectory that ends on the lower
   milestone
  '''
  downward_index = None
  data_file = open(data_file_name, 'r')
  for i, line in enumerate(data_file.readlines()):
    line = line.split()
    if line[0] == destination:
      downward_index = i
      if not last_frame: # then read the first frame
        break 
  data_file.close()
  assert downward_index != None, "FAILURE: no downward forward trajectories detected. Please run additional umbrella sampling and rev/fwd trajectories."
  return downward_index

print "Parse arguments"
if len(sys.argv) != 3:
  print "Usage:\npython dig_deeper.py milestone pickle"
  exit()

which = int(sys.argv[1])
picklename = sys.argv[2]

print "Loading SEEKR calculation."
me = seekr.openSeekrCalc(picklename)

milestone = me.milestones[which]
lower_milestone = me.milestones[which-1]

# define all directories and files
fwd_rev_dir = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev')
lower_temp_equil_dir = os.path.join(me.project.rootdir, lower_milestone.directory, 'md', 'temp_equil')
lower_temp_equil_filename = os.path.join(lower_temp_equil_dir, 'equilibrated.pdb')
lower_milestone_holo = os.path.join(me.project.rootdir, lower_milestone.directory, 'md', 'holo_wet.pdb')
lower_milestone_building = os.path.join(me.project.rootdir, lower_milestone.directory, 'md', 'building')
data_file_name = os.path.join(fwd_rev_dir, 'transition_fwd.dat') # find the index of a successful downward trajectory
prmtop = os.path.join(me.project.rootdir, milestone.directory, 'md', 'building', 'holo.parm7')
new_prmtop = os.path.join(lower_milestone_building, 'holo.parm7')
new_inpcrd = os.path.join(lower_milestone_building, 'holo.rst7')

# figure out which forward to pull out from
print "Attempting to extract the last frame of a successful downward trajectory."
downward_index = read_data_file_transition_down(data_file_name)
downward_fwd_dcd = os.path.join(fwd_rev_dir, 'forward%i_0.dcd' % downward_index)
print "Extracting frame from file:", downward_fwd_dcd

print "Writing new structures and files needed to run umbrella simulation on the lower milestone (milestone %d)" % lower_milestone.index
#last_fwd_frame = mdtraj.load(downward_fwd_dcd, top=prmtop)[-1]
last_fwd_frame = seekr.load_last_mdtraj_frame(downward_fwd_dcd, prmtop)
last_fwd_frame.save_pdb(lower_temp_equil_filename)
last_fwd_frame.save_pdb(lower_milestone_holo)
last_fwd_frame.save_amberrst7(new_inpcrd)

# copy the prmtop to the lower building directory
copyfile(prmtop, new_prmtop)
lower_milestone.openmm.prmtop_filename = new_prmtop
lower_milestone.openmm.inpcrd_filename = new_inpcrd

me.save()
