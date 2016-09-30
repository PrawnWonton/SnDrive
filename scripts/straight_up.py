#Blender command line scripting attempt
#
#--- Changelog ---
#-----------------
# for use with
#--------------------------------------

import bpy
#import numpy as np
import sys
import os
import glob
import optparse

#number of frames per timestep
frame_count = 0

###
o = optparse.OptionParser()
o.set_description('Suck a bag of dicks')

o.add_option('--glob_pat',type=str,
    help='location of .raw files for glob')

opts,args = o.parse_args(sys.argv[5:])
###

file_list = glob.glob(opts.glob_pat)
animate_length = len(file_list)
skip_length = 0

#initialize counters
rot_begin = 0. #frame to start animation on
rot_end = 300. #frame to stop fade on
animation_end_frame = 261
i = 0
raw_file_counter = 20

#print('the last argument is:'+str(sys.argv))

#--- make sure proper textures are displayed
#bpy.data.materials['Material'].use_textures[2] = True
#bpy.data.materials['Material'].use_textures[0] = False

for filepath in sorted(file_list):
    if (i < raw_file_counter):
        print("skipped file (a): {:}, file: {:}".format(i, filepath))
        i += 1
    elif (i >= raw_file_counter and i <= animation_end_frame):
        #--- Render time evolution
        bpy.data.scenes["Scene"].frame_start = i
        bpy.data.scenes["Scene"].frame_end = i#+frame_count
        bpy.data.textures["hydrogen"].voxel_data.filepath = filepath

        #--- Start animating ---#
        #print("rendered frame (b): {:}, file: {:}".format(i, filepath))
        bpy.ops.render.render(animation=True)
        i = i+1
        #raw_file_counter += 1
