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
frame_count = 7

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
#rot_begin = 0. #frame to start animation on
#rot_end = 500. #frame to stop fade on
animation_end_frame = 810

#animation segment counters:
seg_1 = 140 #0-150 initial slow-mo time evolution to step 20
seg_2 = 420 #150-420 move camera to inside.  Then go back and forth w/ other blend file
seg_3 = 570 #420-570 zoom camera back out
seg_4 = 810 #570-810 rotate again a bit w/ time evolution

i = 0
raw_file_counter = 0

#print('the last argument is:'+str(sys.argv))

#--- make sure proper textures are displayed
bpy.data.materials['Material'].use_textures[0] = True
bpy.data.materials['Material'].use_textures[1] = False

for filepath in sorted(file_list):
    if (i < seg_1):
        #--- Render time evolution
        bpy.data.scenes["Scene"].frame_start = i
        bpy.data.scenes["Scene"].frame_end = i#+frame_count
        bpy.data.textures["hydrogen"].voxel_data.filepath = filepath

        #--- Start animating ---#
        #print("rendered frame (a): {:}, file: {:}".format(i, filepath))
        bpy.ops.render.render(animation=True)
        i = i+frame_count
        raw_file_counter += 1
    elif (i >= seg_1 and i < seg_3):
        bpy.data.scenes["Scene"].frame_start = i
        bpy.data.scenes["Scene"].frame_end = seg_3
        bpy.data.textures["hydrogen"].voxel_data.filepath = filepath

        #print("rendered frame (b): {:}, file: {:}".format(i, filepath))
        bpy.ops.render.render(animation=True)
        i = seg_3
#    elif (i >= seg_2 and i < seg_3):
#        bpy.data.scenes["Scene"].frame_start = i
#        bpy.data.scenes["Scene"].frame_end = i
#        bpy.data.textures["hydrogen"].voxel_data.filepath = filepath
#
#        print("rendered frame (b): {:}, file: {:}".format(i, filepath))
#        #bpy.ops.render.render(animation=True)
#        i += 1
    elif (i >= seg_3 and i <= animation_end_frame):
        bpy.data.scenes["Scene"].frame_start = i
        bpy.data.scenes["Scene"].frame_end = i
        #print('rendered frame (c): {:}, file: {:}'.format(i, filepath))
        i += 1
        bpy.ops.render.render(animation=True)
        #sys.exit()
