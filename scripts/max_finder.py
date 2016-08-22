import yt
import numpy as np
import sys, os, glob, optparse

###
o = optparse.OptionParser()
o.set_description('Suck a bag of dicks')

o.add_option('--glob_pat',type=str,
    help='location of .raw files for glob')
o.add_option('--field', type=str,
    help='name of field')

opts,args = o.parse_args(sys.argv[1:])
###
file_list = glob.glob(opts.glob_pat)
field = opts.field

max_list = []
min_list = []

for filepath in sorted(file_list):
    ds = yt.load(filepath)
    cube = ds.smoothed_covering_grid(0, [0,0,0], dims=ds.domain_dimensions)
    data = cube[field]
    max_list.append(data.max())
    min_list.append(data.min())

super_max = np.max(max_list)
super_min = np.min(min_list)

for dataset in sorted(file_list):
    ds = yt.load(dataset)
    cube = ds.smoothed_covering_grid(0, [0,0,0], dims=ds.domain_dimensions)
    data = cube[field]
    data = data.flatten()
    data -= super_min
    data = data/super_max
    pointdata = data.flatten()
    pointdata *= 255

    #output to file
    binfile = open('{:}_{:}.raw'.format(dataset[-4:], field), 'wb')
    pointdata.astype(np.uint8).tofile(binfile)
