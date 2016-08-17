import yt
import numpy as np
import sys


filename = sys.argv[1]
field = sys.argv[2]

ds = yt.load(filename)
cube = ds.smoothed_covering_grid(0, [0,0,0], dims=ds.domain_dimensions)

data = cube[field]

data -= data.min()
data = data/data.max()
pointdata = data.flatten()
pointdata *= 255/pointdata.max()

#output to file
binfile = open('{:}_{:}_raw.raw'.format(filename[-4:], field), 'wb')
pointdata.astype(np.uint8).tofile(binfile)
