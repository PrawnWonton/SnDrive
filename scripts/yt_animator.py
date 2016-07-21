
#-------------------#
#--- Description ---#
#-------------------#
# iterate thru all plot files in a directory, slice (density, temp, etc.)
# centered at the first supernova.  Output plots as pictures with name scheme
# following frame_x.png, etc.
#

import yt
import sys, os
import optparse
import fnmatch
import path

#---

#o = optparse.OptionParser()
#o.set_description('Give as argument file name of first frame')
#o.add_option('--filename', type=string, help='name of file of first frame')
#opts, args = o.parse_args(sys.argv[1:])

#p = path(location)

inF = open('flash.log', 'r')
for line in inF:
    if 'x_sn' in line:
        for word in line.split():
            xvalue = word
        for word in inF.next().split():
            yvalue = word
        for word in inF.next().split():
            zvalue = word
        break
center = [float(xvalue), float(yvalue), float(zvalue)]

for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*0001'):
        first_plt = yt.load(file)
        sp = first_plt.sphere(center, (10., 'pc'))
        dens_zmin = float(sp['density'].min()/1000)
        dens_zmax = float(sp['density'].max()*10)

for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*plt_cnt_0*'):


        dataset = yt.load(file)
        width_pc = int(dataset.domain_width.in_units('pc')[0])
        year = round(dataset.current_time.in_units('Myr'), 5)
        j=int(file[-3:])

        slc_dens = yt.SlicePlot(dataset, 'z', 'density', width = (width_pc, 'pc'), center = center)
        slc_dens.set_log('density', True)
        slc_dens.set_zlim('density', dens_zmin, dens_zmax)
        slc_dens.save('Frame_{:}'.format(j))

        #print(j)
