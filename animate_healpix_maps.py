#!/usr/bin/env python

'''Assemble a number of Healpix FITS maps into a GIF animated file.
If the `--mask' argument is present, a mask will be applied to every
map. Use --column to specify which column of the FITS files to read (default is
0, which usually is the I Stokes parameter).'''

import matplotlib
matplotlib.use('Agg') # Non-GUI backend

import matplotlib.pyplot as plt
import numpy as np
import healpy
import os
import sys
import tempfile
import logging as log
from optparse import OptionParser
from collections import namedtuple

#####################################################################


def parse_command_line():
    "Interpret the command line parameters using an OptionParser object"
    parser = OptionParser(usage="%prog OUTPUT_FILE_NAME MAP1 COEFF1 [MAP2 COEFF2] ...")

    parser.add_option('--mask', '-m', metavar='FILE',
                      default=None, dest='mask',
                      help='Path to the FITS file containing the mask to be applied')
    parser.add_option('--no-common-scale', action='store_true',
                      dest='no_common_scale', default=False,
                      help='Do not use a common scale for the colorbar')
    parser.add_option('--component', '-c', default='I', metavar='STRING',
                      dest='stokes_component',
                      help='Stokes'' component to plot (default: "%default")')
    parser.add_option('--output', '-o', metavar="FILE",
                      dest='output_file_name', default='output.gif',
                      help='Name of the output file that will contain the map'
                      ' (default: "%default")')
    parser.add_option('--delay', '-d', metavar='TIME',
                      dest='frame_delay', type='int', default=75,
                      help='Time to wait between frames, in hundreds of second '
                      '(default: %default)')
    parser.add_option('--power-spectra', action='store_true',
                      dest='plot_spectra', default=False,
                      help='Include a plot of the power spectra in the animation')
    parser.add_option('--pixel-distributions', action='store_true',
                      dest='plot_distributions', default=False,
                      help='Include a histogram of the pixels'' values')

    return parser.parse_args()

#####################################################################


def validate_args(options, args):
    # We want to convert a string of Stokes components like "I,q"
    # into something recognized by healpy.read_map, i.e. (0, 1).
    # The use of "frozenset" removes duplicates.
    component_map = {'I': 0, 'Q': 1, 'U': 2}
    try:                                                                                                                                                                                                     
        stokes_component = component_map[options.stokes_component.upper()]
    except KeyError:
        log.fatal('Unknown Stokes component %s in string "%s" '
                  '(available choices are %s)',
                  sys.exc_value,
                  options.stokes_components,
                  ', '.join(['"%s"' % x for x in component_map.keys()]))
        sys.exit(1)

    # Now overwrite options.stokes_components: we do not need the
    # user-provided string any longer
    options.stokes_component = stokes_component
    
#####################################################################


def create_namedtuple_for_map(options):    
    tuple_fields = ['title', 'pixels']
    if options.plot_distributions:
        tuple_fields.append('histogram')
    
    return namedtuple('Map', tuple_fields)
    
#####################################################################


def hist_x_axis_points(bins):
    '''If BINS is a list with N+1 elements containing the bin edges
    of a histogram computed with hp.histogram, return a N-element
    list with the bin mid-points.'''
    return (bins[:-1] + bins[1:]) * 0.5
    
    
#####################################################################

log.basicConfig(level = log.DEBUG,
                format = '[%(asctime)s %(levelname)s] %(message)s')

OPTIONS, ARGS = parse_command_line()
if len(ARGS) < 2:
    sys.stderr.write('Error: at least a map and a title are expected '
                     'on the command line\n')
    sys.exit(1)

validate_args(OPTIONS, ARGS)

map_title_name_pairs = zip(*[iter(ARGS)]*2)
Map = create_namedtuple_for_map(OPTIONS)

maps_read = []
for cur_title, cur_map_file in map_title_name_pairs:
    try:
        log.info('reading map `%s\'...' % cur_map_file)
        cur_map = healpy.read_map(cur_map_file, field=OPTIONS.stokes_component)
        cur_map[cur_map < -1.6e+30] = np.NaN

        if type(OPTIONS.mask) is not None:
            cur_map[OPTIONS.mask == 0] = np.NaN

        cur_entry = dict(title=cur_title,
                         pixels=cur_map)
        if OPTIONS.plot_distributions:
            cur_entry['histogram'] = np.histogram(cur_map[~np.isnan(cur_map)],
                                                  bins=50)
        maps_read.append(Map(**cur_entry))
        log.info('...ok')
    except:
        log.info('...error, skipping this map')
        raise

min_temperature = None
max_temperature = None
if not OPTIONS.no_common_scale:
    min_temperature = np.nanmin([np.nanmin(x.pixels) for x in maps_read])
    max_temperature = np.nanmax([np.nanmax(x.pixels) for x in maps_read])

IMAGE_SCALE = 2.5
png_file_names = []

num_of_plots = 1
if OPTIONS.plot_spectra:
    num_of_plots += 1
if OPTIONS.plot_distributions:
    num_of_plots += 1

for cur_entry in maps_read:
    try:
        log.info('plotting the map using "%s" as title',
                 cur_entry.title)
        fig = plt.figure(1, (3 * IMAGE_SCALE,
                             2 * IMAGE_SCALE * num_of_plots))
        healpy.mollview(cur_entry.pixels,
                        title=cur_entry.title,
                        fig=1,
                        sub=(num_of_plots * 100 + 10 + 1),
                        min=min_temperature,
                        max=max_temperature)
        
        next_plot = 2
        if OPTIONS.plot_distributions:
            plt.subplot(num_of_plots, 1, next_plot)
            num_of_bins = 50
            for i in maps_read:
                n, bins = i.histogram
                log.debug("#bins = %d, #x = %d, #y = %d",
                          len(bins),
                          len(hist_x_axis_points(bins)),
                          len(n))
                log.debug("n = %s", str(n))
                plt.plot(hist_x_axis_points(bins), n,
                         label=i.title)
            # Fill the histogram for the current map
            plt.fill_between(hist_x_axis_points(cur_entry.histogram[1]),
                             cur_entry.histogram[0],
                             alpha=0.5)
            plt.xlabel('Value of the pixel')
            plt.ylabel('Fractions of pixels in the map')
            plt.title('Distribution of the pixel values')
            plt.legend()
            plt.grid()

        file_name = tempfile.mktemp() + '.png'
        log.info('saving the map in temporary file `%s\'' % file_name)
        fig.savefig(file_name)
        plt.clf()
        png_file_names.append(file_name)
    except:
        sys.stderr.write('Unable to produce a map for `%s'', skipping\n'
                         % cur_map_file)
        raise

convert_cmd_line = ' '.join(['convert',
                             '-delay {0}'.format(OPTIONS.frame_delay),
                             '-loop 0',
                             ' '.join(png_file_names),
                             OPTIONS.output_file_name])
log.info('running command `%s\'' % convert_cmd_line)
os.system(convert_cmd_line)

for cur_file in png_file_names:
    try:
        log.info('trying to remove file `%s\'...' % cur_file)
        os.unlink(cur_file)
        log.info('...done')
    except:
        log.info('...file not found')
