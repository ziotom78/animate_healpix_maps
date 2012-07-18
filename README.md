# animate_healpix_maps

This Python script takes a set of Healpix maps (in FITS format) and
produces an animated .gif file with Mollweide projections of the maps.

## Requirements

* Python 2.7 (maybe 2.6 will work too?)
* Matplotlib
* NumPy
* [Healpy](https://github.com/healpy/healpy)
* [Imagemagick](http://www.imagemagick.org) or
  [Graphicsmagick](http://www.graphicsmagick.org/)
  
## Usage

The tool must be used from the command line. Each FITS file must be
preceded by a title, usually enclosed within single/double quotes:

    $ animate_healpix_maps.py 'First' 1.fits 'Second' 2.fits
	
To see a list of available options, use `--help`.
