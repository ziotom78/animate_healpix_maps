# animate_healpix_maps

This Python script takes a set of Healpix maps (in FITS format) and
produces an animated .gif file with Mollweide projections of the maps.

## Requirements

* Python 3.6 (maybe older versions 3.x will work too?)
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

## Example

The following animation has been produced using the maps delivered by
the COBE/DMR team for the two 31 GHz channels (download them from
[here](http://lambda.gsfc.nasa.gov/product/cobe/dmr_4year_skymaps_get.cfm)).
The command used to produce the .gif file was

    $ animate_healpix_maps.py \
	    -s 0.75 \
		--pixel-distributions \
	    -o dmr_31.gif \
		'DMR(31a)' dmr_31a_imap_4yr.fits \
		'DMR(31b)' dmr_31b_imap_4yr.fits

![Animation][example]

[example]: https://github.com/ziotom78/animate_healpix_maps/raw/master/examples/dmr_31.gif
