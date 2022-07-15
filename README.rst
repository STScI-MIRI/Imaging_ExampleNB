Example Notebooks MIRI Imaging Data Reduction
=============================================

Example notebooks and supporting python code for reducing MIRI Imaging
observations.  Extensively tested on MIRI Commissioning and Early
Release Observations (EROs).  Uses the
`official jwst pipeline <https://github.com/spacetelescope/jwst>`_.

`Presentation <https://speakerdeck.com/karllark/jwst-miri-imaging-processing-lessons-from-eros>`_
given to the MIRI Imaging Team 20 June 2022.

Extra processing
----------------

Extra steps were done to enhance the data.

1. After the Detector1 processing, the rate file was recreated from the rateints
file to correct a bug.  This bug affects the reported rate value in regions that
saturate in after the 3rd group.  This bug will be fixed in a future version of
the official jwst pipeline.

2. The astrometry (WCS, pointing) between mosaics tiles is updated based on
an analysis of the "tweakreg" results in the Image3 stage of the pipeline.
Experience shows that the relative astrometry between dithers using the same
guide star is quite good, but there can be significant (up to 2 arcsec)
updates needed between mosaic tiles where different guide stars are used.
The analysis results from the shortest wavelength filter are used for all
filters.  This is done as often there are not enough sources at longer
wavelengths to provide a good update to the astrometry.

3. In the case of sufficiently sparse fields, a sky/background image can be
created and subtracted from all the images.  This sky image is created by
sigma clipping all the images aligned in instrument coordinates.  This can
and does inject residual small scale variations in the background.  This
step was done for the SMACS 0723 and Stephan's Qunitet EROs and the LMC Flat
region.

LMCFlat Region (PID: 1040)
--------------------------

Portion of the LMC astrometric field used as part of
creating the Imaging flat fields.  A five mosaic tile observation with multiple
dithers in each tile.  All 9 Imaging filters observed.

Notebooks:

* LMCFlat_F560W: Reduction of F560W data.

* LMCFlat_F560W_skysub: Creation and subtraction of a sky/background image.

* LMCFlat_F770W: Reduction of F770W data.  Note that the updates to the
  astronomery used the results from the F560W tweakreg analysis.

* LMCFlat_F770W_skysub: Creation and subtraction of a sky/background image.

Contributors
------------
Karl Gordon

License
-------

This code is licensed under a 3-clause BSD style license (see the
``LICENSE`` file).
