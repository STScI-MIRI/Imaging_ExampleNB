import os.path

from jwst.pipeline import calwebb_detector1
from jwst.pipeline import calwebb_image2
from jwst.pipeline import calwebb_image3


def miri_detector1(
    miri_uncal_files,
    output_dir,
    jump_sigma=5.0,
    maskfile=None,
    resetfile=None,
    darkfile=None,
    linfile=None,
    save_jump_info=False,
    firstframe=True,
    darksub=True,
    onlyifnorate=False,
    cpufraction="half",
    after_jump=False,
):
    """
    Run CALWEBB_DETECTOR1 pipeline on uncal files.
    """
    for miri_uncal_file in miri_uncal_files:
        print(miri_uncal_file)

        miri_rate_file = miri_uncal_file.replace("stage0", "stage1").replace("uncal", "rate")
        if onlyifnorate and os.path.isfile(miri_rate_file):
            print(f"{miri_rate_file} already exists, skipping detector1 for this file")
            pass
        else:

            # Using the run() method:
            # Instantiate the pipeline
            miri1 = calwebb_detector1.Detector1Pipeline()

            # Save the final output of the pipeline
            miri1.output_dir = output_dir
            miri1.save_results = True

            # override the mask file
            # miri1.dq_init.override_mask = "./RefFiles/MIRI_MIRIMAGE_MASK_09.00.04.fits"

            if not firstframe:
                miri1.firstframe.skip = True

            # override the rscd file
            # miri1.rscd.override_rscd = "./RefFiles/MIRI_MIRIMAGE_RSCD_09.03.00.fits"

            # Use the flight saturation file
            # miri1.saturation.override_saturation = "./RefFiles/MIRIMAGE_SATURATION_09.00.00.fits"

            # Use the flight reset anomaly file
            if resetfile is not None:
                miri1.reset.override_reset = f"./RefFiles/{resetfile}"

            # Use the flight reset anomaly file
            if linfile is not None:
                miri1.linearity.override_linearity = f"./RefFiles/{linfile}"

            # Use the flight dark file
            if darkfile is not None:
                miri1.dark_current.override_dark = f"./RefFiles/{darkfile}"
            if not darksub:
                miri1.dark_current.skip = True

            # Save the output from the jump detecion step,
            # and let's use a more stringent limit for cosmic ray
            # detection
            # miri1.jump.save_results = True
            miri1.jump.rejection_threshold = jump_sigma
            miri1.jump.maximum_cores = cpufraction

            if after_jump:
                miri1.jump.after_jump_flag_dn1 = 10.0
                miri1.jump.after_jump_flag_time1 = 15 * 2.8
                miri1.jump.after_jump_flag_dn2 = 1000.0
                miri1.jump.after_jump_flag_time2 = 1e4 * 2.8
                # miri1.jump.after_jump_flag_time2 = 0

            # Save the linearized ramp
            # miri1.linearity.save_results = True

            # Save the results of the ramp fitting
            if save_jump_info:
                miri1.save_calibrated_ramp = True
                miri1.ramp_fit.save_results = True
                miri1.ramp_fit.save_opt = True
            miri1.ramp_fit.maximum_cores = cpufraction

            # Call the run() method
            miri1.run(miri_uncal_file)


def miri_image2(miri_rate_files, output_dir,
                flatfile=None,
                photomfile=None):
    """
    Run CALWEBB_IMAGE2 pipeline on rate files.
    """
    # Using the run() method
    # Create an instance of the pipeline class
    miri_im2 = calwebb_image2.Image2Pipeline()

    # Save the pipeline output, and specify the output directory
    miri_im2.output_dir = output_dir
    miri_im2.save_results = True

    # use a flat field that has high/low responsive pixels removed
    if flatfile is not None:
        miri_im2.flat_field.override_flat = flatfile

    # Use the flight reset anomaly file
    if photomfile is not None:
        miri_im2.photom.override_photom = f"./RefFiles/{photomfile}"

    # Set any step-related paramters
    miri_im2.resample.pixfrac = 1.0

    # Call the run() method and provide the MIRI association file as input
    miri_im2.run(miri_rate_files)


def miri_image3(miri_asn_file, output_dir, minobj=5, snr=5, fwhm=None,
                crval=None, crpix=None, rotation=None, output_shape=None,
                align_to_gaia=False, tweakreg=True,
                matchbkg=True, pixel_scale=0.11):
    """
    Run CALWEBB_IMAGE3 pipeline on asn file.
    """
    # Create an instance of the pipeline class
    im3 = calwebb_image3.Image3Pipeline()

    im3.output_dir = output_dir

    # parameters based on Mattia's work

    sigma = 3.  # clipping limit, in sigma units, used when performing fit, default=3
    fit_geom = "shift"  # ftype of affine transformation to be considered when fitting catalogs, default='general'
    use2dhist = True  # boolean indicating whether to use 2D histogram to find initial offset, default=True

    if not tweakreg:
        im3.tweakreg.skip = True
    if fwhm is not None:
        im3.tweakreg.kernel_fwhm = fwhm
    im3.tweakreg.save_results = True
    im3.tweakreg.snr_threshold = snr
    im3.tweakreg.minobj = minobj
    im3.tweakreg.sigma = sigma
    im3.tweakreg.fitgeometry = fit_geom
    im3.tweakreg.use2dhist = use2dhist
    # # im3.source_catalog.kernel_fwhm = fwhm
    im3.source_catalog.snr_threshold = snr
    # im3.skymatch.save_results = True
    im3.skymatch.lower = 0.0
    # im3.outlier_detection.save_results = True
    # im3.resample.save_results = True
    # im3.source_catalog.save_results = True
    im3.save_results = True

    # currently does not work always or anytime?
    im3.tweakreg.align_to_gaia = align_to_gaia
    im3.tweakreg.searchrad = 5.0
    # im3.tweakreg.separation = 2.0
    # im3.tweakreg.tolerance = 2.0

    if crval is not None:
        im3.resample.crval = crval
    if crpix is not None:
        im3.resample.crpix = crpix
    im3.resample.kernel = "square"
    im3.resample.pixel_scale = pixel_scale
    if rotation is not None:
        im3.resample.rotation = rotation
    if output_shape is not None:
        im3.resample.output_shape = output_shape

    if not matchbkg:
        im3.skymatch.skip = True

    # Call the run() method
    im3.run(miri_asn_file)
