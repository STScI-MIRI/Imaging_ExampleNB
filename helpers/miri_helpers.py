import os.path

from jwst.pipeline import calwebb_detector1
from jwst.pipeline import calwebb_image2
from jwst.pipeline import calwebb_image3


def miri_detector1(
    miri_uncal_files,
    output_dir,
    jump_sigma=5.0,
    maskfile=None,
    reset=True,
    resetfile=None,
    darkfile=None,
    linfile=None,
    save_jump_info=False,
    firstframe=True,
    darksub=True,
    onlyifnorate=False,
    cpufraction="half",
    after_jump=False,
    logfile=None,
):
    """
    Run CALWEBB_DETECTOR1 pipeline on uncal files.
    """
    # setup the directory of step paramters
    det1_dict = {}
    det1_dict["ipc"] = {"skip": True}
    if not firstframe:
        det1_dict["firstframe"] = {"skip": True}

    # Use the flight reset anomaly file
    if not reset:
        det1_dict["reset"] = {"skip": True}
    if resetfile is not None:
        det1_dict["reset"] = {"override_reset": f"./RefFiles/{resetfile}"}

    # Use the flight reset anomaly file
    if linfile is not None:
        det1_dict["linearity"] = {"override_linearity": f"./RefFiles/{linfile}"}

    # Use the flight dark file
    if not darksub:
        det1_dict["dark_current"] = {"skip": True}
    if darkfile is not None:
        det1_dict["dark_current"] = {"override_dark": f"./RefFiles/{darkfile}"}

    # Save the output from the jump detecion step,
    # miri1.jump.save_results = True
    det1_dict["jump"] = {"rejection_threshold": jump_sigma,
                         "maximum_cores": cpufraction}
    if after_jump:
        det1_dict["jump"]["after_jump_flag_dn1"] = 10.0
        det1_dict["jump"]["after_jump_flag_time1"] = 15 * 2.8
        det1_dict["jump"]["after_jump_flag_dn2"] = 1000.0
        det1_dict["jump"]["after_jump_flag_time2"] = 1e4 * 2.8

    # Save the linearized ramp
    # det1_dict["linearity"] = {"save_results": True}

    # Save the results of the ramp fitting
    det1_dict["ramp_fit"] = {"maximum_cores": cpufraction}
    if save_jump_info:
        det1_dict["ramp_fit"]["save_results"] = True
        det1_dict["ramp_fit"]["save_opt"] = True
        save_calibrated_ramp = True
    else:
        save_calibrated_ramp = False

    for miri_uncal_file in miri_uncal_files:
        print(miri_uncal_file)

        miri_rate_file = miri_uncal_file.replace("stage0", "stage1").replace("uncal", "rate")
        if onlyifnorate and os.path.isfile(miri_rate_file):
            print(f"{miri_rate_file} already exists, skipping detector1 for this file")
            pass
        else:
            calwebb_detector1.Detector1Pipeline.call(miri_uncal_file,
                                                     steps=det1_dict,
                                                     output_dir=output_dir,
                                                     save_results=True,
                                                     save_calibrated_ramp=save_calibrated_ramp,
                                                     logcfg=logfile)


def miri_image2(miri_rate_files, output_dir,
                flatfile=None,
                photomfile=None,
                logfile=None):
    """
    Run CALWEBB_IMAGE2 pipeline on rate files.
    """
    # setup the directory of step paramters
    im2_dict = {}
    if flatfile is not None:
        im2_dict["flat_field"] = {"override_flat": flatfile}
    if photomfile is not None:
        im2_dict["photom"] = {"override_photom": f"./RefFiles/{photomfile}"}
    im2_dict["resample"] = {"pixfrac": 1.0}

    calwebb_image2.Image2Pipeline.call(miri_rate_files,
                                       steps=im2_dict,
                                       output_dir=output_dir,
                                       save_results=True,
                                       logcfg=logfile)


def miri_image3(miri_asn_file, output_dir, minobj=5, snr=5, fwhm=None,
                crval=None, crpix=None, rotation=None, output_shape=None,
                align_to_gaia=False, tweakreg=False,
                matchbkg=False, pixel_scale=0.11, sourcecat=True,
                logfile=None):
    """
    Run CALWEBB_IMAGE3 pipeline on asn file.
    """
    # setup the directory of step paramters
    im3_dict = {}

    # parameters based on Mattia's work

    sigma = 3.  # clipping limit, in sigma units, used when performing fit, default=3
    fit_geom = "shift"  # ftype of affine transformation to be considered when fitting catalogs, default='general'
    use2dhist = True  # boolean indicating whether to use 2D histogram to find initial offset, default=True

    if not tweakreg:
        im3_dict["tweakreg"] = {"skip": True}
    else:
        im3_dict["tweakreg"] = {"save_results": True,
                                "snr_threshold": snr,
                                "minobj": minobj,
                                "sigma": sigma,
                                "fitgeometry": fit_geom,
                                "use2dhist": use2dhist,
                                # "align_to_gaia": align_to_gaia,
                                "searchrad": 5.0}
        if fwhm is not None:
            im3_dict["tweakreg"]["kernel_fwhm"] = fwhm

    im3_dict["source_catalog"] = {"snr_threshold": snr}

    im3_dict["resample"] = {"kernel": "square",
                            "pixel_scale": pixel_scale}
    if crval is not None:
        im3_dict["resample"]["crval"] = crval
    if crpix is not None:
        im3_dict["resample"]["crpix"] = crpix
    if rotation is not None:
        im3_dict["resample"]["rotation"] = rotation
    if output_shape is not None:
        im3_dict["resample"]["output_shape"] = output_shape

    if not matchbkg:
        im3_dict["skymatch"] = {"skip": True}
    else:
        im3_dict["skymatch"] = {"lower": 0.0}

    if not sourcecat:
        im3_dict["source_catalog"] = {"skip": True}

    calwebb_image3.Image3Pipeline.call(miri_asn_file,
                                       steps=im3_dict,
                                       output_dir=output_dir,
                                       save_results=True,
                                       logcfg=logfile)
