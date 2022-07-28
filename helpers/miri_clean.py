import copy
import warnings
import numpy as np


from astropy.modeling import models, fitting

from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.io import fits
from astropy.wcs import WCS

# from astropy.wcs.utils import proj_plane_pixel_scales
from regions import Regions
from astropy.coordinates import SkyCoord

from jwst import datamodels
from jwst.datamodels import dqflags
from tweakwcs import JWSTgWCS


def column_clean(mfile, exclude_above=None):
    """
    Remove the median of each column to suppress residual detector artifacts

    Requires rate, rateints, and cal versions of each file to exist

    Parameters
    ----------
    mfile : str
        filename with a MIRI rate image (i.e., xxx_rate.fits)
    exclude_above : float
        value above which to exclude data from calculating the column median
    """
    # Create kernel
    g = Gaussian1DKernel(stddev=20)

    # read in the final rate image
    rifile = mfile.replace("rate", "rateints")
    cfile = mfile.replace("rate", "wcs_cal").replace("stage1", "stage2")

    rdata = datamodels.open(mfile)
    ridata = datamodels.open(rifile)
    cdata = datamodels.open(cfile)

    # use the cal file dq flags as only after flat fielding are the outside the
    # FOV regions flagged
    bdata = cdata.dq & dqflags.pixel["DO_NOT_USE"] > 0

    # where to save the residual column data subtracted
    ridata_cols = copy.copy(ridata)

    colimage = np.zeros(rdata.data.shape)
    for k in range(ridata.data.shape[0]):
        cimage = copy.copy(ridata.data[k, :, :])
        # mask all the do_not_use data with NaNs
        cimage[bdata] = np.NaN
        # also remove zeros due to 2nd+ integration bug
        cimage[cimage == 0.0] = np.NaN
        # mask data above a threshold
        if exclude_above is not None:
            cimage[cimage > exclude_above] = np.NaN
        # compute the median of each column
        with warnings.catch_warnings():
            warnings.filterwarnings(
                action="ignore", message="All-NaN slice encountered"
            )
            colmeds = np.nanmedian(cimage, axis=0)
        # create a smoothed version to avoid removing large scale structure
        colmeds_smooth = convolve(colmeds - np.nanmedian(colmeds), g)
        # remove large scale structure from column medians
        colmeds_sub = colmeds - colmeds_smooth
        # make the 2D image version
        for j in range(cimage.shape[0]):
            colimage[j, :] = colmeds_sub
        # NaN all the no data pixels so they are not included in the median
        colimage[bdata] = np.NaN
        # subtarct the mean as we only want to remove residuals
        colimage -= np.nanmedian(colimage)
        # zero all the no data pixels
        colimage[bdata] = 0.0
        ridata.data[k, :, :] = ridata.data[k, :, :] - colimage
        ridata_cols.data[k, :, :] = colimage

    # median average new rateints to get the new rate image
    rdata.data = np.nanmedian(ridata.data, axis=0)

    # save the new rateints and rate results
    nfile = mfile.replace("rate.fits", "ccrateints.fits")
    ridata.save(nfile)

    nfile = mfile.replace("rate.fits", "ccrate.fits")
    rdata.save(nfile)

    # save the columns subtracted
    nfile = mfile.replace("rate.fits", "rateints_cols.fits")
    ridata_cols.save(nfile)


def row_clean(mfile, exclude_above=None):
    """
    Remove the median of each row to suppress residual detector artifacts

    Requires rate, rateints, and cal versions of each file to exist

    Parameters
    ----------
    mfile : str
        filename with a MIRI rate image (i.e., xxx_rate.fits)
    exclude_above : float
        value above which to exclude data from calculating the column median
    """
    # Create kernel
    g = Gaussian1DKernel(stddev=20)

    # read in the final rate image
    rifile = mfile.replace("ccrate", "ccrateints")
    cfile = mfile.replace("ccrate", "wcs_cal").replace("stage1", "stage2")

    rdata = datamodels.open(mfile)
    ridata = datamodels.open(rifile)
    cdata = datamodels.open(cfile)

    # use the cal file dq flags as only after flat fielding are the outside the
    # FOV regions flagged
    bdata = cdata.dq & dqflags.pixel["DO_NOT_USE"] > 0

    # where to save the residual column data subtracted
    ridata_rows = copy.copy(ridata)

    rowimage = np.zeros(rdata.data.shape)
    for k in range(ridata.data.shape[0]):
        rimage = copy.copy(ridata.data[k, :, :])
        # mask all the do_not_use data with NaNs
        rimage[bdata] = np.NaN
        # also remove zeros due to 2nd+ integration bug
        rimage[rimage == 0.0] = np.NaN
        # mask data above a threshold
        if exclude_above is not None:
            rimage[rimage > exclude_above] = np.NaN
        # exclude everything to the left of the imager FOV (basically the Lyot)
        rimage[:, 0:325] = np.NaN

        # compute the median of each column
        with warnings.catch_warnings():
            warnings.filterwarnings(
                action="ignore", message="All-NaN slice encountered"
            )
            rowmeds = np.nanmedian(rimage, axis=1)
        # create a smoothed version to avoid removing large scale structure
        rowmeds_smooth = convolve(rowmeds - np.nanmedian(rowmeds), g)
        # remove large scale structure from column medians
        rowmeds_sub = rowmeds - rowmeds_smooth
        # make the 2D image version
        for i in range(rimage.shape[1]):
            rowimage[:, i] = rowmeds_sub
        # NaN all the no data pixels so they are not included in the median
        rowimage[bdata] = np.NaN
        # subtarct the mean as we only want to remove residuals
        rowimage -= np.nanmedian(rowimage)
        # zero all the no data pixels
        rowimage[bdata] = 0.0
        ridata.data[k, :, :] = ridata.data[k, :, :] - rowimage
        ridata_rows.data[k, :, :] = rowimage

    # median average new rateints to get the new rate image
    rdata.data = np.nanmedian(ridata.data, axis=0)

    # save the new rateints and rate results
    nfile = mfile.replace("ccrate.fits", "cccrrateints.fits")
    ridata.save(nfile)

    nfile = mfile.replace("ccrate.fits", "cccrrate.fits")
    rdata.save(nfile)

    # save the columns subtracted
    nfile = mfile.replace("rate.fits", "rateints_rows.fits")
    ridata_rows.save(nfile)


def cal_column_clean(mfile, exclude_above=None):
    """
    Remove the median of each column to suppress residual detector artifacts

    works on cal images

    Parameters
    ----------
    mfile : str
        filename with a MIRI cal image (i.e., xxx_cal.fits)
    exclude_above : float
        value above which to exclude data from calculating the column median
    """
    # Create kernel
    g = Gaussian1DKernel(stddev=20)

    # read in the final rate image
    rdata = datamodels.open(mfile)
    rimage = copy.deepcopy(rdata.data)

    # use the cal file dq flags as only after flat fielding are the outside the
    # FOV regions flagged
    bdata = rdata.dq & dqflags.pixel["DO_NOT_USE"] > 0

    colimage = np.zeros(rimage.shape)

    # mask all the do_not_use data with NaNs
    rimage[bdata] = np.NaN

    # compute the median of each column
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", message="All-NaN slice encountered")
        colmeds = np.nanmedian(rimage, axis=0)
    # create a smoothed version to avoid removing large scale structure
    colmeds_smooth = convolve(colmeds - np.nanmedian(colmeds), g)
    # remove large scale structure from column medians
    colmeds_sub = colmeds - colmeds_smooth
    # make the 2D image version
    for j in range(rdata.shape[0]):
        colimage[j, :] = colmeds_sub
    # NaN all the no data pixels so they are not included in the median
    colimage[bdata] = np.NaN
    # subtarct the mean as we only want to remove residuals
    colimage -= np.nanmedian(colimage)
    # zero all the no data pixels
    colimage[bdata] = 0.0

    rdata.data -= colimage

    # save the new rateints and rate results
    nfile = mfile.replace("cal.fits", "cccal.fits")
    rdata.save(nfile)


def cal_row_clean(mfile, exclude_above=None):
    """
    Remove the median of each row to suppress residual detector artifacts

    works on cal images

    Parameters
    ----------
    mfile : str
        filename with a MIRI cal image (i.e., xxx_cal.fits)
    exclude_above : float
        value above which to exclude data from calculating the column median
    """
    # Create kernel
    g = Gaussian1DKernel(stddev=20)

    # read in the final cal image
    rdata = datamodels.open(mfile)

    # use the cal file dq flags as only after flat fielding are the outside the
    # FOV regions flagged
    bdata = rdata.dq & dqflags.pixel["DO_NOT_USE"] > 0

    rimage = copy.deepcopy(rdata.data)

    rowimage = np.zeros(rdata.data.shape)

    # mask all the do_not_use data with NaNs
    rimage[bdata] = np.NaN
    # also remove zeros due to 2nd+ integration bug
    rimage[rimage == 0.0] = np.NaN
    # mask data above a threshold
    if exclude_above is not None:
        rimage[rimage > exclude_above] = np.NaN
    # exclude everything to the left of the imager FOV (basically the Lyot)
    rimage[:, 0:325] = np.NaN

    # compute the median of each column
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", message="All-NaN slice encountered")
        rowmeds = np.nanmedian(rimage, axis=1)
    # create a smoothed version to avoid removing large scale structure
    rowmeds_smooth = convolve(rowmeds - np.nanmedian(rowmeds), g)
    # remove large scale structure from column medians
    rowmeds_sub = rowmeds - rowmeds_smooth
    # make the 2D image version
    for i in range(rimage.shape[1]):
        rowimage[:, i] = rowmeds_sub
    # NaN all the no data pixels so they are not included in the median
    rowimage[bdata] = np.NaN
    # subtarct the mean as we only want to remove residuals
    rowimage -= np.nanmedian(rowimage)
    # zero all the no data pixels
    rowimage[bdata] = 0.0

    rdata.data -= rowimage

    # save the new rateints and rate results
    nfile = mfile.replace("cccal.fits", "cccrcal.fits")
    rdata.save(nfile)


def flat_column_clean(rdata):
    """
    Remove the median of each column to suppress residual detector artifacts

    Assumes all bad data has already been NaNed
    Written to clean flat fields

    Parameters
    ----------
    rdata : 2D float ndarray
        image of flat field
    """
    # Create kernel
    g = Gaussian1DKernel(stddev=20)

    # where to save the residual column data subtracted
    rdata_cols = copy.copy(rdata)
    bdata = rdata == np.NaN

    colimage = np.zeros(rdata.data.shape)

    # compute the median of each column
    colmeds = np.nanmedian(rdata, axis=0)
    # create a smoothed version to avoid removing large scale structure
    colmeds_smooth = convolve(colmeds - np.nanmedian(colmeds), g)
    # remove large scale structure from column medians
    colmeds_sub = colmeds - colmeds_smooth
    # make the 2D image version
    for j in range(rdata.shape[0]):
        colimage[j, :] = colmeds_sub
    # NaN all the no data pixels so they are not included in the median
    colimage[bdata] = np.NaN
    # subtarct the mean as we only want to remove residuals
    colimage -= np.nanmedian(colimage)
    # zero all the no data pixels
    colimage[bdata] = 0.0
    rdata_cols -= colimage

    return rdata_cols


def fix_rateints_to_rate(mfile):
    """
    Fix the computation fo the rate file file the rateints file

    Issue with averaging the rateints file *and* in removing data with only
    one good measurement

    Parameters
    ----------
    mfile : str
        filename with a MIRI rate image (i.e., xxx_rate.fits)
    """
    # read in the final rate and rateint images
    rifile = mfile.replace("rate", "rateints")

    rdata = datamodels.open(mfile)
    ridata = datamodels.open(rifile)

    ridata.data[ridata.data == 0.0] = np.NaN
    # ridata.data[ridata.data > 8000.0] = np.NaN

    # median average new rateints to get the new rate image
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", message="Mean of empty slice")
        rdata.data = np.nanmean(ridata.data, axis=0)

    ndata = np.isnan(rdata.data)
    rdata.data[ndata] = 0.0
    rdata.dq[ndata] = 3
    ridata.data[np.isnan(ridata.data)] = 0.0

    # save the new rateints and rate results
    nfile = mfile.replace("rate.fits", "fixed_rateints.fits")
    ridata.save(nfile)

    nfile = mfile.replace("rate.fits", "fixed_rate.fits")
    rdata.save(nfile)


def shift_rate_wcs(mfile, shifts):
    """
    Shift the WCS of a rate file by the input shifts in arcsec

    Parameters
    ----------
    mfile : str
       filename with a MIRI rate image (i.e., xxx_rate.fits)
    shifts : [float, float]
       V2/V3 shifts in arcsec to apply
    """
    rate = fits.open(mfile)

    rate["SCI"].header["V2_REF"] += shifts[0]
    rate["SCI"].header["V3_REF"] += shifts[1]

    rate.writeto(
        mfile.replace("fixed_rate.fits", "fixed_wcs_rate.fits"), overwrite=True
    )
    rate.close()


def shift_cal_wcs(mfile, shifts):
    """
    Shift the WCS of a cal file by the input shifts in arcsec in V2, V3

    Parameters
    ----------
    mfile : str
       filename with a MIRI cal image (i.e., xxx_cal.fits)
    shifts : [float, float]
       V2/V3 shifts in arcsec to apply
    """
    image_model = datamodels.open(mfile)

    # no rotation or skew
    matrix = [[1.0, 0.0], [0.0, 1.0]]

    # create JWST WCS corrector object:
    wcs_corrector = JWSTgWCS(
        wcs=image_model.meta.wcs, wcsinfo=image_model.meta.wcsinfo.instance
    )
    wcs_corrector.set_correction(
        matrix=matrix,
        shift=shifts,
        ref_tpwcs=wcs_corrector,
    )
    image_model.meta.wcs = wcs_corrector.wcs

    image_model.write(mfile.replace("_fixed_cal.fits", "_fixed_wcs_cal.fits"))


def make_sky(
    files,
    subfiles=None,
    scalebkg=False,
    exclude_above=None,
    exclude_delta=None,
    ds9regions=None,
):
    """
    Make sky background by sigma clipping in image coordinates and subtract it
    from all the input files.

    Parameters
    ----------
    files : strs
       Array of cal files to use to create the sky image
    scalebkg : boolean
       Scale each image by its median to the average value [default=False]
    exclude_above : float
       Exclude data above this value from the sky creation
    exclude_delta : float
       Exclude data above the median bkg + this value from sky creation
    ds9regions : ds9 region file
       Exclude pixels inside ds9 regions from sky creation
    """
    if ds9regions is not None:
        ereg = Regions.read(ds9regions, format="ds9")
        # for creg in ereg:
        #     creg.radius *= 0.5

    istack = None
    for k, cfile in enumerate(files):
        print(f"processing {cfile}")
        cdata = datamodels.open(cfile)
        if istack is None:
            isize = cdata.data.shape
            istack = np.empty((isize[0], isize[1], len(files)))
            istackmed = np.empty((len(files)))
        tdata = cdata.data

        # remove all the non imager data
        # bdata = cdata.dq & dqflags.pixel["DO_NOT_USE"] > 0
        # tdata[bdata] = np.NaN

        if exclude_above is not None:
            tdata[tdata > exclude_above] = np.NaN

        if ds9regions is not None:
            # radeg, decdeg = cdata.meta.wcs([500., 600.], [500., 400.])
            # skycoord = SkyCoord(radeg, decdeg, unit='deg')
            # print(skycoord)
            # get standard WCS info from FITS header
            # with warnings.catch_warnings():
            #     warnings.simplefilter("ignore")
            #     t = fits.open(cfile)
            #     cwcs = WCS(t[1].header)
            # cwcs = cdata.meta.wcs.to_fits()

            fits_header, fits_hdulist = cdata.meta.wcs.to_fits()
            cwcs = WCS(fits_header)  # <-- "astropy" wcs

            pixx = np.arange(isize[1])
            pixy = np.arange(isize[0])
            imagex, imagey = np.meshgrid(pixx, pixy)
            imagera, imagedec = cwcs.wcs_pix2world(imagex, imagey, 0)
            # imagera, imagedec = cwcs.pixel_to_world(imagex, imagey, 0)
            skycoord = SkyCoord(imagera, imagedec, unit="deg")
            for creg in ereg:
                inoutimage = creg.contains(skycoord, cwcs)
                tdata[inoutimage] = np.NaN
            cdata.data = tdata
            cdata.write(cfile.replace("cal.fits", "cal_mask.fits"))
            # fits.writeto("test.fits", inoutimage * 1., overwrite=True)

        istackmed[k] = np.nanmedian(tdata)
        print(f"median sky = {istackmed[k]}")

        if exclude_delta is not None:
            tdata[tdata > istackmed[k] + exclude_delta] = np.NaN

        istack[:, :, k] = tdata

    # adjust the levels to the median
    # allows for data taken at different times with different backgrounds
    medsky = np.mean(istackmed)
    if scalebkg:
        for k in range(len(files)):
            istack[:, :, k] += medsky - istackmed[k]
            print(k, np.nanmedian(istack[:, :, k]))
    else:
        print("Not scaling individual images to median bkg")

    skyflat_mean, skyflat_median, skyflat_std = sigma_clipped_stats(
        istack, sigma_lower=3, sigma_upper=1, axis=2
    )

    # subtract the sky properly adjusted from the data
    if subfiles is None:
        subfiles = files
    for k, cfile in enumerate(subfiles):
        cdata = datamodels.open(cfile)
        cdata.data -= skyflat_mean
        if scalebkg:
            print(cfile, medsky - istackmed[k])
            cdata.data += medsky - istackmed[k]
        else:
            print(cfile)
        ndata = cdata.data == np.NaN
        cdata.data[ndata] = 0.0
        cdata.dq[ndata] = cdata.dq[ndata] & dqflags.pixel["DO_NOT_USE"]
        cdata.write(cfile.replace("_cal.fits", "_skysub_cal.fits"))

    return skyflat_mean


def get_subregion(pix, region, dsize, xwidth=60, ywidth=60, xoffset=20, yoffset=20):
    """
    Get the subregion for col/row pull up/down work

    Note that x,y are swapped from standard display due to python swapping them

    Parameters
    ----------
    pix : [float, float]
       pixel coordinates defining the center of col/row pull up/down
    region : str
       region to for which to create subregion [left, right, top, bottom]
    dsize : [int, int]
       image size in pixels
    xwidth : int
       x size of subregion
    ywidth : int
       y size of subregion
    xoffset : int
       x offset from pix to start region
    yoffset : int
       y offset from pix to start region

    Returns
    -------
    subregion : 4 tuple
       x1, x2, y1, y2 tuple giving the boundaries of the col/row affected region
    """
    if region == "bottom":
        x1 = 10
        x2 = pix[0] - xoffset
        y1 = int(max(pix[1] - 0.5 * ywidth, 0))
        y2 = int(min(pix[1] + 0.5 * ywidth, dsize[0]))
    elif region == "top":
        x1 = pix[0] + xoffset
        x2 = dsize[1] - 10
        y1 = int(max(pix[1] - 0.5 * ywidth, 0))
        y2 = int(min(pix[1] + 0.5 * ywidth, dsize[0]))
    elif region == "right":
        x1 = int(max(pix[0] - 0.5 * xwidth, 0))
        x2 = int(min(pix[0] + 0.5 * xwidth, dsize[1]))
        y1 = pix[1] + yoffset
        y2 = dsize[0] - 10
    elif region == "left":
        x1 = int(max(pix[0] - 0.5 * xwidth, 0))
        x2 = int(min(pix[0] + 0.5 * xwidth, dsize[1]))
        y1 = 10
        y2 = pix[1] - yoffset

    return x1, x2, y1, y2


def fix_rowcol_pull_updown(
    cfile, ds9_regfile, cortype=None, xwidth=60, ywidth=60, xoffset=20, yoffset=20
):
    """
    Correct the row/col pull up/down due to bright point sources

    Parameters
    ----------
    cfile : str
       MIRI cal file with col/row issues
    ds9_regfile : str
       ds9 region file that gives *point* regions at the peak of the sources
       causing row/col pull up/down
    cortype : array of strs
       on string for each point giving the regions to correct
       [None, all, topbotright, bottom]
    xwidth : int
       x size of subregion
    ywidth : int
       y size of subregion
    xoffset : int
       x offset from pix to start region
    yoffset : int
       y offset from pix to start region

    Returns
    -------
    corrected cal datamodel
    """
    cdata = datamodels.open(cfile)
    dsize = cdata.data.shape
    fits_header, fits_hdulist = cdata.meta.wcs.to_fits()
    cwcs = WCS(fits_header)  # <-- "astropy" wcs

    # get the points needing corrections
    regions = Regions.read(ds9_regfile, format="ds9")

    for k, creg in enumerate(regions):
        # print(f"correcting row/col pull up/down for source {k+1}")

        pix = np.flip(np.rint(cwcs.world_to_pixel(creg.center))).astype(int)

        dsize = cdata.data.shape

        # has to be on the main imager
        if (pix[0] >= 300) & (pix[0] < dsize[0]) & (pix[1] >= 0) & (pix[1] < dsize[1]):
            print("on image", cfile)

            if cortype is None:
                fixreg = ["bottom", "top", "right", "left"]
            else:
                if cortype[k] == "all":
                    fixreg = ["bottom", "top", "right", "left"]
                elif cortype[k] == "topbotright":
                    fixreg = ["bottom", "top", "right"]
                else:
                    fixreg = ["bottom"]
            for creg in fixreg:
                x1, x2, y1, y2 = get_subregion(
                    pix,
                    creg,
                    dsize,
                    xwidth=xwidth,
                    ywidth=ywidth,
                    xoffset=xoffset,
                    yoffset=yoffset,
                )

                # print(x1, x2, y1, y2)
                if (x1 < x2) & (y1 < y2):
                    if pix[0] < 300:
                        print(x1, x2, y1, y2)
                        print(pix)
                        print(cfile, "weird region")
                    # print("yes")
                    subimage = copy.copy(cdata.data[x1:x2, y1:y2])
                    if creg in ["bottom", "top"]:
                        caxis = 0
                    else:
                        caxis = 1
                        # remove all data outside the imager FOV
                        if creg == "left":
                            subimage[:, 0:350] = np.NaN
                    cvals = np.nanmedian(subimage, axis=caxis)

                    # fit a line
                    fit = fitting.LinearLSQFitter()
                    line_init = models.Linear1D()
                    x = np.arange(len(cvals))
                    y = cvals
                    fitted_line = fit(line_init, x, y)

                    # subtract it and make the image to subtract
                    diffy = y - fitted_line(x)
                    colimage = np.tile(diffy, (subimage.shape[caxis], 1))
                    if creg in ["right", "left"]:
                        colimage = np.transpose(colimage)

                    cdata.data[x1:x2, y1:y2] -= colimage

    cdata.write(cfile.replace(".fits", "_colrow.fits"))

    return cdata
