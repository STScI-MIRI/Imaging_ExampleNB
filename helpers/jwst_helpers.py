import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, ManualInterval, SqrtStretch


def show_image(
    data_2d, vmin, vmax, xpixel=None, ypixel=None, title=None, dmap="binary",
):
    """Function to generate a 2D, log-scaled image of the data,
    with an option to highlight a specific pixel (with a red dot).

    Parameters
    ----------
    data_2d : numpy.ndarray
        Image to be displayed

    vmin : float
        Minimum signal value to use for scaling

    vmax : float
        Maximum signal value to use for scaling

    xpixel : int
        X-coordinate of pixel to highlight

    ypixel : int
        Y-coordinate of pixel to highlight

    title : str
        String to use for the plot title
    """
    norm = ImageNormalize(
        data_2d, interval=ManualInterval(vmin=vmin, vmax=vmax), stretch=SqrtStretch()
    )
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(data_2d, origin="lower", norm=norm, cmap=plt.get_cmap(dmap))

    if xpixel and ypixel:
        plt.plot(xpixel, ypixel, marker="o", color="red", label="Selected Pixel")

    fig.colorbar(im, label="DN")
    plt.xlabel("Pixel column")
    plt.ylabel("Pixel row")
    if title:
        plt.title(title)


def overlay_catalog(
    data_2d,
    catalog,
    flux_limit=0,
    vmin=0,
    vmax=10,
    title=None,
    units="MJy/str",
    dmap="binary",
):
    """Function to generate a 2D image of the data,
    with sources overlaid.

    data_2d : numpy.ndarray
        2D image to be displayed

    catalog : astropy.table.Table
        Table of sources

    flux_limit : float
        Minimum signal threshold to overplot sources from catalog.
        Sources below this limit will not be shown on the image.

    vmin : float
        Minimum signal value to use for scaling

    vmax : float
        Maximum signal value to use for scaling

    title : str
        String to use for the plot title

    units : str
        Units of the data. Used for the annotation in the
        color bar
    """
    norm = ImageNormalize(
        data_2d, interval=ManualInterval(vmin=vmin, vmax=vmax), stretch=SqrtStretch()
    )
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(data_2d, origin="lower", norm=norm, cmap=plt.get_cmap(dmap))

    for row in catalog:
        if row["aper_total_flux"].value > flux_limit:
            plt.plot(
                row["xcentroid"],
                row["ycentroid"],
                marker="o",
                markersize="3",
                color="red",
            )

    plt.xlabel("Pixel column")
    plt.ylabel("Pixel row")

    fig.colorbar(im, label=units)
    fig.tight_layout()
    plt.subplots_adjust(left=0.15)

    if title:
        plt.title(title)
