# third-party imports
import numpy as np

# plotting imports
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.feature as cfeature

# cartopy local applications
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.io.shapereader import Reader
from cartopy.crs import PlateCarree


class GeoAxes:
    """
    This is a class for formatting and customizing geographical plots.
    """
    def __init__(self, ax: plt.Axes):
        """
        Initializes GeoAxes object.

        Args:
            ax (plt.Axes): Axes to the current figure.
        """
        self.ax = ax

    def geosettings(self,
                    xlocator=None,
                    ylocator=None,
                    left_labels=True,
                    right_labels=False,
                    top_labels=False,
                    bottom_labels=True,
                    xlines=True,
                    ylines=True,
                    ticks_size=20,
                    coastlines_color='#DDDDDD',
                    extent=None,
                    gridlines=True,
                    shapefile=None) -> None:
        """
        Set the settings for the geographical plot.

        Args:
            xlocator (array-like, optional, Defaults to None):
                Longitude tick locations.

            ylocator (array-like, optional, Defaults to None):
                Latitude tick locations.

            left_labels (bool, optional, Defaults to True):
                If True, place coordinates on the left-corner of the axes.

            right_labels (bool, optional, Defaults to False):
                If True, place coordinates on the right-corner of the axes.

            top_labels (bool, optional, Defaults to False):
                If True, annotate coordinates on the top of the axes

            bottom_labels (bool, optional, Defaults to True):
                If True, annotate coordinates at the bottom of the axes.

            xlines (bool, optional, Defaults to True):
                If True, add horizontal (longitudinal) gridlines.

            ylines (bool, optional, Defaults to True):
                If True, add vertical (latitudinal) gridlines.

            ticks_size (int, optional, Defaults to 20):
                Size of the coordinates ticks.

            coastlines_color (str, optional, Defaults to '#DDDDDD'):
                Color of the cartopy's coastlines contors.

            extent (list, optional, Defaults to None):
                Geographical extent of the map.
                If None, then the extent is automatically set from the dataset.
                If you want to set a specific extent, please pass a list
                with this format : ``[x0, x1, y0, y1]``, where:
                ``x0 ``is the most westerly gridpoint and ``x1`` is the most easterly one.
                ``y0`` is the northernmost gripoint and ``y1`` the southernmost one.

            gridlines (bool, optional):
                If True, add gridlines to the plot.

            shapefile (shapely Geometry, optional, Defaults to None):
                A collection of shapely geometries.
                Pass a shapefile to plot on axes.
        """
        if gridlines is True:
            # Adds gridlines to the axes, in the given coordinate system (crs)
            # In this case, the given crs is PlateCarree()
            gl = self.ax.gridlines(
                        crs=PlateCarree(),
                        draw_labels=True,
                        linewidth=0.5,
                        color='black',
                        alpha=0.3,
                        linestyle='--'
                    )
            # Configuring gridlines formatting
            gl.top_labels = top_labels
            gl.bottom_labels = bottom_labels
            gl.left_labels = left_labels
            gl.right_labels = right_labels
            gl.xlines = xlines
            gl.ylines = ylines

            # The Formatter to use for the x and y labels.
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER

            # Configuring x and y-ticks locating and formatting
            if xlocator is None:
                gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 10))
            else:
                gl.xlocator = mticker.FixedLocator(xlocator)
            if ylocator is None:
                gl.ylocator = mticker.FixedLocator(np.arange(-180, 180, 5))
            else:
                gl.ylocator = mticker.FixedLocator(ylocator)

            # Styling of the text labels
            gl.xlabel_style = {'size': ticks_size}
            gl.ylabel_style = {'size': ticks_size}

        # Add Cartopy's feature instances
        self.ax.coastlines()
        self.ax.add_feature(cfeature.LAND)
        self.ax.add_feature(cfeature.OCEAN)
        self.ax.add_feature(cfeature.COASTLINE, edgecolor=coastlines_color)
        try:
            self.ax.add_feature(cfeature.ShapelyFeature(
                                                Reader(shapefile).geometries(),
                                                crs=PlateCarree(),
                                                ls='-',
                                                edgecolor='#161616',
                                                facecolor='none'),
                                                lw=1.5
                                            )
        except TypeError:
            pass
        # Get the extent (x0, x1, y0, y1) of the map in the given crs
        if extent is None:
            self.ax.get_extent()
        # If no extent (x0, x1, y0, y1) is provided
        # Set the extent of the map automatically from axes
        else:
            self.ax.set_extent(extent)
