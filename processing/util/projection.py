# Copyright 2017 Tom Eulenfeld, MIT license
import cartopy.crs as ccrs
import numpy as np
from osgeo import gdal


EPSG = 102003  # see https://epsg.io/102003

WKT = """PROJCS["USA_Contiguous_Albers_Equal_Area_Conic",
    GEOGCS["GCS_North_American_1983",
      DATUM["North_American_Datum_1983",
        SPHEROID["GRS_1980",6378137,298.257222101]],
      PRIMEM["Greenwich",0],
      UNIT["Degree",0.017453292519943295]],
    PROJECTION["Albers_Conic_Equal_Area"],
    PARAMETER["False_Easting",0],
    PARAMETER["False_Northing",0],
    PARAMETER["longitude_of_center",-96],
    PARAMETER["Standard_Parallel_1",29.5],
    PARAMETER["Standard_Parallel_2",45.5],
    PARAMETER["latitude_of_center",37.5],
    UNIT["Meter",1],
    AUTHORITY["EPSG","102003"]]
    """

PROJ4 = ('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 '
         '+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs')

AEA = ccrs.AlbersEqualArea(
    central_longitude=-96.0, central_latitude=37.5,
    false_easting=0.0, false_northing=0.0,
    standard_parallels=(29.5, 45.5),
    globe=ccrs.Globe(datum='NAD83', ellipse='GRS80'))

PC = ccrs.PlateCarree()


def scale_data(data):
    exp = 16
    max_ = 2**exp-2
    noval = 2**exp-1
    scale = (np.nanmax(data) - np.nanmin(data)) / max_
    offset = np.nanmin(data)
    dtype = np.uint16
    data = np.round((data - offset) / scale)
    data[np.isnan(data)] = noval
    data = np.array(data, dtype=dtype)
    return data, scale, offset, noval, dtype


def write_tif(fname, data, geotransform, description=''):
    data, scale, offset, noval, dtype = scale_data(data)
    dtype = {np.uint8: gdal.GDT_Byte, np.uint16: gdal.GDT_UInt16,
             np.int16: gdal.GDT_Int16}[dtype]
    # Byte, UInt16, Int16, UInt32, Int32, Float32, Float64
    (x, y) = data.shape
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(fname, y, x, 1, dtype)
    ds.SetMetadataItem('TIFFTAG_DOCUMENTNAME', fname.split('/')[-1])
    ds.SetMetadataItem('TIFFTAG_IMAGEDESCRIPTION', description)
    ds.SetMetadataItem('TIFFTAG_ARTIST', 'Tom Eulenfeld')
    copyr = ('Creative Commons Attribution 4.0 International License '
             '(CC BY 4.0), see http://creativecommons.org/licenses/by/4.0/')
    ds.SetMetadataItem('TIFFTAG_COPYRIGHT', copyr)
    ds.SetProjection(WKT)
    ds.SetGeoTransform(geotransform)
    band = ds.GetRasterBand(1)
    band.SetScale(scale)
    band.SetOffset(offset)
    band.SetNoDataValue(noval)
    band.WriteArray(data)
    # band.GetStatistics(0,1)
    ds = None


def load_tif(fname):
    ds = gdal.Open(fname)
    band = ds.GetRasterBand(1)
    zi = band.ReadAsArray()
    nan_ind = zi == band.GetNoDataValue()
    zi = zi * band.GetScale() + band.GetOffset()
    zi[nan_ind] = np.nan
    trans = ds.GetGeoTransform()
    extent = (trans[0], trans[0] + ds.RasterXSize*trans[1],
              trans[3] + ds.RasterYSize*trans[5], trans[3])
    ds = None
    return zi, extent


def add_scale(ax, length, loc, crs=ccrs.PlateCarree(), lw=1,
              cap=2, label=True, size=None, vpad=2):
    bx, by = ax.projection.transform_point(loc[0], loc[1], crs)
    bx1, bx2 = bx - 500 * length, bx + 500 * length
    ax.plot((bx1, bx2), (by, by), color='k', linewidth=lw)
    if cap:
        kw = {'xycoords': 'data', 'textcoords': 'offset points', 'arrowprops':
              {'arrowstyle': '-', 'connectionstyle': 'arc3',
               'shrinkA': 0, 'shrinkB': 0, 'linewidth': lw}}
        ax.annotate('', (bx1, by), (0, cap), **kw)
        ax.annotate('', (bx1, by), (0, -cap), **kw)
        ax.annotate('', (bx2, by), (0, cap), **kw)
        ax.annotate('', (bx2, by), (0, -cap), **kw)
    if label:
        ax.annotate(str(length) + ' km', (bx, by), (0, vpad), size=size,
                    textcoords='offset points', ha='center', va='bottom')
