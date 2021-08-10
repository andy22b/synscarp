"""
To deform DEM based on info from fault
"""
import rasterio
import rasterio.mask
from rasterio.enums import Resampling
import geopandas as gpd
import numpy as np
from affine import Affine
import rioxarray
import rioxarray.exceptions
from shapely.geometry import box
import os


def deform_dem(dem: str, out_tif: str, shifted_hw: gpd.GeoSeries, clipped_fw: gpd.GeoSeries, holes: gpd.GeoSeries,
               east_shift: float, north_shift: float, vertical_shift: float):
    """

    :param dem:
    :param out_tif:
    :param shifted_hw:
    :param clipped_fw:
    :param holes:
    :param east_shift:
    :param north_shift:
    :param vertical_shift:
    :return:
    """

    # Names of temporary files:
    temp_hw = "hw_shifted.tif"
    temp_fw = "clipped_fw.tif"
    assert os.path.exists(dem)
    # Extract data in hanging wall
    with rasterio.open(dem) as src:
        out_image = src.read(1)
        out_transform = src.transform
        out_meta = src.meta

    # To account for different structure in different tifs (makes sure positive north is up)
    if out_transform.e < 0.:
        north_shift *= -1.
    # Update transform with shift
    out_meta.update({"driver": "GTiff",
                     "transform": out_transform * Affine.translation(east_shift, north_shift)})

    out_image[out_image != 0.] += vertical_shift

    # Write shifted data to temporary tif for later handling with rioxarray (could probably improve)
    with rasterio.open(temp_hw, "w", **out_meta) as dest:
        dest.write(out_image, 1)

    # Extract data in clipped footwall polygon
    with rasterio.open(dem) as src:
        out_image, out_transform = rasterio.mask.mask(src, clipped_fw.geometry, crop=True)
        out_meta = src.meta

    # Update metadata to account for possible cropping
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})

    # Write to (another) temporary file
    with rasterio.open(temp_fw, "w", **out_meta) as dest:
        dest.write(out_image)

    # Read temporary tifs back in as x arrays
    hwx = rioxarray.open_rasterio(temp_hw)
    fwx = rioxarray.open_rasterio(temp_fw)
    # Reproject/interpolate shifted hanging wall onto same grid as footwall
    hwx_projected = hwx.rio.reproject_match(fwx, resampling=Resampling.bilinear)
    hwx.close()
    hwx_clipped = hwx_projected.rio.clip(shifted_hw.geometry)
    # hwx_clipped.rio.to_raster("hwx_clipped.tif")
    hwx_projected.close()
    hwx_padded = hwx_clipped.rio.pad_box(*fwx.rio.bounds())
    # hwx_padded.rio.to_raster("hwx_padded_early.tif")
    hwx_clipped.close()

    # Set NaNs to zero
    hwx_padded.data[np.abs(hwx_padded.data) > 10000.] = 0.

    # Combine shifted hangingwall and footwall datasets
    fwx.data[np.abs(fwx.data) > 10000.] = 0.
    hwx_padded.data += fwx.data
    fwx.close()
    # hwx_padded.rio.to_raster("hwx_padded.tif")

    # To deal with holes
    # Set insides of holes to nan
    clipped_projected = hwx_padded.rio.clip(holes.geometry, invert=True)
    hwx_padded.close()

    if holes is not None:
        # Empty array to store data from filled holes
        hole_filling = np.zeros(clipped_projected.shape)
        clipped_holes = gpd.clip(holes, box(*clipped_projected.rio.bounds()))
        for hole in list(clipped_holes.explode().geometry):
            # Cut to a box around hole of interest( might need to make box slightly wider in future

            clipped_hole = clipped_projected.rio.clip_box(*hole.bounds)
            # Fill hole using scipy linear interpolation
            filled_hole = clipped_hole.rio.interpolate_na(method="nearest")
            # Set data around hole to nans
            filled_hole_clipped = filled_hole.rio.clip([hole])
            # make same size as combined hw-fw dataset
            padded_hole = filled_hole_clipped.rio.pad_box(*clipped_projected.rio.bounds())
            # Set nans to zeros
            padded_hole.data[np.abs(padded_hole.data) > 10000.] = 0.
            # Add to storage dataset
            hole_filling += padded_hole.data
            padded_hole.close()
            filled_hole_clipped.close()
            filled_hole.close()
            clipped_hole.close()



        # Set nans to zero to add filled holes
        clipped_projected.data[clipped_projected.data > 10000.] = 0.
        clipped_projected.data += hole_filling
        # Set back to nans for writing
        clipped_projected.data[clipped_projected.data == 0.] = np.nan
    clipped_projected.rio.to_raster(out_tif)
    clipped_projected.close()
    os.remove(temp_fw)
    os.remove(temp_hw)
