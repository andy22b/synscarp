# %% [markdown]
# Import modules requires for running

# %%
import geopandas as gpd
from synthetic_offsets.dem import deform_dem
from synthetic_offsets.fault import Fault

# %% [markdown]
# Supply names of DEM, fault trace shapefile (both NZTM) and fault source parameters

# %%
dem = "BU23.tif"
trace = "simplified_kakapo.geojson"

fault = Fault(trace=trace, dip=80., boundary=dem, slip=-10., rake=180., dip_direction=340.)

# %% [markdown]
# Create deformed DEM (may throw up a warning but should still run)

# %%
hw, fw = fault.footwall_and_hangingwall_polygons(fault.clipped_trace)
e, n, u = fault.slip_and_rake
hw_shifted, clipped_fw = fault.moved_footwall_hangingwall_polygons(fault.clipped_trace, x_shift=e, y_shift=n)
holes = gpd.GeoSeries(hw.difference(hw_shifted), crs=2193)
deform_dem(dem, "kakapo_scarp.tif", hw_shifted, clipped_fw, holes, e, n, u)

# %%
len(list(holes.geometry[0].geoms))

# %%
import geopandas as gpd
across_lines = gpd.read_file("across.geojson")
parallel_south = gpd.read_file("parallel_south.geojson")
parallel_north = gpd.read_file("parallel_north.geojson")

# %%
from synthetic_offsets.io.array_operations import profile_tiff_along_linestring

# %%
u

# %%
profile_tiff_along_linestring("kakapo_scarp.tif", across_lines.geometry.iloc[0], spacing=1)

# %%
from synthetic_offsets.io.array_operations import read_tiff, read_grid



# %%
line = across_lines.geometry.iloc[1]
x, y, z = read_tiff("kakapo_scarp.tif", window=line.bounds, make_y_ascending=True, nan_threshold=1.e4)

# %%
z

# %%
import numpy as np
from scipy.interpolate import RegularGridInterpolator
distances = np.arange(0., line.length, 1.)
sample_points = [line.interpolate(d) for d in distances]
sample_x, sample_y = np.array([[point.x, point.y] for point in sample_points]).T

interp_func = RegularGridInterpolator((y, x), z, bounds_error=False, fill_value=np.nan)
profile_z = interp_func(np.array([sample_y, sample_x]).T)

from matplotlib import pyplot as plt
plt.plot(distances, profile_z)
plt.xlabel("Distance along line (m)")
plt.ylabel("Elevation (m)")
# plt.ylim(613,615)
plt.show()

# %%
read_grid("BU23.tif")

# %%
import rasterio
raster = rasterio.open("kakapo_scarp.tif")
raster.close()




