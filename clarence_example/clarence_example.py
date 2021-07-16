import geopandas as gpd
from synthetic_offsets.dem import deform_dem
from synthetic_offsets.fault import Fault

dem = "clarence2012.tif"
trace = "simplified_papatea.shp"

fault = Fault(trace=trace, dip=60., boundary=dem, slip=10., rake=45., dip_direction=270.)

hw, fw = fault.footwall_and_hangingwall_polygons(fault.clipped_trace)
e, n, u = fault.slip_and_rake
hw_shifted, clipped_fw = fault.moved_footwall_hangingwall_polygons(fault.clipped_trace, x_shift=e, y_shift=n)
holes = gpd.GeoSeries(hw.difference(hw_shifted), crs=2193)
deform_dem(dem, "papatea_scarp.tif", hw_shifted, clipped_fw, holes, e, n, u)
