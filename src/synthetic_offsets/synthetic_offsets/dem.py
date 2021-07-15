from synthetic_offsets.io.array_operations import read_tiff
import rasterio
import geopandas


class Dem:
    def __init__(self, tiff: str, region_file: str, Fault):

