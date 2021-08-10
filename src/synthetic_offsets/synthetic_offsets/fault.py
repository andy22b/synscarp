from typing import Union
import os
import numpy as np
from shapely.geometry import LineString, Polygon, box
from synthetic_offsets.utilities import extend_trace, split_polygon_by_line, calculate_dip_direction
from synthetic_offsets.utilities import reverse_bearing, normalize_bearing, move_polygon_xy
import geopandas as gpd
import rasterio


class Fault:
    """
    To handle fault trace
    """
    def __init__(self, trace: Union[LineString, str], dip: Union[float, int], boundary: Union[Polygon, str] = None, slip: Union[float, int] = None,
                 rake: Union[float, int] = None, dip_direction: float = None, flip_dip_direction: bool = False, crs=2193):
        self._trace, self._dip, self._dip_dir_str, self.dip_dir_float = (None,) * 4
        self._boundary, self._dip_direction = (None,) * 2
        self._clipped_trace, self._clipped_moved_trace, self._moved_trace = (None,) * 3

        print("Initializing Fault object...")
        if dip_direction is not None:
            self.dip_direction = normalize_bearing(dip_direction)
        else:
            if flip_dip_direction:
                self.dip_direction = reverse_bearing(calculate_dip_direction(trace))
            else:
                self.dip_direction = calculate_dip_direction(trace)

        assert isinstance(trace, (LineString, str))
        if isinstance(trace, str):
            assert os.path.exists(trace)
            gdf = gpd.read_file(trace)
            self.trace = list(gdf.geometry.explode())[0]
        else:
            self.trace = trace

        self.dip = dip
        if isinstance(boundary, str):
            assert os.path.exists(boundary)
            with rasterio.open(boundary) as src:
                self.boundary = box(*src.bounds)
        else:
            self.boundary = boundary
        self.clipped_trace = self.trace_within_boundary(self.trace)
        self.crs = crs

        self.rake = rake
        self.slip = slip



    @property
    def trace(self):
        return self._trace

    @trace.setter
    def trace(self, trace: LineString):
        self._trace = trace

    @property
    def boundary(self):
        return self._boundary

    @boundary.setter
    def boundary(self, poly: Polygon):
        assert isinstance(poly, Polygon)
        if self.trace is not None:
            assert any([poly.contains(self.trace), poly.intersects(self.trace)])
        self._boundary = poly

    @property
    def dip(self):
        return self._dip

    @dip.setter
    def dip(self, value: Union[float, int]):
        assert isinstance(value, (float, int))
        assert value > 0., "Dip must be positive"
        assert value < 90., "Dip must be shallower than 90, to avoid ambiguity about slip vector."

        self._dip = value

    @property
    def dip_direction(self):
        return self._dip_direction

    @dip_direction.setter
    def dip_direction(self, value):
        print("Dip direction (azimuth is {:.2f})".format(float(value)))
        print("Change by overriding or using flip_dip_direction option...")
        self._dip_direction = normalize_bearing(value)

    @property
    def strike(self):
        return normalize_bearing(self.dip_direction - 90.)
    
    @property
    def down_dip_vector(self):
        horizontal_component = np.cos(np.radians(self.dip))
        vertical_component = -1. * np.sin(np.radians(self.dip))
        dd_vec = np.array([horizontal_component * np.sin(np.radians(self.dip_direction)),
                           horizontal_component * np.cos(np.radians(self.dip_direction)),
                           vertical_component])
        return dd_vec

    @property
    def up_dip_vector(self):
        # udv = self.down_dip_vector[:]
        # udv[-1] = -1 / udv[-1]
        # return udv / np.linalg.norm(udv)
        return -1. * self.down_dip_vector
    
    @property
    def along_strike_vector(self):
        return np.array([np.sin(np.radians(self.dip_direction - 90.)),
                         np.cos(np.radians(self.dip_direction - 90.)), 0.])

    @property
    def clipped_trace(self):
        return self._clipped_trace
    
    def slip_rake_to_3d(self, slip: float, rake: float):
        ud_cont = slip * np.sin(np.radians(rake)) * self.up_dip_vector
        ss_cont = slip * np.cos(np.radians(rake)) * self.along_strike_vector
        return ud_cont + ss_cont
    
    @property
    def slip_and_rake(self):
        if not any([x is None for x in (self.slip, self.rake)]):
            return self.slip_rake_to_3d(self.slip, self.rake)
        else:
            return None, None, None

    @clipped_trace.setter
    def clipped_trace(self, trace: LineString):
        self._clipped_trace = self.trace_within_boundary(trace)

    @property
    def clipped_trace_array(self):
        return np.array(self.clipped_trace.coords)

    @property
    def moved_clipped_trace(self):
        return self.moved_clipped_trace

    def move_clipped_trace(self, east_shift: float, north_shift: float):
        coord_array = np.array(self.trace.coords)
        coord_array[:, 0] += east_shift
        coord_array[:, 1] += north_shift
        moved_trace = LineString(coord_array)
        self._moved_trace = moved_trace
        self._clipped_moved_trace = self.trace_within_boundary(moved_trace)
    
    def trace_within_boundary(self, trace: LineString):
        """
        :param trace:
        :return:
        """
        extended_trace = extend_trace(trace)
        trace_gdf = gpd.GeoDataFrame(geometry=[extended_trace])
        boundary_gdf = gpd.GeoDataFrame(geometry=[self.boundary])

        clipped_trace = gpd.clip(trace_gdf, boundary_gdf)
        if len(clipped_trace.geometry) > 1:
            raise ValueError("Part of fault trace leaves DEM boundary")
        else:
            return list(clipped_trace.geometry)[0]

    def footwall_and_hangingwall_polygons(self, clipped_trace: LineString):
        polygons = split_polygon_by_line(self.boundary, clipped_trace)
        assert len(polygons) == 2, "Something gone wrong splitting DEM into two"
        hwfw_dict = {}
        for i, poly in enumerate(polygons):
            hwfw_dict[i] = 0.
            diff = poly.exterior.difference(clipped_trace)
            diff_array = np.array(diff.coords)
            for coord in np.array(clipped_trace.coords):
                diff_vectors = diff_array - coord
                distances = np.dot(diff_vectors, self.down_dip_vector[:-1])
                hwfw_dict[i] += np.sum(distances)


        if hwfw_dict[0] > hwfw_dict[1]:
            hanging_wall = polygons[0]
            footwall = polygons[1]

        elif hwfw_dict[0] < hwfw_dict[1]:
            hanging_wall = polygons[1]
            footwall = polygons[0]
        else:
            raise ValueError("Cannot tell hanging wall and footwall apart!")

        return gpd.GeoSeries(hanging_wall, crs=self.crs), gpd.GeoSeries(footwall, crs=self.crs)

    def moved_footwall_hangingwall_polygons(self, clipped_trace: LineString, x_shift: float, y_shift: float):
        hw, fw = self.footwall_and_hangingwall_polygons(clipped_trace)
        hw_shifted = move_polygon_xy(list(hw.geometry)[0], x_shift, y_shift)
        overriding = hw_shifted.intersection(fw)
        clipped_fw = fw.difference(overriding)

        return gpd.GeoSeries(hw_shifted, crs=self.crs), gpd.GeoSeries(clipped_fw, crs=self.crs)







