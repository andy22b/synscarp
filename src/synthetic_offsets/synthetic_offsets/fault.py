from typing import Union
import numpy as np
from shapely.geometry import LineString, Polygon, Point
from synthetic_offsets.utilities import extend_trace, split_polygon_by_line, calculate_dip_direction, calculate_strike
from synthetic_offsets.utilities import reverse_bearing, normalize_bearing, move_polygon_xy
import geopandas as gpd





class Fault:
    """
    To handle fault trace
    """
    def __init__(self, trace: LineString, dip: Union[float, int], boundary: Polygon = None, slip: Union[float, int] = None,
                 rake: Union[float, int] = None, flip_dip_direction: bool = False):
        self._trace, self._dip, self._dip_dir_str, self.dip_dir_float = (None,) * 4
        self._boundary, self._dip_direction = (None,) * 2
        self._clipped_trace, self._clipped_moved_trace, self._moved_trace = (None,) * 3

        print("Initializing Fault object...")
        if flip_dip_direction:
            self.dip_direction = reverse_bearing(calculate_dip_direction(trace))
        else:
            self.dip_direction = calculate_dip_direction(trace)

        self.trace = trace
        self.dip = dip
        self.boundary = boundary

        self.clipped_trace = self.trace_within_boundary(self.trace)


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
    def clipped_trace(self):
        return self._clipped_trace

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

        return hanging_wall, footwall

    def moved_footwall_hangingwall_polygons(self, clipped_trace: LineString, x_shift: float, y_shift: float):
        hw, fw = self.footwall_and_hangingwall_polygons(clipped_trace)
        hw_shifted = move_polygon_xy(hw, x_shift, y_shift)
        return hw_shifted, fw







