from typing import Union
import numpy as np
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import polygonize, unary_union, linemerge
import geopandas as gpd


def bounds_tuple_to_polygon(bounds: tuple):
    assert len(bounds) == 4


def smallest_difference(value1, value2):
    """
    Finds smallest angle between two bearings
    :param value1:
    :param value2:
    :return:
    """
    abs_diff = abs(value1 - value2)
    if abs_diff > 180:
        smallest_diff = 360 - abs_diff
    else:
        smallest_diff = abs_diff

    return smallest_diff


def normalize_bearing(bearing: Union[float, int]):
    """
    change a bearing (in degrees) so that it is an azimuth between 0 and 360.
    :param bearing:
    :return:
    """
    while bearing < 0:
        bearing += 360.

    while bearing >= 360.:
        bearing -= 360.

    return bearing


def bearing_leq(value: Union[int, float], benchmark: Union[int, float], tolerance: Union[int, float] = 0.1):
    """
    Check whether a bearing (value) is anticlockwise of another bearing (benchmark)
    :param value:
    :param benchmark:
    :param tolerance: to account for rounding errors etc
    :return:
    """
    smallest_diff = smallest_difference(value, benchmark)
    if smallest_diff > tolerance:
        compare_value = normalize_bearing(value + smallest_diff)
        return abs(compare_value - normalize_bearing(benchmark)) <= tolerance
    else:
        return False


def bearing_geq(value: Union[int, float], benchmark: Union[int, float], tolerance: Union[int, float] = 0.1):
    """
    Check whether a bearing (value) is clockwise of another bearing (benchmark)
    :param value:
    :param benchmark:
    :param tolerance: to account for rounding errors etc
    :return:
    """
    smallest_diff = smallest_difference(value, benchmark)
    if smallest_diff > tolerance:
        compare_value = normalize_bearing(value - smallest_diff)
        return abs(compare_value - normalize_bearing(benchmark)) <= tolerance
    else:
        return False


def reverse_bearing(bearing: Union[int, float]):
    """
    180 degrees from supplied bearing
    :param bearing:
    :return:
    """
    assert isinstance(bearing, (float, int))
    assert 0. <= bearing <= 360.
    new_bearing = bearing + 180.

    # Ensure strike is between zero and 360 (bearing)
    return normalize_bearing(new_bearing)


def reverse_line(line: LineString):
    """
    Change the order that points in a LineString object are presented.
    Important for OpenSHA, I think
    :param line:
    :return:
    """
    assert isinstance(line, LineString)
    x, y = line.xy
    x_back = x[-1::-1]
    y_back = y[-1::-1]
    new_line = LineString([[xi, yi] for xi, yi in zip(x_back, y_back)])
    return new_line


def calculate_strike(line: LineString):
    """
    Calculate the strike of a shapely linestring object with coordinates in NZTM,
    then adds 90 to get dip direction. Dip direction is always 90 clockwise from strike of line.
    :param line: Linestring object
    :return:
    """
    # Get coordinates
    x, y = line.xy
    x, y = np.array(x), np.array(y)
    # Calculate gradient of line in 2D
    p = np.polyfit(x, y, 1)
    gradient = p[0]
    # Gradient to bearing
    bearing = 90.0 - np.degrees(np.arctan2(gradient, 1))
    bearing_vector = np.array([np.sin(np.radians(bearing)), np.cos(np.radians(bearing))])

    # Determine whether line object fits strike convention
    relative_x = x - x[0]
    relative_y = y - y[0]

    distances = np.matmul(np.vstack((relative_x, relative_y)).T, bearing_vector)
    num_pos = np.count_nonzero(distances >= 0)
    num_neg = np.count_nonzero(distances < 0)

    if num_neg > num_pos:
        bearing += 180.

    while bearing >= 360.:
        bearing -= 360.
    while bearing < 0.:
        bearing += 0.

    return bearing


def calculate_dip_direction(line: LineString):
    """
    Calculate the strike of a shapely linestring object with coordinates in NZTM,
    then adds 90 to get dip direction. Dip direction is always 90 clockwise from strike of line.
    :param line: Linestring object
    :return:
    """
    bearing = calculate_strike(line)

    dip_direction = bearing + 90.
    # Ensure strike is between zero and 360 (bearing)
    while dip_direction < 0:
        dip_direction += 360.

    while dip_direction >= 360.:
        dip_direction -= 360.

    return dip_direction


def two_point_vector(p1: tuple, p2: tuple):
    """
    Find vector from P1 to P2
    :param p1:
    :param p2:
    :return:
    """
    assert not all([p1[0] == p2[0], p1[1] == p2[1]]), "Points are the same!"
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]

    vector = np.array([dx, dy])
    normalized_vector = vector / np.linalg.norm(vector)
    return normalized_vector


def extend_two_point_line(p1: tuple, p2: tuple, distance: float):
    """
    Extend a line between points by a specified distance
    :param p1:
    :param p2:
    :param distance:
    :return: Point at end of new line
    """
    assert distance > 0., "Distance should be > 0!"
    vector = two_point_vector(p1, p2)
    return Point(p2[0] + vector[0] * distance, p2[1] + vector[1] * distance)


def extend_trace(trace: LineString, distance: float = 1.e5):
    """
    Extend line with two or more points at both ends
    :param trace:
    :param distance:
    :return:
    """
    coords = list(trace.coords)
    new_start = extend_two_point_line(coords[1], coords[0], distance)
    new_end = extend_two_point_line(coords[-2], coords[-1], distance)
    new_coords = [new_start] + coords + [new_end]
    return LineString(new_coords)


def move_line(line: LineString, azimuth: Union[float, int], distance: Union[float, int]):
    """

    :param line:
    :param azimuth:
    :param distance:
    :return:
    """
    azimuth_vector = np.array([np.sin(np.radians(azimuth)), np.cos(np.radians(azimuth))])
    coords = list(line.coords)
    new_coords = []
    for point in coords:
        new_point = Point(point[0] + azimuth_vector[0] * distance, point[1] + azimuth_vector[1] * distance)
        new_coords.append(new_point)

    return LineString(new_coords)


def move_line_xy(line: LineString, xshift: float, yshift: float):
    """
    Shift shapely polygon by (xshift, yshift)
    :param polygon:
    :param xshift:
    :param yshift:
    :return:
    """
    coord_array = np.array(line.coords)
    coord_array[:, 0] += xshift
    coord_array[:, 1] += yshift

    return LineString(coord_array)

def move_polygon(polygon:Polygon, azimuth: Union[float, int], distance: Union[float, int]):
    moved_poly = move_line(polygon.exterior, azimuth, distance)
    return Polygon(moved_poly)

# def move_polygon_gpd(polygon:Polygon, azimuth: Union[float, int], distance: Union[float, int], crs: int = 2193):

def move_polygon_xy(polygon: Polygon, xshift: float, yshift: float):
    moved_poly = move_line_xy(polygon.exterior, xshift, yshift)
    return gpd.GeoSeries(Polygon(moved_poly), crs=2193)




def between_moved_lines(line1: LineString, line2: LineString):
    """

    :param line1:
    :param line2:
    :return:
    """
    line2_reversed = reverse_line(line2)
    l1_coords = list(line1.coords)
    combined_coords = LineString(list(line1.coords) + list(line2_reversed.coords) + [l1_coords[0]])
    mls = unary_union(combined_coords)
    return polygonize(mls)


def split_polygon_by_line(polygon: Polygon, line: LineString):
    """

    :param polygon:
    :param line:
    :return:
    """
    merged = linemerge([polygon.boundary, line])
    borders = unary_union(merged)
    polygons = polygonize(borders)
    return list(polygons)
