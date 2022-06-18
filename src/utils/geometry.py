import numba
import numpy as np
from numba import jit, prange
from numba import types as t


def cal_boundary_length(length: int):
    return [[-length / 2.0, length / 2.0], [-length / 2.0, length / 2.0]]


def find_out_of_bounds(x_co_ordinates, y_co_ordinates, border_path):
    in_bound, out_bound = get_inbound_outbound(x_co_ordinates, y_co_ordinates, border_path)
    return out_bound


def get_inbound_outbound(x_co_ordinates, y_co_ordinates, border_path):
    bound_status = _get_bound_status(x_co_ordinates, y_co_ordinates, border_path)
    in_bound = np.where(bound_status == True)[0]
    out_bound = np.where(bound_status == False)[0]
    return in_bound, out_bound


# TODO: Use shapely
def _get_bound_status(x_co_ordinates, y_co_ordinates, border_path):
    return border_path.contains_points(np.vstack((x_co_ordinates, y_co_ordinates)).T)


@jit(t.UniTuple(t.float64, 2)(t.float64, t.float64))
def rand_cos_norm(norm_x, norm_y):
    """
    Returns a vector selected from a cosine distribution relative to an input normal vector.
    :param norm_x:
    :param norm_y:
    :return:
    """
    # An angle between -pi/2 and pi/2 following cosine distribution
    theta = np.pi * (np.random.rand() - 0.5)
    # TODO: ????
    while np.random.rand() > np.cos(theta):
        theta = np.pi * (np.random.rand() - 0.5)
    v_x_new = -norm_x * np.cos(theta) - norm_y * np.sin(theta)
    v_y_new = -norm_y * np.cos(theta) + norm_x * np.sin(theta)
    return v_x_new, v_y_new


# Expectation is len(x_co_ordinates) == len(y_co_ordinates)
@jit(t.UniTuple(t.float64[:], 3)(t.float64[:], t.float64[:]), nopython=True, parallel=True)
def poly_norms(x_co_ordinates, y_co_ordinates):
    dx = x_co_ordinates[1:] - x_co_ordinates[:-1]
    dy = y_co_ordinates[1:] - y_co_ordinates[:-1]

    lengths = np.sqrt(dx ** 2 + dy ** 2)
    norm_x = -dy / lengths
    norm_y = dx / lengths

    return norm_x, norm_y, lengths


@jit(t.UniTuple(t.float64, 2)(t.float64, t.float64, t.float64, t.float64))
def mirror_norm(norm_x, norm_y, v_x_in, v_y_in):
    vec_proj = norm_x * v_x_in + norm_y * v_y_in
    v_x_out = v_x_in - 2 * norm_x * vec_proj
    v_y_out = v_y_in - 2 * norm_y * vec_proj

    return v_x_out, v_y_out


@jit(t.boolean(t.float64, t.float64, t.float64, t.float64, t.float64, t.float64, t.float64, t.float64), nopython=True)
def seg_intersect(x1, y1, x2, y2, x3, y3, x4, y4):
    """
    Checks if two points [x3, y3] and [x4, y4] intersects a line/segment with endpoints [x1, y1] and [x2, y2].
    """
    x13 = x1 - x3
    y13 = y1 - y3
    x21 = x2 - x1
    y21 = y2 - y1
    x43 = x4 - x3
    y43 = y4 - y3

    ts = (x13 * y21 - y13 * x21) / (x43 * y21 - y43 * x21)
    us = (x13 * y43 - y13 * x43) / (x43 * y21 - y43 * x21)

    return (ts >= 0) and (ts <= 1) and (us >= 0) and (us <= 1)


@jit(t.int64[:](t.float64[:], t.float64[:], t.float64[:], t.float64[:]), nopython=True)
def seg_cross_poly(poly_x, poly_y, p1, p2):
    crossed = []
    total_edges = len(poly_x) - 1

    for edge_index in range(total_edges):
        if seg_intersect(poly_x[edge_index], poly_y[edge_index], poly_x[edge_index + 1],
                         poly_y[edge_index + 1], p1[0], p1[1], p2[0], p2[1]):
            crossed.append(edge_index)
    return np.array(crossed)


# TODO: Implement array-wise
# def find_crossed_edges(border_x, border_y, current_x, current_y, prev_x, prev_y):
#     x13 = border_x - current_x


# @jit(t.boolean(t.float64, t.float64, t.float64[:], t.float64[:]), nopython=True)
def check_bounds(x, y, poly_x, poly_y):
    n = len(poly_x)
    inside = False
    p2x = 0.0
    p2y = 0.0
    xints = 0.0
    p1x, p1y = poly_x[0], poly_y[0]
    for i in range(n + 1):
        p2x, p2y = poly_x[i % n], poly_y[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    # print(f"{p1x}, {p2x}, {x}")
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside


@jit(t.boolean(t.float64[:, :], t.float64[:]), nopython=True)
def is_inside_sm(polygon, point):
    length = len(polygon) - 1
    dy2 = point[1] - polygon[0][1]
    intersections = 0
    ii = 0
    jj = 1

    while ii < length:
        dy = dy2
        dy2 = point[1] - polygon[jj][1]

        # consider only lines which are not completely above/bellow/right from the point
        if dy * dy2 <= 0.0 and (point[0] >= polygon[ii][0] or point[0] >= polygon[jj][0]):

            # non-horizontal line
            if dy < 0 or dy2 < 0:
                F = dy * (polygon[jj][0] - polygon[ii][0]) / (dy - dy2) + polygon[ii][0]

                if point[0] > F:  # if line is left from the point - the ray moving towards left, will intersect it
                    intersections += 1
                elif point[0] == F:  # point on line
                    return 2

            # point on upper peak (dy2=dx2=0) or horizontal line (dy=dy2=0 and dx*dx2<=0)
            elif dy2 == 0 and (point[0] == polygon[jj][0] or (
                    dy == 0 and (point[0] - polygon[ii][0]) * (point[0] - polygon[jj][0]) <= 0)):
                return 2

            # there is another posibility: (dy=0 and dy2>0) or (dy>0 and dy2=0). It is skipped
            # deliberately to prevent break-points intersections to be counted twice.

        ii = jj
        jj += 1

    # print 'intersections =', intersections
    return intersections & 1


@jit(t.UniTuple(t.int64[:], 2)(t.float64[:], t.float64[:], t.float64[:], t.float64[:]), nopython=True, parallel=True)
def find_bound_status(x_co_ordinates, y_co_ordinates, poly_x, poly_y):
    status = np.empty(len(x_co_ordinates), dtype=numba.boolean)
    polygon = np.vstack((poly_x, poly_y)).T
    points = np.vstack((x_co_ordinates, y_co_ordinates)).T
    for i in prange(len(x_co_ordinates)):
        status[i] = is_inside_sm(polygon, points[i])
    in_bound = np.where(status)[0]
    out_bound = np.where(status == False)[0]
    return in_bound, out_bound
