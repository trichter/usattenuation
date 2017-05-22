# Copyright 2017 Tom Eulenfeld, MIT license
import math

import numpy as np
from scipy.spatial import Delaunay
import shapely.geometry as geometry
from shapely.ops import cascaded_union, polygonize


# alpha_shape from
# http://nbviewer.ipython.org/github/dwyerk/boundaries/blob/master/
# concave_hulls.ipynb
def alpha_shape(points, alpha, allow_holes=True):
    """
    Compute the alpha shape (concave hull) of a set of points.

    @param points: Iterable container of points.
    @param alpha: alpha value to influence the gooeyness of the border. Smaller
                  numbers don't fall inward as much as larger numbers.
                  Too large, and you lose everything!
    """
    if len(points) < 4:
        return geometry.MultiPoint(list(points)).convex_hull

    def add_edge(edges, edge_points, coords, i, j):
        if (i, j) not in edges and (j, i) not in edges:
            edges.add((i, j))
            edge_points.append(coords[[i, j]])

    coords = np.array([point.coords[0] for point in points])
    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    hole_edges = set()
    hole_points = []
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
        a = math.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = math.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = math.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
        s = (a + b + c)/2.0
        area = math.sqrt(s*(s-a)*(s-b)*(s-c))
        circum_r = a*b*c/(4.0*area)
        if circum_r < 1.0/alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)
        else:
            add_edge(hole_edges, hole_points, coords, ia, ib)
            add_edge(hole_edges, hole_points, coords, ib, ic)
            add_edge(hole_edges, hole_points, coords, ic, ia)
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    if allow_holes:
        for tr in triangles[:]:
            if len(tr.boundary.coords.xy[0]) > 4:
                triangles.remove(tr)
    return cascaded_union(triangles), edge_points


def concave_hull(x, y, alpha, buf=0):
    points = [geometry.Point(*p) for p in zip(x, y)]
    hull, edge_points = alpha_shape(points, alpha=alpha)
    hull = hull.buffer(buf)
    return hull


def points_outside_network(x, y, xi, yi, alpha, buf=0):
    hull = concave_hull(x, y, alpha, buf=buf)
    index = np.zeros(shape=xi.shape, dtype=bool)
    it = np.nditer(index, flags=['multi_index'])
    while not it.finished:
        ind = it.multi_index
        index[ind] = not hull.contains(geometry.Point(xi[ind], yi[ind]))
        it.iternext()
    return index
