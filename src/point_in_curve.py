import numpy as np


def point_in_curve(point, tuples_of_points):
    return point in zip(*tuples_of_points)


def point_in_function(point, function, tolerance=0.001):
    return np.isclose(point[1], function(point[0])[1], atol=tolerance)
