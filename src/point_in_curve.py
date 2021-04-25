import numpy as np


def point_in_function(point, function, tolerance=0.001):
    """
    Verifies if a given point coordinates are in a given function domain/image

    Parameters
    ----------
    point : tuple of pint quantities
    function : function
    tolerance : float, optional

    Returns
    -------
    bool
    """
    return np.isclose(point[1].magnitude, function(point[0]).magnitude, atol=tolerance, rtol=0)
