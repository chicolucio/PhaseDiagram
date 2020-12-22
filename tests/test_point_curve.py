from phase_diagram.phase_diagram import PhaseDiagram
from src.point_in_curve import point_in_curve, point_in_function
import numpy as np

water = PhaseDiagram('H2O')


def straight_line(x0, range_=0, points=1, a=1, b=0):
    x = np.linspace(x0, x0 + range_, points)
    y = a * x + b
    return x, y


def test_point_in_straight_line_one_point():
    assert point_in_curve((3, 3), straight_line(3))


def test_point_not_in_straight_line_one_point():
    assert not point_in_curve((3, 4), straight_line(3))


def test_point_in_straight_line_two_points():
    assert point_in_curve((3, 3), straight_line(3, 4, 2))


def test_point_not_in_straight_line_two_points():
    assert not point_in_curve((3, 4), straight_line(3, 4, 2))


def test_point_in_straight_line_two_points_not_in_return_tuple():
    assert not point_in_curve((3.1, 3.1), straight_line(3, 4, 2))


def test_point_in_straight_line_function():
    assert point_in_function((3, 3), straight_line)


def test_point_not_in_straight_line_function():
    assert not point_in_function((3, 4), straight_line)


def test_point_in_straight_line_function_tolerance():
    assert point_in_function((3, 3.001), straight_line)


def test_point_not_in_straight_line_function_tolerance():
    assert not point_in_function((3, 3.002), straight_line)
