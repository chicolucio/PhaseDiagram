from phase_diagram.phase_diagram import PhaseDiagram
from src.point_in_curve import point_in_curve, point_in_function
import numpy as np
# from phase_diagram import ureg

water = PhaseDiagram('H2O')
ureg = water.ureg
Q_ = ureg.Quantity


@ureg.wraps((ureg.m, ureg.m), (ureg.m, ureg.m, None, ureg.m, ureg.m), strict=False)
def straight_line(x0, range_=0, points=1, a=1, b=0):
    x = np.linspace(x0, x0 + range_, points)
    y = a * x + b
    return x, y


@ureg.wraps(ureg.m, (ureg.m, ureg.m, None, ureg.m, ureg.m), strict=False)
def straight_line_y(x0, range_=0, points=1, a=1, b=0):
    x = np.linspace(x0, x0 + range_, points)
    y = a * x + b
    return y


def test_point_in_straight_line_one_point():
    assert point_in_curve((Q_('3 m'), Q_('3 m')), straight_line(3))


def test_point_not_in_straight_line_one_point():
    assert not point_in_curve((Q_('3 m'), Q_('4 m')), straight_line(3))


def test_point_in_straight_line_two_points():
    assert point_in_curve((Q_('3 m'), Q_('3 m')), straight_line(3, 4, 2))


def test_point_not_in_straight_line_two_points():
    assert not point_in_curve((Q_('3 m'), Q_('4 m')), straight_line(3, 4, 2))


def test_point_in_straight_line_two_points_not_in_return_tuple():
    assert not point_in_curve((Q_('3.1 m'), Q_('3.1 m')), straight_line(3, 4, 2))


def test_point_in_straight_line_function():
    assert point_in_function((Q_('3 m'), Q_('3 m')), straight_line_y)


def test_point_not_in_straight_line_function():
    assert not point_in_function((Q_('3 m'), Q_('4 m')), straight_line_y)


def test_point_in_straight_line_function_tolerance():
    assert point_in_function((Q_('3 m'), Q_('3.001 m')), straight_line_y)


def test_point_not_in_straight_line_function_tolerance():
    assert not point_in_function((Q_('3 m'), Q_('3.002 m')), straight_line_y)
