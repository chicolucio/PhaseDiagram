from phase_diagram.phase_diagram import PhaseDiagram
from src.point_in_curve import point_in_curve, point_in_function
import numpy as np
from functools import partial
# from phase_diagram import ureg

water = PhaseDiagram('H2O')
ureg = water.ureg
Q_ = ureg.Quantity

water_clapeyron_sv = partial(water._clapeyron_sv_lv, curve='sv')
water_clapeyron_lv = partial(water._clapeyron_sv_lv, curve='lv')

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


def test_point_water_01():
    assert point_in_function((Q_('400 K'),Q_('246493.814 Pa') ), water._antoine_lv)


def test_point_water_02():
    assert not point_in_function((Q_('400 K'),Q_('246493.813 Pa') ), water._antoine_lv)


def test_point_water_03():
    assert point_in_function((Q_('100 K'),Q_('1.64e-12 Pa')), water_clapeyron_sv)


def test_point_water_04():
    assert not point_in_function((Q_('100 K'),Q_('1.64e-3 Pa')), water_clapeyron_sv)


def test_point_water_05():
    assert point_in_function((Q_('200 K'),Q_('0.875 Pa')), water_clapeyron_lv)


def test_point_water_06():
    assert not point_in_function((Q_('200 K'),Q_('0.877 Pa')), water_clapeyron_lv)


def test_point_water_07():
    assert point_in_function((Q_('270 K'),Q_('44150144.527 Pa')), water._clapeyron_sl)


def test_point_water_08():
    assert not point_in_function((Q_('270 K'),Q_('44150144.529 Pa')), water._clapeyron_sl)
