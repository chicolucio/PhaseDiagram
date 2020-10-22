from phase_diagram.phase_diagram import PhaseDiagram
import numpy as np


def test_water_clapeyron_antoine_array_methods():
    water = PhaseDiagram('water')
    assert np.allclose(water.clapeyron_sl(),
                       np.loadtxt("tests/golden_files/golden_water_clapeyron_sl.csv", delimiter=","))
    assert np.allclose(water.clapeyron_sv(),
                       np.loadtxt("tests/golden_files/golden_water_clapeyron_sv.csv", delimiter=","))
    assert np.allclose(water.clapeyron_lv(),
                       np.loadtxt("tests/golden_files/golden_water_clapeyron_lv.csv", delimiter=","))
    assert np.allclose(water.antoine_lv(),
                       np.loadtxt("tests/golden_files/golden_water_antoine_lv.csv", delimiter=","))


def test_co2_clapeyron_antoine_array_methods():
    co2 = PhaseDiagram('CO2')
    assert np.allclose(co2.clapeyron_sl(),
                       np.loadtxt("tests/golden_files/golden_co2_clapeyron_sl.csv", delimiter=","))
    assert np.allclose(co2.clapeyron_sv(),
                       np.loadtxt("tests/golden_files/golden_co2_clapeyron_sv.csv", delimiter=","))
    assert np.allclose(co2.clapeyron_lv(),
                       np.loadtxt("tests/golden_files/golden_co2_clapeyron_lv.csv", delimiter=","))
    assert np.allclose(co2.antoine_lv(),
                       np.loadtxt("tests/golden_files/golden_co2_antoine_lv.csv", delimiter=","))


def test_water_antoine_si():
    water = PhaseDiagram('water')
    assert water.antoine_si == (273.15999999999997, 647.13, 10.18063302013294, 1723.6425, -40.069999999999965)


def test_water_formula():
    water = PhaseDiagram('water')
    assert water.format_formula() == r'$\mathregular{H_2O}$'
