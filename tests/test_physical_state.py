from phase_diagram.phase_diagram import PhaseDiagram
from src.point_in_curve import point_in_function
import numpy as np

water = PhaseDiagram('H2O')
ureg_water = water.ureg
Q_water = ureg_water.Quantity

carbon_dioxide = PhaseDiagram('CO2')
ureg_carbon = carbon_dioxide.ureg
Q_carbon = ureg_carbon.Quantity


def test_supercritical_fluid():
    assert water.physical_state((Q_water('700 K'), Q_water('1e8 Pa'))) == 'supercritical fluid'

def test_gas():
    assert water.physical_state((Q_water('700 K'), Q_water('1e5 Pa'))) == 'gas'

def test_vapour_close_to_lv_curve():
    assert water.physical_state((Q_water('500 K'), Q_water('1e5 Pa'))) == 'vapour'

def test_vapour_close_to_sv_curve():
    assert water.physical_state((Q_water('250 K'), Q_water('1e1 Pa'))) == 'vapour'

def test_solid():
    assert water.physical_state((Q_water('250 K'), Q_water('1e3 Pa'))) == 'solid'

def test_liquid():
    assert water.physical_state((Q_water('400 K'), Q_water('1e7 Pa'))) == 'liquid'

def test_supercritical_fluid_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_carbon('350 K'), Q_carbon('1e10 Pa'))) == 'supercritical fluid'

def test_gas_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_carbon('350 K'), Q_carbon('1e6 Pa'))) == 'gas'

def test_vapour_close_to_lv_curve_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_carbon('260 K'), Q_carbon('1e6 Pa'))) == 'vapour'

def test_vapour_close_to_sv_curve_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_carbon('200 K'), Q_carbon('1e4 Pa'))) == 'vapour'

def test_solid_dioxide_carbon():
    assert carbon_dioxide.physical_state((Q_carbon('180 K'), Q_carbon('1e6 Pa'))) == 'solid'

def test_liquid_dioxide_carbon():
    assert carbon_dioxide.physical_state((Q_carbon('240 K'), Q_carbon('1e7 Pa'))) == 'liquid'

def test_liquid_vapour_balance_water():
    assert water.physical_state((Q_water('647.1 K'), Q_water('21936363.03431385 Pa'))) == 'liquid-vapour balance'

def test_solid_liquid_balance_water():
    assert water.physical_state((Q_water('268.16 K'), Q_water('70096120.44156641 Pa'))) == 'solid-liquid balance'

def test_solid_vapour_balance_water():
    assert water.physical_state((Q_water('213.16000000000003 K'), Q_water('2.619617119846615 Pa'))) == 'solid-vapour balance'

def test_liquid_vapour_balance_dioxide_carbon():
    assert carbon_dioxide.physical_state((Q_carbon('304.13 K'), Q_carbon('7371844.047303765 Pa'))) == 'liquid-vapour balance'

def test_solid_liquid_balance_dioxide_carbon():
    assert carbon_dioxide.physical_state((Q_carbon('221.58 K'), Q_carbon('23011236.108132083 Pa'))) == 'solid-liquid balance'

def test_solid_vapour_balance_dioxide_carbon():
    assert carbon_dioxide.physical_state((Q_carbon('156.58 K'), Q_carbon('2.619617119846615 Pa'))) == 'solid-vapour balance'
