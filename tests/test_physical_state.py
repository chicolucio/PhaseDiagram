from phase_diagram.phase_diagram import PhaseDiagram
from phase_diagram.phase_diagram import ureg

Q_ = ureg.Quantity

water = PhaseDiagram('H2O')
carbon_dioxide = PhaseDiagram('CO2')


def test_supercritical_fluid():
    assert water.physical_state((Q_('700 K'), Q_('1e8 Pa'))) == 'supercritical fluid'


def test_gas():
    assert water.physical_state((Q_('700 K'), Q_('1e5 Pa'))) == 'gas'


def test_vapour_close_to_lv_curve():
    assert water.physical_state((Q_('500 K'), Q_('1e5 Pa'))) == 'vapour'


def test_vapour_close_to_sv_curve():
    assert water.physical_state((Q_('250 K'), Q_('1e1 Pa'))) == 'vapour'


def test_solid():
    assert water.physical_state((Q_('250 K'), Q_('1e3 Pa'))) == 'solid'


def test_liquid():
    assert water.physical_state((Q_('400 K'), Q_('1e7 Pa'))) == 'liquid'


def test_supercritical_fluid_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_('350 K'), Q_('1e10 Pa'))) == 'supercritical fluid'


def test_gas_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_('350 K'), Q_('1e6 Pa'))) == 'gas'


def test_vapour_close_to_lv_curve_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_('260 K'), Q_('1e6 Pa'))) == 'vapour'


def test_vapour_close_to_sv_curve_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_('200 K'), Q_('1e4 Pa'))) == 'vapour'


def test_solid_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_('180 K'), Q_('1e6 Pa'))) == 'solid'


def test_liquid_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_('240 K'), Q_('1e7 Pa'))) == 'liquid'


def test_liquid_vapour_balance_water():
    assert water.physical_state((Q_('643.32282828 K'), Q_('21056478.669068832 Pa'))) == 'liquid-vapour curve'


def test_solid_liquid_balance_water():
    assert water.physical_state((Q_('268.16 K'), Q_('70096120.44156641 Pa'))) == 'solid-liquid curve'


def test_solid_vapour_balance_water():
    assert water.physical_state((Q_('213.16000000000003 K'), Q_('2.619617119846615 Pa'))) == 'solid-vapour curve'


def test_liquid_vapour_balance_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_('300 K'), Q_('6733410.220019728 Pa'))) == 'liquid-vapour curve'


def test_solid_liquid_balance_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_('221.58 K'), Q_('23011236.108132083 Pa'))) == 'solid-liquid curve'


def test_solid_vapour_balance_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_('156.58 K'), Q_('2431.4640881597606 Pa'))) == 'solid-vapour curve'


def test_temperature_triple_point_low_pressure_water():
    assert water.physical_state((Q_('273.16 K'), Q_('1E2 Pa'))) == 'vapour'


def test_temperature_triple_point_high_pressure_water():
    assert water.physical_state((Q_('273.16 K'), Q_('1E9 Pa'))) == 'liquid'


def test_temperature_triple_point_low_pressure_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_('216.58 K'), Q_('1E2 Pa'))) == 'vapour'


def test_temperature_triple_point_high_pressure_carbon_dioxide():
    assert carbon_dioxide.physical_state((Q_('216.58 K'), Q_('1E7 Pa'))) == 'solid'


def test_temperature_critical_point_low_pressure_water():
    assert water.physical_state((Q_('647.1 K'), Q_('1E2 Pa'))) == 'vapour'


def test_temperature_critical_point_high_pressure_water():
    assert water.physical_state((Q_('647.1 K'), Q_('1E9 Pa'))) == 'liquid'
