from phase_diagram.phase_diagram import PhaseDiagram
from src.point_in_curve import point_in_function
import numpy as np

water = PhaseDiagram('H2O')
ureg = water.ureg
Q_ = ureg.Quantity


def test_supercritical_fluid():
    assert water.physical_state((Q_('700 K'), Q_('1e8 Pa'))) == 'supercritical fluid'

def test_gas():
    assert water.physical_state((Q_('700 K'), Q_('1e5 Pa'))) == 'gas'

def test_vapour_close_to_lv_curve():
    assert water.physical_state((Q_('500 K'), Q_('1e5 Pa'))) == 'vapour'

def test_vapour_close_to_sv_curve():
    assert water.physical_state((Q_('250 K'), Q_('1e1 Pa'))) == 'vapour'

# TODO repeat with CO2
