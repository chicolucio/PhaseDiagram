from core.core import PhaseDiagram, state_index, compound_index
import pint


def test_water():
    water_data = PhaseDiagram('water')
    assert isinstance(water_data, PhaseDiagram)
    assert water_data.idx == 1
    assert water_data.name == 'water'
    assert water_data.cas == '7732-18-5'
    assert water_data.formula == 'H2O'
    assert water_data.molar_mass == pint.Quantity(18.0153, 'gram/mole')
    assert water_data.alternative_names == (None, None, None)
    assert isinstance(water_data.density_solid, pint.Quantity)
    assert water_data.density_solid == pint.Quantity(0.9167, 'gram/cm**3')
    assert water_data.density_liquid == pint.Quantity(0.9970474, 'gram/cm**3')
    assert water_data.antoine == (0.01, 373.98, 8.05573, 1723.6425, 233.08)
    assert water_data.boiling_point == (pint.Quantity(373.124, 'kelvin'), pint.Quantity(101325.0, 'pascal'))
    assert water_data.melting_point == (pint.Quantity(273.15, 'kelvin'), pint.Quantity(101325.0, 'pascal'))
    assert water_data.triple_point == (pint.Quantity(273.16, 'kelvin'), pint.Quantity(611.657, 'pascal'))
    assert water_data.critical_point == (pint.Quantity(647.1, 'kelvin'), pint.Quantity(2.206E7, 'pascal'))
    assert water_data.enthalpy_fusion == pint.Quantity(6.009, 'kJ/mole')
    assert water_data.enthalpy_sublimation == pint.Quantity(44.0, 'kJ/mole')
    assert water_data.enthalpy_vaporization == pint.Quantity(40.66, 'kJ/mole')
    assert water_data.volume_change_fusion == pint.Quantity(-1.634, 'cm**3/mole')

def test_carbon_dioxide():
    carbon_dioxide_data = PhaseDiagram('carbon dioxide')
    assert isinstance(carbon_dioxide_data, PhaseDiagram)
    assert carbon_dioxide_data.idx == 2
    assert carbon_dioxide_data.name == 'carbon dioxide'
    assert carbon_dioxide_data.cas == '124-38-9'
    assert carbon_dioxide_data.formula == 'CO2'
    assert carbon_dioxide_data.molar_mass == pint.Quantity(44.0095, 'gram/mole')
    assert carbon_dioxide_data.alternative_names == ('carbonic anhydride', None, None)
    assert isinstance(carbon_dioxide_data.density_solid, pint.Quantity)
    assert carbon_dioxide_data.density_solid == pint.Quantity(1.562, 'gram/cm**3')
    assert carbon_dioxide_data.density_liquid == pint.Quantity(1.179, 'gram/cm**3')
    assert carbon_dioxide_data.antoine == (-56.57, 31.04, 7.58828, 861.8205, 271.88)
    assert carbon_dioxide_data.boiling_point == (pint.Quantity(194.69, 'kelvin'), pint.Quantity(101325.0, 'pascal'))
    assert carbon_dioxide_data.melting_point == (pint.Quantity(216.59, 'kelvin'), pint.Quantity(101325.0, 'pascal'))
    assert carbon_dioxide_data.triple_point == (pint.Quantity(216.58, 'kelvin'), pint.Quantity(518500.0, 'pascal'))
    assert carbon_dioxide_data.critical_point == (pint.Quantity(304.13, 'kelvin'), pint.Quantity(7.375E6, 'pascal'))
    assert carbon_dioxide_data.enthalpy_fusion == pint.Quantity(9.02, 'kJ/mole')
    assert carbon_dioxide_data.enthalpy_sublimation == pint.Quantity(25.2, 'kJ/mole')
    assert carbon_dioxide_data.enthalpy_vaporization == pint.Quantity(15.8, 'kJ/mole')
    assert carbon_dioxide_data.volume_change_fusion == pint.Quantity(-1.672, 'cm**3/mole')


def test_state_index():
    assert state_index('solid') == 1
    assert state_index('liquid') == 2


def test_compound_index():
    assert compound_index('water') == 1
    assert compound_index('CO2') == 2
    assert compound_index('iodine') == 3
