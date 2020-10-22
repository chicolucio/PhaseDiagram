import numpy as np
from scipy import constants

from src.helpers import compound_index, compound_identification, compound_names, density_table, density, \
    antoine, point, enthalpy, volume_change_fusion
from . import ureg
import re
from collections import namedtuple


gas_constant = constants.gas_constant * ureg.J/(ureg.mol*ureg.K)


class PhaseDiagram:
    def __init__(self, compound):
        self.compound = compound
        self.idx = compound_index(self.compound)
        self.cas = compound_identification(self.compound).cas
        self.formula = compound_identification(self.compound).formula
        self.molar_mass = compound_identification(self.compound).molar_mass * ureg('gram/mole')
        self.name = compound_names(self.compound).name
        self.alternative_names = (compound_names(self.compound).alt_name1,
                                  compound_names(self.compound).alt_name2,
                                  compound_names(self.compound).alt_name3)
        self.density_solid = density(self.compound, 'solid')
        self.density_liquid = density(self.compound, 'liquid')
        self.antoine = antoine(self.compound)
        self.boiling_point = point(self.compound, 'boiling_point')
        self.melting_point = point(self.compound, 'melting_point')
        self.triple_point = point(self.compound, 'triple_point')
        self.critical_point = point(self.compound, 'critical_point')
        self.enthalpy_fusion = enthalpy(self.compound, 'fusion')
        self.enthalpy_sublimation = enthalpy(self.compound, 'sublimation')
        self.enthalpy_vaporization = enthalpy(self.compound, 'vaporization')
        self.volume_change_fusion = volume_change_fusion(self.compound)
        self.density_table = density_table(self.compound)

    def __repr__(self):
        return f'{self.__class__.__name__}(name= {self.name}, CAS= {self.cas}, formula= {self.formula})'

    def __str__(self):
        return f'Phase diagram data for compound {self.name}, CAS {self.cas}, formula {self.formula}'

    def clapeyron_sl(self, temp_range=5):
        """Clausius-Clapeyron solid-liquid line data
        Parameters
        ----------
        temp_range : int, optional
            Temperature range around the triple point, by default 5
        Returns
        -------
        tuple
            Tuple of arrays (temperature, pressure)
        """
        # P(T) = P' + (H_melt / V_melt) ln(T / T') where T' is TP_temperature
        V_melt = self.volume_change_fusion

        if V_melt > 0:
            temp_range = -temp_range

        T_arr = np.linspace(self.triple_point.temperature.magnitude,
                            self.triple_point.temperature.magnitude-temp_range,
                            100) * ureg.kelvin
        cte = self.enthalpy_fusion / V_melt
        P_arr = self.triple_point.pressure + cte * np.log(T_arr / self.triple_point.temperature)
        return T_arr, P_arr

    def clapeyron_sv(self, temp_range=60):
        """Clausius-Clapeyron solid-vapor line data
        Parameters
        ----------
        temp_range : int, optional
            Temperature range around the triple point, by default 60
        Returns
        -------
        tuple
            Tuple of arrays (temperature, pressure)
        """
        # P(T) = P' exp[ (H_sub / R) (1 / T' - 1 / T) ] where T' is triple_point[0]
        T_arr = np.linspace(self.triple_point.temperature.magnitude - temp_range,
                            self.triple_point.temperature.magnitude,
                            100) * ureg.K
        cte = self.enthalpy_sublimation / gas_constant
        P_arr = self.triple_point.pressure * np.exp(cte * (1/self.triple_point.temperature - 1/T_arr))
        return T_arr, P_arr

    def clapeyron_lv(self):
        """Clausius-Clapeyron liquid-vapor line data
        Returns
        -------
        tuple
            Tuple of arrays (temperature, pressure)
        """
        # P(T) = P' exp[ (H_vap / R) (1 / T' - 1 / T) ] where T' is TP_temperature
        T_arr = np.linspace(self.triple_point.temperature.magnitude,
                            self.critical_point.temperature.magnitude,
                            100) * ureg.kelvin

        H_vap = self.enthalpy_vaporization

        cte = H_vap / gas_constant
        P_arr = self.triple_point.pressure * np.exp(cte * (1/self.triple_point.temperature - 1/T_arr))
        return T_arr, P_arr

    @property
    def antoine_si(self):
        Tmin = self.antoine.Tmin + 273.15
        Tmax = self.antoine.Tmax + 273.15
        A = self.antoine.A + np.log10(101325/760)
        B = self.antoine.B
        C = self.antoine.C - 273.15
        Antoine = namedtuple("antoine_si", ["Tmin", "Tmax", "A", "B", "C"])
        return Antoine(Tmin, Tmax, A, B, C)

    def antoine_lv(self):
        """Antoine liquid-vapor line data
        Returns
        -------
        tuple
            A, B and C for SI units. Temperature range (Tmin and Tmax) in Kelvin
            (temperature array, pressure array, A, B, C, Tmin, Tmax)
        """
        # log10(P) = A - (B / (C + T))
        T_arr = np.linspace(self.triple_point.temperature.magnitude,
                            self.critical_point.temperature.magnitude, 100) * ureg.K

        Tmin, Tmax, A, B, C = self.antoine_si

        right_side = A - (B / (C + T_arr.magnitude))

        P_arr = 10**right_side * ureg.Pa

        return T_arr, P_arr

    def format_formula(self):
        """ Display chemical formulas in a proper way
        Returns
        -------
        string
            LaTeX code to display chemical formulas in a proper way
        """
        label_formula = re.sub("([0-9])", "_\\1", self.formula)
        label_formula = r'$\mathregular{'+label_formula+'}$'
        return label_formula


