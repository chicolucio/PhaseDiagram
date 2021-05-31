from functools import partial

import numpy as np
from scipy import constants

from src.helpers import compound_index, compound_identification, compound_names, density_table, density, \
    antoine, point, enthalpy, volume_change_fusion
from src.plot import Plot
from src.point_in_curve import point_in_function
from . import ureg
import re
from collections import namedtuple
import matplotlib.pyplot as plt


gas_constant = constants.gas_constant * ureg.J/(ureg.mol*ureg.K)


class PhaseDiagram:
    def __init__(self, compound):
        """
        Instantiates a PhaseDiagram object
        Parameters
        ----------
        compound : str
            compound name, formula or CAS
        """
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
        self.ureg = ureg
        self.number_of_points = 100

    def __repr__(self):
        return f'{self.__class__.__name__}(name= {self.name}, CAS= {self.cas}, formula= {self.formula})'

    def __str__(self):
        return f'Phase diagram data for compound {self.name}, CAS {self.cas}, formula {self.formula}'

    def _clapeyron_sl(self, temperature):
        cte = self.enthalpy_fusion / self.volume_change_fusion
        return self.triple_point.pressure + cte * np.log(temperature / self.triple_point.temperature)

    def _clapeyron_sv_lv(self, temperature, curve):
        if curve == 'sv':
            cte = self.enthalpy_sublimation / gas_constant
        if curve == 'lv':
            cte = self.enthalpy_vaporization / gas_constant
        return self.triple_point.pressure * np.exp(cte * (1/self.triple_point.temperature - 1/temperature))

    def _antoine_lv(self, temperature):
        Tmin, Tmax, A, B, C = self.antoine_si
        right_side = A - (B / (C + temperature.magnitude))
        return 10**right_side * ureg.Pa

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

        if self.volume_change_fusion > 0:
            temp_range = -temp_range

        T_arr = np.linspace(self.triple_point.temperature.magnitude,
                            self.triple_point.temperature.magnitude-temp_range,
                            self.number_of_points) * ureg.kelvin
        P_arr = self._clapeyron_sl(T_arr)
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
        if self.triple_point.temperature.magnitude - temp_range < 0:
            temp_range = int(abs(self.triple_point.temperature.magnitude - temp_range))
        T_arr = np.linspace(self.triple_point.temperature.magnitude - temp_range,
                            self.triple_point.temperature.magnitude,
                            self.number_of_points) * ureg.K
        P_arr = self._clapeyron_sv_lv(T_arr, curve='sv')
        return T_arr, P_arr

    def clapeyron_lv(self):
        """Clausius-Clapeyron liquid-vapor line data

        Returns
        -------
        tuple
            Tuple of arrays (temperature, pressure)
        """
        T_arr = np.linspace(self.triple_point.temperature.magnitude,
                            self.critical_point.temperature.magnitude,
                            self.number_of_points) * ureg.kelvin
        P_arr = self._clapeyron_sv_lv(T_arr, curve='lv')
        return T_arr, P_arr

    @property
    def antoine_si(self):
        """
        A, B and C in SI units.
        Returns
        -------
        tuple
            Minimum and maximum temperature for A, B and C in SI units.
        """
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
            temperature array, pressure array
        """
        T_arr = np.linspace(self.triple_point.temperature.magnitude,
                            self.critical_point.temperature.magnitude, self.number_of_points) * ureg.K

        P_arr = self._antoine_lv(T_arr)

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

    def plot(self, ax=None, T_unit='K', P_unit='Pa', scale_log=True, legend=True, title=True, title_text='',
             clapeyron_lv=False):
        """
        Plots the phase diagram

        Parameters
        ----------
        ax : matplotlib axis, optional
            axis where the plot will be shown. If None, one will be created
        T_unit : str
            pint unit of the temperature
        P_unit : str
            pint unit of the pressure
        scale_log : bool, default=True
            if the y-axis will have a log scale
        legend : bool, default=True
            if a legend will be shown
        title : bool, default=True
            if the plot will have a title
        title_text : str, default=''
            title text
        clapeyron_lv : bool, default=False
            if the Clapeyron liquid-vapour curve will be plotted along the Antoine one
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 8), facecolor=(1.0, 1.0, 1.0))
        ax.set_axisbelow(True)

        if title_text == '':
            title_text = f'Calculated phase diagram - {self.format_formula()}'

        graph = Plot(ax=ax, x_label='Temperature', y_label='Pressure', x_unit=T_unit, y_unit=P_unit, legend=legend,
                     scale_log=scale_log, title=title, title_text=title_text)
        linewidth = 3
        marker_size = 100
        graph.plot_arrays(self.clapeyron_sl(), limit=self.critical_point.pressure, label='Clapeyron S-L',
                          linewidth=linewidth, zorder=1)
        graph.plot_arrays(self.clapeyron_sv(), label='Clapeyron S-V', linewidth=linewidth, zorder=1)
        graph.plot_arrays(self.antoine_lv(), label='Antoine L-V', linewidth=linewidth, zorder=1)

        if clapeyron_lv:
            graph.plot_arrays(self.clapeyron_lv(), label='Clapeyron L-V', linewidth=linewidth, linestyle='--', zorder=1)

        graph.plot_point(self.triple_point, label='Triple Point', color='red', s=marker_size, zorder=2)
        graph.plot_point(self.critical_point, label='Critical Point', color='purple', s=marker_size, zorder=2)

    @staticmethod
    def plot_custom(curves=None, points=None, ax=None, T_unit='K', P_unit='Pa', scale_log=True, legend=True,
                    title=True, title_text=''):
        """
        Allow the creation of a custom plot

        Parameters
        ----------
        curves : tuple
            tuple with two arrays (x array, y array) with pint units
        points : tuple
            tuple with two values (x value, y value) with pint units
        ax : matplotlib axis, optional
            axis where the plot will be shown. If None, one will be created
        T_unit : str
            pint unit of the temperature
        P_unit : str
            pint unit of the pressure
        scale_log : bool, default=True
            if the y-axis will have a log scale
        legend : bool, default=True
            if a legend will be shown
        title : bool, default=True
            if the plot will have a title
        title_text : str, default=''
            title text
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 8), facecolor=(1.0, 1.0, 1.0))
            ax.set_axisbelow(True)

        graph = Plot(ax=ax, x_label='Temperature', y_label='Pressure', x_unit=T_unit, y_unit=P_unit, legend=legend,
                     scale_log=scale_log, title=title,
                     title_text=title_text)
        if curves:
            for curve in curves:
                graph.plot_arrays(curve['data_tuple'], label=curve['label'], **curve['kwargs'])

        if points:
            for point in points:
                graph.plot_point(point['data_tuple'], label=point['label'], **point['kwargs'])

    def physical_state(self, point):
        """
        Returns the physical state for a given point(temperature, pressure)
        Parameters
        ----------
        point
            tuple with temperature and pressure

        Returns
        -------
        string
            physical state
        """
        state = ''
        clapeyron_sv = partial(self._clapeyron_sv_lv, curve='sv')

        # triple point temperature

        if point[0] == self.triple_point[0]:
            if point[1] < self.triple_point[1]:
                state = 'vapour'
            else:
                if self.volume_change_fusion < 0:
                    state = 'liquid'
                else:
                    state = 'solid'

        # critical point temperature

        elif point[0] == self.critical_point[0]:
            if point[1] < self.critical_point[1]:
                state = 'vapour'
            else:
                state = 'liquid'

        # point on curve

        elif point_in_function(point, self._antoine_lv):
            state = 'liquid-vapour curve'
        elif point_in_function(point, self._clapeyron_sl):
            state = 'solid-liquid curve'
        elif point_in_function(point, clapeyron_sv):
            state = 'solid-vapour curve'

        # regions

        elif point[0] > self.critical_point.temperature:
            if point[1] > self.critical_point.pressure:
                state = 'supercritical fluid'
            else:
                state = 'gas'
        elif (point[0] > self.triple_point.temperature) and (point[1] < self._antoine_lv(point[0])):
            state = 'vapour'
        elif (point[0] < self.triple_point.temperature) and (point[1] < self._clapeyron_sv_lv(point[0], curve='sv')):
            state = 'vapour'
        elif self.volume_change_fusion > 0:
            if (point[0] < self.triple_point.temperature) and (point[1] > self._clapeyron_sv_lv(point[0], curve='sv')):
                state = 'solid'
            else:
                state = 'liquid'
        elif self.volume_change_fusion < 0:
            if (point[0] < self.triple_point.temperature) and (point[1] > self._clapeyron_sv_lv(point[0], curve='sv')):
                state = 'solid'
            else:
                state = 'liquid'

        return state
