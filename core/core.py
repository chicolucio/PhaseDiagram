import sqlite3
import pandas as pd
import numpy as np
from scipy import constants
from . import ureg, Q_
import re
import matplotlib.pyplot as plt

DB = 'data/data.db'

gas_constant = constants.gas_constant * ureg.J/(ureg.mol*ureg.K)


def database_dict(database):
    with sqlite3.connect(database) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")
        tables = cursor.fetchall()
        d = {}
        for t in tables:
            table_name = t[0]
            d[table_name] = pd.read_sql(f"select * from '{table_name}'", conn)
    return d


d = database_dict(DB)


def compound_index(compound):
    compound_name_idx = d['names'].loc[d['names'].isin([compound]).any(axis=1)].index.tolist()
    compound_formula_cas_idx = d['compounds'].loc[d['compounds'].isin([compound]).any(axis=1)].index.tolist()
    try:
        compound_idx = (compound_formula_cas_idx + compound_name_idx)[0]
    except IndexError as err:
        print(f'{str(err).capitalize()}. Not a valid compound.')
    else:
        return d['compounds'].loc[compound_idx, 'id']


def compound_identification(compound):
    compound_idx = compound_index(compound)
    return list(d['compounds'].loc[(d['compounds']['id'] == compound_idx)].itertuples(index=False, name=None))[0]


def state_index(state):
    try:
        return d['phys_states'].loc[d['phys_states']['state'] == state, 'id'].item()
    except ValueError:
        print('Invalid state')


def density_table(compound):
    compound_idx = compound_index(compound)
    try:
        return d['density'].loc[(d['density']['id'] == compound_idx)]
    except IndexError:
        print('Invalid state or compound')


@ureg.wraps('gram/cm**3', [None, None, None])
def density(compound, state, value_index=0):
    table = density_table(compound)
    state_idx = state_index(state)
    try:
        return list(table.loc[(table['state'] == state_idx), 'value'])[value_index]
    except IndexError:
        print('Invalid state or compound')


def antoine_table(compound):
    compound_idx = compound_index(compound)
    try:
        return d['antoine'].loc[(d['antoine']['id'] == compound_idx)]
    except IndexError:
        print('Invalid state or compound')


def antoine(compound, value_index=0):
    table = antoine_table(compound)
    try:
        return list(table.loc[:, ['t_min', 't_max', 'A', 'B', 'C']].itertuples(index=False, name=None))[value_index]
    except IndexError:
        print('Invalid compound')


def point_table(compound, point_name):
    compound_idx = compound_index(compound)
    try:
        return d[point_name].loc[d[point_name]['id'] == compound_idx]
    except KeyError:
        print('Invalid point name')


@ureg.wraps(('kelvin', 'pascal'), [None, None, None])
def point(compound, point_name, value_index=0):
    table = point_table(compound, point_name)
    return list(table.loc[:, ['temperature', 'pressure']].itertuples(index=False, name=None))[value_index]


def enthalpy_table(compound, enthalpy_name):
    compound_idx = compound_index(compound)
    name_dict = {'fusion': 'h_melt',
                 'sublimation': 'h_sub',
                 'vaporization': 'h_vap_boil'}
    try:
        return d[name_dict[enthalpy_name]].loc[d[name_dict[enthalpy_name]]['id'] == compound_idx]
    except KeyError:
        print('Invalid enthalpy name')


@ureg.wraps('kJ/mole', [None, None, None])
def enthalpy(compound, enthalpy_name, value_index=0):
    table = enthalpy_table(compound, enthalpy_name)
    return list(table.loc[:, 'value'])[value_index]


@ureg.wraps(None, ['gram/(cm**3)', 'gram/(cm**3)', 'gram/mol'])
def volume_change_fusion_calc(density_liquid, density_solid, molar_mass):
    return (1/density_liquid - 1/density_solid) * molar_mass


@ureg.wraps('(cm**3)/mole', [None, None, None, None])
def volume_change_fusion(compound, value_index=0, calc=True, calc_values_index=[0, 0]):
    compound_idx = compound_index(compound)
    if calc:
        d_sol = density(compound, 'solid', calc_values_index[0])
        d_liq = density(compound, 'liquid', calc_values_index[1])
        MM = compound_identification(compound)[3] * ureg('gram/mole')
        return volume_change_fusion_calc(d_liq, d_sol, MM)
    try:
        return list(d['v_melt'].loc[(d['v_melt']['id'] == compound_idx), 'value'])[value_index]
    except IndexError:
        print('Not valid. Try calculated value.')


class PhaseDiagram:
    def __init__(self, compound):
        self.compound = compound
        self.idx = compound_index(self.compound)
        self.cas = compound_identification(self.compound)[1]
        self.formula = compound_identification(self.compound)[2]
        self.molar_mass = compound_identification(self.compound)[3] * ureg('gram/mole')
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

        T_arr = np.linspace(self.triple_point[0].magnitude,
                            self.triple_point[0].magnitude-temp_range,
                            100) * ureg.kelvin
        cte = self.enthalpy_fusion / V_melt
        P_arr = self.triple_point[1] + cte * np.log(T_arr / self.triple_point[0])
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
        T_arr = np.linspace(self.triple_point[0].magnitude - temp_range,
                            self.triple_point[0].magnitude,
                            100) * ureg.K
        cte = self.enthalpy_sublimation / gas_constant
        P_arr = self.triple_point[1] * np.exp(cte * (1/self.triple_point[0] - 1/T_arr))
        return T_arr, P_arr

    def clapeyron_lv(self):
        """Clausius-Clapeyron liquid-vapor line data

        Returns
        -------
        tuple
            Tuple of arrays (temperature, pressure)
        """
        # P(T) = P' exp[ (H_vap / R) (1 / T' - 1 / T) ] where T' is TP_temperature
        T_arr = np.linspace(self.triple_point[0].magnitude,
                            self.critical_point[0].magnitude,
                            100) * ureg.kelvin

        H_vap = self.enthalpy_vaporization

        cte = H_vap / gas_constant
        P_arr = self.triple_point[1] * np.exp(cte * (1/self.triple_point[0] - 1/T_arr))
        return T_arr, P_arr

    def antoine_lv(self):
        """Antoine liquid-vapor line data

        Returns
        -------
        tuple
            A, B and C for SI units. Temperature range (Tmin and Tmax) in Kelvin
            (temperature array, pressure array, A, B, C, Tmin, Tmax)
        """
        # log10(P) = A - (B / (C + T))
        T_arr = np.linspace(self.triple_point[0].magnitude,
                            self.critical_point[0].magnitude, 100) * ureg.K

        A = self.antoine[2] + np.log10(101325/760)
        B = self.antoine[3]
        C = self.antoine[4] - 273.15
        Tmin = self.antoine[0] + 273.15
        Tmax = self.antoine[1] + 273.15

        right_side = A - (B / (C + T_arr.magnitude))

        P_arr = 10**right_side * ureg.Pa

        return T_arr, P_arr, A, B, C, Tmin, Tmax

    def format_formula(self):
        """ Display chemical formulas in a proper way

        Returns
        -------
        string
            LaTeX code to display chemical formulas in a proper way
        """
        label_formula = re.sub("([0-9])", "_\\1", self.formula)
        label_formula = '$\mathregular{'+label_formula+'}$'
        return label_formula

    def _plot_params(self, ax=None):
        """Internal function for plot parameters.

        Parameters
        ----------
        ax : Matplotlib axes, optional
            axes where the graph will be plotted, by default None
        """
        linewidth = 2
        size = 12

        # grid and ticks settings
        ax.minorticks_on()
        ax.grid(b=True, which='major', linestyle='--',
                linewidth=linewidth - 0.5)
        ax.grid(b=True, which='minor', axis='both',
                linestyle=':', linewidth=linewidth - 1)
        ax.tick_params(which='both', labelsize=size+2)
        ax.tick_params(which='major', length=6, axis='both')
        ax.tick_params(which='minor', length=3, axis='both')

        # labels and size
        ax.xaxis.label.set_size(size+4)
        ax.yaxis.label.set_size(size+4)
        # ax.title.set_fontsize(size+6)  # not working, don't know why...

        return

    def plot(self, parts=(1, 1, 0, 1), size=(10, 8), ax=None, T_unit='K',
             P_unit='Pa', scale_log=True, legend=False, title=True,
             title_text=''):
        """Plot function

        Parameters
        ----------
        parts : tuple, optional
            which lines will be plotted, by default (1, 1, 0, 1)
            By default, the solid-liquid, solid-vapor and liquid-vapor from
            Antoine equation lines are plotted. This can be changed with 0 and
            1's in a tuple`. 0 means turn off and 1 means turn on. The order in
            the tuple is:
            (solid-liquid Clausius-Clapeyron, solid-vapor Clausius-Clapeyron,
            liquid-vapor Clausius-Clapeyron, liquid-vapor Antoine)
        size : tuple, optional
            plot size, by default (10, 8)
        ax : Matplotlib axes, optional
            axes where the graph will be plotted, by default None
        T_unit : str, optional
            temperature unit, by default 'K'
        P_unit : str, optional
            pressure unit, by default 'Pa'
        scale_log : bool, optional
            logarithmic scale, by default True
        legend : bool, optional
            If a legend will be shown, by default False
        title : bool, optional
            If the plot will have a title, by default True
        title_text : str, optional
            Title text, by default ''

        Returns
        -------
        Matplotlib axes
            axes where the graph will be plotted
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=size, facecolor=(1.0, 1.0, 1.0))

        self._plot_params(ax)

        linewidth = 3.0

        if parts[0] == 1:
            T_clapeyron_sl, P_clapeyron_sl = self.clapeyron_sl()

            # ax.plot(T_clapeyron_sl.to(T_unit), P_clapeyron_sl.to(P_unit))
            # in order to avoid long SL lines, limit the pressure values to
            # those lower than self.CP_pressure
            P_clapeyron_sl = P_clapeyron_sl[P_clapeyron_sl < self.critical_point[1]]
            ax.plot(T_clapeyron_sl[:len(P_clapeyron_sl)].to(T_unit),
                    P_clapeyron_sl.to(P_unit),
                    'k-', label='SL boundary', linewidth=linewidth)

        if parts[1] == 1:
            T_clapeyron_sv, P_clapeyron_sv = self.clapeyron_sv()
            ax.plot(T_clapeyron_sv.to(T_unit),
                    P_clapeyron_sv.to(P_unit),
                    'b-', label='SV boundary', linewidth=linewidth)

        if parts[2] == 1:
            T_clapeyron_lv, P_clapeyron_lv = self.clapeyron_lv()
            ax.plot(T_clapeyron_lv.to(T_unit),
                    P_clapeyron_lv.to(P_unit),
                    'g--', label='LV boundary', linewidth=linewidth)

        if parts[3] == 1:
            T_antoine_lv, P_antoine_lv, *_ = self.antoine_lv()
            ax.plot(T_antoine_lv.to(T_unit),
                    P_antoine_lv.to(P_unit),
                    'r-', label='LV boundary - Antoine', linewidth=linewidth)

        if parts[2] == 1 or parts[3] == 1:
            ax.scatter(self.critical_point[0], self.critical_point[1],
                       s=100, label='Critical Point',
                       facecolors='orange', edgecolors='orange', zorder=3)

        ax.scatter(self.triple_point[0], self.triple_point[1],
                   s=100, label='Triple Point',
                   facecolors='m', edgecolors='m', zorder=3)

        if scale_log:
            ax.set_yscale('log')
            ax.set_ylabel('log(Pressure / {:~P})'.format(ureg(P_unit).units))
        else:
            # setting the y-axis to scientific notation and
            # getting the order of magnitude
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.yaxis.major.formatter._useMathText = True
            ax.figure.canvas.draw()  # Update the text
            order_magnitude = ax.yaxis.get_offset_text().get_text().replace('\\times', '')
            ax.yaxis.offsetText.set_visible(False)
            ax.set_ylabel('Pressure / ' + order_magnitude +
                          ' {:~P}'.format(ureg(P_unit).units))

        ax.set_xlabel('Temperature / {:~P}'.format(ureg(T_unit).units))

        if legend:
            ax.legend(loc='best', fontsize=14,
                      title=self.format_formula(), title_fontsize=14)

        if not title:
            pass
        elif title_text == '':
            ax.set_title('Calculated phase diagram - ' + self.format_formula(),
                         fontsize=18)
        else:
            ax.set_title(title_text, fontsize=18)

        return ax
