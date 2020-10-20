import sqlite3
import pandas as pd
import numpy as np
from scipy import constants
from . import ureg
import re
from collections import namedtuple
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


def interaction_table():
    names = d['names'].drop(['id', 'name_alt1', 'name_alt2', 'name_alt3'], axis=1)
    compound = d['compounds'].drop(['name'], axis=1)
    compound['name'] = names
    compounds = compound[['id', 'name', 'formula', 'cas', 'molar mass']]
    compound_table = compounds.style.hide_index()
    return compound_table


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
    compound_ident_tuple = list(d['compounds'].loc[(d['compounds']['id'] == compound_idx)].itertuples(index=False,
                                                                                                      name=None))[0]
    Identification = namedtuple("comp_id", ["index", "cas", "formula", "molar_mass", "name_ref"])
    return Identification(*compound_ident_tuple)


def compound_names(compound):
    compound_idx = compound_index(compound)
    compound_names_tuple = list(d['names'].loc[(d['names']['id'] == compound_idx)].itertuples(index=False,
                                                                                              name=None))[0]
    Names = namedtuple("names", ["index", "name", "alt_name1", "alt_name2", "alt_name3"])
    return Names(*compound_names_tuple)


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
        antoine_tuple = list(table.loc[:, ['t_min', 't_max', 'A', 'B', 'C']].itertuples(index=False,
                                                                                        name=None))[value_index]
        Antoine = namedtuple("antoine", ["Tmin", "Tmax", "A", "B", "C"])
        return Antoine(*antoine_tuple)
    except IndexError:
        print('Invalid compound')


def point_table(compound, point_name):
    compound_idx = compound_index(compound)
    try:
        return d[point_name].loc[d[point_name]['id'] == compound_idx]
    except KeyError:
        print('Invalid point name')


@ureg.wraps(('kelvin', 'pascal'), [None, None, None])
def _point(compound, point_name, value_index=0):
    table = point_table(compound, point_name)
    return list(table.loc[:, ['temperature', 'pressure']].itertuples(index=False, name=None))[value_index]


def point(compound, point_name, value_index=0):
    point_with_units = _point(compound, point_name, value_index)
    Point = namedtuple("point", ["temperature", "pressure"])
    return Point(*point_with_units)


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
def volume_change_fusion(compound, value_index=0, calc=True, calc_values_index=(0, 0)):
    compound_idx = compound_index(compound)
    if calc:
        d_sol = density(compound, 'solid', calc_values_index[0])
        d_liq = density(compound, 'liquid', calc_values_index[1])
        molar_mass = compound_identification(compound)[3] * ureg('gram/mole')
        return volume_change_fusion_calc(d_liq, d_sol, molar_mass)
    try:
        return list(d['v_melt'].loc[(d['v_melt']['id'] == compound_idx), 'value'])[value_index]
    except IndexError:
        print('Not valid. Try calculated value.')


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


class Plot:
    def __init__(self, x_unit, y_unit, ax=None, scale_log=True, legend=False, title=True, title_text=''):
        self.x_unit = x_unit
        self.y_unit = y_unit
        self.ax = ax
        self.scale_log = scale_log
        self.legend = legend
        self.title = title
        self.title_text = title_text

        if self.ax is None:
            fig, self.ax = plt.subplots(figsize=(10, 8), facecolor=(1.0, 1.0, 1.0))

    def plot_customization(self):

        linewidth = 2
        size = 12

        # grid and ticks settings
        self.ax.minorticks_on()
        self.ax.grid(b=True, which='major', linestyle='--',
                     linewidth=linewidth - 0.5)
        self.ax.grid(b=True, which='minor', axis='both',
                     linestyle=':', linewidth=linewidth - 1)
        self.ax.tick_params(which='both', labelsize=size + 2)
        self.ax.tick_params(which='major', length=6, axis='both')
        self.ax.tick_params(which='minor', length=3, axis='both')

        # labels and size
        self.ax.xaxis.label.set_size(size + 4)
        self.ax.yaxis.label.set_size(size + 4)
        # ax.title.set_fontsize(size+6)  # not working, don't know why...

        if self.scale_log:
            self.ax.set_yscale('log')
            self.ax.set_ylabel('log(Pressure / {:~P})'.format(ureg(self.y_unit).units))
        else:
            # setting the y-axis to scientific notation and
            # getting the order of magnitude
            self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            self.ax.yaxis.major.formatter._useMathText = True
            self.ax.figure.canvas.draw()  # Update the text
            order_magnitude = self.ax.yaxis.get_offset_text().get_text().replace('\\times', '')
            self.ax.yaxis.offsetText.set_visible(False)

            self.ax.set_ylabel('Pressure / ' + order_magnitude +
                          ' {:~P}'.format(ureg(self.y_unit).units))

        self.ax.set_xlabel('Temperature / {:~P}'.format(ureg(self.x_unit).units))

        if self.legend:
            self.ax.legend(loc='best', fontsize=14,
                           title=self.format_formula(), title_fontsize=14)

        if not self.title:
            pass
        elif self.title_text == '':
            # TODO: resolver o format_formula
            # self.ax.set_title('Calculated phase diagram - ' + self.format_formula(),
            #                   fontsize=18)
            pass
        else:
            self.ax.set_title(self.title_text, fontsize=18)

        return self.ax

    def plot_arrays(self, tuple_two_arrays, limit=None, **kwargs):
        x, y = tuple_two_arrays
        try:
            y = y[y < limit]
            x = x[:len(y)]
        except:
            pass
        self.ax.plot(x.to(self.x_unit), y.to(self.y_unit), **kwargs)
        self.plot_customization()

    def plot_point(self, tuple_point, **kwargs):
        self.ax.scatter(tuple_point[0].to(self.x_unit), tuple_point[1].to(self.y_unit), **kwargs)
        self.plot_customization()
