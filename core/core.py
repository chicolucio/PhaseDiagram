import sqlite3
import pandas as pd
import pint

DB = 'data/data.db'

ureg = pint.UnitRegistry()


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


def state_index(state):
    try:
        return d['phys_states'].loc[d['phys_states']['state'] == state, 'id'].item()
    except ValueError:
        print('Invalid state')


@ureg.wraps('gram/cm**3', [None, None])
def density_value(compound, state):
    compound_idx = compound_index(compound)
    state_idx = state_index(state)
    try:
        return list(d['density'].loc[((d['density']['state'] == state_idx) & (d['density']['id'] == compound_idx)),
                                     'value'])[0]
    except IndexError:
        print('Invalid state or compound')


def antoine(compound):
    compound_idx = compound_index(compound)
    try:
        return list(d['antoine'].loc[(d['antoine']['id'] == compound_idx),
                                     ['t_min', 't_max', 'A', 'B', 'C']].itertuples(index=False, name=None))[0]
    except IndexError:
        print('Invalid compound')


@ureg.wraps('kelvin', [None])
def boiling_point(compound):
    compound_idx = compound_index(compound)
    try:
        return list(d['boiling_point'].loc[(d['boiling_point']['id'] == compound_idx),
                                           'temperature'])[0]
    except IndexError:
        print('Invalid compound')


@ureg.wraps('kelvin', [None])
def melting_point(compound):
    compound_idx = compound_index(compound)
    try:
        return list(d['melting_point'].loc[(d['melting_point']['id'] == compound_idx),
                                           'temperature'])[0]
    except IndexError:
        print('Invalid compound')


@ureg.wraps(('kelvin', 'pascal'), [None])
def triple_point(compound):
    compound_idx = compound_index(compound)
    try:
        return list(d['triple_point'].loc[d['triple_point']['id'] == compound_idx,
                                           ['temperature', 'pressure']].itertuples(index=False, name=None))[0]
    except IndexError:
        print('Invalid compound')


@ureg.wraps(('kelvin', 'pascal'), [None])
def critical_point(compound):
    compound_idx = compound_index(compound)
    try:
        return list(d['critical_point'].loc[d['critical_point']['id'] == compound_idx,
                                            ['temperature', 'pressure']].itertuples(index=False, name=None))[0]
    except IndexError:
        print('Invalid compound')


@ureg.wraps('kJ/mole', [None])
def enthalpy_fusion(compound):
    compound_idx = compound_index(compound)
    try:
        return list(d['h_melt'].loc[(d['h_melt']['id'] == compound_idx), 'value'])[0]
    except IndexError:
        print('Invalid compound')


@ureg.wraps('kJ/mole', [None])
def enthalpy_sublimation(compound):
    compound_idx = compound_index(compound)
    try:
        return list(d['h_sub'].loc[(d['h_sub']['id'] == compound_idx), 'value'])[0]
    except IndexError:
        print('Invalid compound')


@ureg.wraps('kJ/mole', [None])
def enthalpy_vaporization(compound):
    compound_idx = compound_index(compound)
    try:
        return list(d['h_vap_boil'].loc[(d['h_vap_boil']['id'] == compound_idx), 'value'])[0]
    except IndexError:
        print('Invalid compound')


@ureg.wraps('(cm**3)/mole', [None])
def volume_change_fusion(compound):
    compound_idx = compound_index(compound)
    try:
        return list(d['v_melt'].loc[(d['v_melt']['id'] == compound_idx), 'value'])[0]
    except IndexError:
        print('Invalid compound')


class Data:
    def __init__(self, compound):
        self.compound = compound
        self.idx = compound_index(self.compound)
        self.density_solid = density_value(self.compound, 'solid')
        self.density_liquid = density_value(self.compound, 'liquid')
        self.antoine = antoine(self.compound)
        self.boiling_point = boiling_point(self.compound)
        self.melting_point = melting_point(self.compound)
        self.triple_point = triple_point(self.compound)
        self.critical_point = critical_point(self.compound)
        self.enthalpy_fusion = enthalpy_fusion(self.compound)
        self.enthalpy_sublimation = enthalpy_sublimation(self.compound)
        self.enthalpy_vaporization = enthalpy_vaporization(self.compound)
        self.volume_change_fusion = volume_change_fusion(self.compound)
