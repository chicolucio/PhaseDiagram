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
