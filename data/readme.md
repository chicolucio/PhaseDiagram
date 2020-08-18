# Data file structure and info

## Column names and units

| Column           | Meaning                                  | Unit          |
| ---------------- | ---------------------------------------- | ------------- |
| 'Name'           | compound name                            | -             |
| 'Formula'        | compound formula                         | -             |
| 'CAS_number      | CAS number                               | -             |
| 'Molar_mass'     | compound molar mass                      | gram/mol      |
| 'd_solid'        | solid density                            | gram/cm^3     |
| 'd_liquid'       | liquid density                           | gram/cm^3     |
| 'd_gas'          | gas density                              | gram/cm^3     |
| 'TP_temperature' | Triple point temperature                 | Kelvin        |
| 'TP_pressure'    | Triple point pressure                    | Pascal        |
| 'CP_temperature  | Critical point temperature               | Kelvin        |
| 'CP_pressure'    | Critical point pressure                  | Pascal        |
| 'CP_density'     | Critical point density                   | mol/liter     |
| 'BP_temperature' | Normal boiling point temperature         | Kelvin        |
| 'MP_temperature' | Normal melting point temperature         | Kelvin        |
| 'A'              | Antoine parameter A                      | mmHg, Celsius |
| 'B'              | Antoine parameter B                      | mmHg, Celsius |
| 'C'              | Antoine parameter C                      | mmHg, Celsius |
| 'Tmin'           | Minimum temperature for A, B and C       | Celsius       |
| 'Tmax'           | Maximum temperature for A, B and C       | Celsius       |
| 'H_melt'         | Enthalpy change of fusion                | kJ/mol        |
| 'S_melt'         | Entropy change of fusion                 | J/(mol * K)   |
| 'V_melt'         | Volume change on fusion                  | cm^3 / mol    |
| 'V_melt_calc'    | Calculated volume change on fusion       | cm^3 / mol    |
| 'H_vap'          | Standard enthalpy change of vaporization | kJ/mol        |
| 'H_vap_boil'     | Enthalpy change of vaporization at BP    | kJ/mol        |
| 'S_vap'          | Entropy change of vaporization           | J/(mol * K)   |
| 'V_vap'          | Volume change on vaporization            | cm^3 / mol    |
| 'V_vap_calc'     | Calculated volume change on vaporization | cm^3 / mol    |
| 'H_sub'          | Enthalpy change of sublimation           | kJ/mol        |
| 'S_sub'          | Entropy change of sublimation            | J/(mol * K)   |
| 'V_sub'          | Volume change on sublimation             | cm^3 / mol    |
| 'V_sub_calc'     | Calculated volume change on sublimation  | cm^3 / mol    |

## References

1. W. M. Haynes, CRC Handbook of Chemistry and Physics, 97th edition, 2017
2. John A. Dean, Lange's Handbook of Chemistry, 15th ed, 1999
3. NIST Chemistry WebBook, available at https://webbook.nist.gov/chemistry/
4. Carl L. Yaws, The Yaws Handbook of Vapor Pressure - Antoine Coefficients, 2nd edition, 2015