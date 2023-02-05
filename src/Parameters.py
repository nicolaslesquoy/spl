"""
Notations used in this program:

CONSTANT = VALUE # [Units - metric system] Comment
[1] <=> no units.
"""

# ATMOSPHERIC CONDITIONS -> ISA Conditions @ sea level.
M0 = 0.78
T0 = 288 # [K]
P0 = 101325 # [Pa]

# CONSTANTS
GAMMA = 1.4 # [1] Heat capacity ratio.
GAMMA_C = 1.33 # [1] Heat capacity ratio of the mixture after the combustion.
R = 287.1 # [J.kg−1.K−1] *Specific* gas constant for dry air.
R_C = 291.6 # [J.kg−1.K−1] *Specific* gas constant for the mixture after the combustion.
CP = GAMMA*R/(GAMMA - 1) # [J.kg−1.K−1] *Specific* heat capacity at constant pressure.
CP_C = GAMMA_C*R_C/(GAMMA_C - 1) # [J.kg−1.K−1] *Specific* heat capacity at constant pressure of the mixture after the combustion.

# ENGINE CHARACTERISTICS
ETA_POL_COMP = 0.9 # [1] Polytropic efficiency of the compressor stage.
ETA_POL_TURB = 0.9 # [1] Polytropic efficiency of the turbine stage.
ETA_MECA = 0.98 # [1] Mechanical efficiency of the transmission.
ETA_COMB = 0.99 # [1] Combustion efficiency.
XI_E = 0.98 # [1] Pressure loss in the air-intake := Pt2/Pt1
XI_C = 0.95 # [1] Pressure loss inside the combustion chamber := Pt4/Pt3.
XI_T = 0.98 # [1] Pressure loss in the nozzle := Pt9/Pt5
PK = 40e6 # [J.kg-1] Lower heating value (LHV) of the fuel.
F = 50 # Targeted thrust value @ 0 ft.