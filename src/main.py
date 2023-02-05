# -*- coding: utf-8 -*-
# Parametric study of a turbojet engine

# Librairies
import math

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

class Simulation:
    """
    OPR # [1] Overall pressure ratio := Pt3/Pt2.
    TET # [K] Turbine inlet temperature := total temperature @ CC exit (Tt4).
    """

    def __init__(
        self,
        OPR:float,
        TET:float) -> None:
        """
        Class initialization method.
        """
        self.OPR = OPR
        self.TET = TET

    def pressure_temperature_exit_air_intake(
        self,
        Pt1: float,
        Tt1: float) -> float:
        """
        Pt2,Tt2.
        """
        Pt2 = XI_E*Pt1
        Tt2 = Tt1*XI_E**((GAMMA - 1)/GAMMA)
        return [Pt2,Tt2]
    
    def pressure_temperature_exit_compressor(
        self,
        Pt2: float,
        Tt2: float) -> float:
        Pt3 = self.OPR*Pt2
        Tt3 = Tt2*self.OPR*((GAMMA - 1)/(GAMMA*ETA_POL_COMP))
        return [Pt3,Tt3]
    
    def pressure_temperature_exit_CC(
        self,
        Pt3: float,
        Tt3: float) -> float:
        Pt4 = XI_C*Pt3
        Tt4 = self.TET
        return [Pt4,Tt4]
    



# def pression_totale(P,M):
#     return P*(1 + (GAMMA - 1)/2*M**2)**(GAMMA/(GAMMA - 1))

# def temperature_totale(T,M):
#     return T*(1 + (GAMMA - 1)/2*M**2)

# def vitesse_vol(T,M):
#     return M*math.sqrt(GAMMA*r*T)

def pression_temperature_sortie_compresseur(P0,T0,M0):
    Pt3 = PI*xi*pression_totale(P0,M0)
    Tt3 = temperature_totale(T0,M0)*PI**((GAMMA - 1)/(GAMMA*eta_pol))
    return [Pt3,Tt3]

def pression_temperature_sortie_CC(Pt3,Tt3):
    Tt4 = 1600
    Pt4 = Pt3 * xi_c
    return [Pt4,Tt4]

def alpha(Tt4,Tt3):
    return (Cp_c*Tt4-Cp*Tt3)/(eta_comb*Pk-Cp_c*Tt4)

def pression_temperature_sortie_turbine(Pt4,Tt4,Tt3,Tt0):
    Tt5 = Tt4 - (Tt3 - Tt0)/((1+alpha(Tt4,Tt3) ) * eta_meca)
    Pt5 = Pt4*(Tt5/Tt4)**((GAMMA/((GAMMA - 1))*eta_pol))
    return [Pt5,Tt5]

def m9(Pt9,P9):
    return math.sqrt(2/(GAMMA - 1)*((Pt9/P9)**((GAMMA - 1)/GAMMA) - 1))

def pression_temperature_sortie_tuyere(Pt5,Tt5):
    return [xi*Pt5,Tt5]

def vitesse_son(T):
    return math.sqrt(GAMMA*r*T)

def t9(Pt9,Tt9,P9):
    return Tt9/(1 + (GAMMA - 1)/2*m9(Pt9,P9)**2)

def V9(Pt9,Tt9,P9):
    return m9(Pt9,P9)*vitesse_son(t9(Pt9,Tt9,P9))

def F_m_point(V9,V0,ALPHA):
    return V9*(1 + ALPHA) - V0

def M_point(F_spe,F):
    return F/F_spe

def W_spec_chimique(ALPHA):
    return ALPHA*Pk

def W_spec_cylce(V9,V0,ALPHA):
    return (1+ALPHA)*1/2*V9**2 - 1/2*V0**2

def W_spec_pr(F_spec,V0):
    return F_spec*V0

def eta_prop(F_spec,V0,W_cycle):
    return F_spec*V0/W_cycle

def eta_th(Wcy,Wchim):
    return Wcy/Wchim

def eta_prop(Wpr,Wcy):
    return Wpr/Wcy

def eta(ETA_PROP,ETA_TH0):
    return ETA_PROP * ETA_TH0

def diametre(m_point,rho,V2):
    S = m_point/(rho * V2)
    D = 2*math.sqrt(S/(3.14159265*(1-(rapport)**2)))
    return D



#### Trace graphique ####

def cycle_entier(P0,T0,M0,PI,Tt4,T):

    eta_pol = 0.9
    xi = 0.52
    eta_meca = 0.98

    L3 = pression_temperature_sortie_compresseur(P0,T0,M0)

    L4 = pression_temperature_sortie_CC(L3[0],L3[1])

    ALPHA= alpha(L4[1],L3[1])

    print(ALPHA)

    L5 = pression_temperature_sortie_turbine(L4[0],L4[1],L3[1],temperature_totale(T0,M0))

    print(L5[0],L5[1])

    L9 = pression_temperature_sortie_tuyere(L5[0],L5[1])
    
    M9 = m9(L9[0],P0)

    T9 = t9(L9[0],L9[1],P0)

    print(M9)

    V9 = vitesse_vol(T9,M9)

    print(V9)

    V0 = vitesse_vol(T,M0)

    F_SPE = F_m_point(V9,V0,ALPHA)

    print(F_SPE)

    W_SPEC_CHIM = W_spec_chimique(ALPHA)
    
    W_SPEC_CYCLE = W_spec_cylce(V9,V0,ALPHA)
    
    W_SPEC_PR = W_spec_pr(F_SPE,V0)

    ETA_PROP = eta_prop(W_SPEC_PR,W_SPEC_CYCLE)
    
    ETA_TH = eta_th(W_SPEC_CYCLE,W_SPEC_CHIM)

    ETA = eta(ETA_PROP,ETA_TH)

    F_SPE = F_m_point(V9,V0,ALPHA)

    M_POINT = M_point(F_SPE,F)

    RHO = P0/(r * T)

    DIAMETRE = diametre(M_POINT,RHO,vitesse_vol(T,0.6))

    print("W_spe_cycle" + str(W_SPEC_CYCLE) +  " W_spe_prop" + str(W_SPEC_PR) + "  " + str(M_POINT) + " diametre" + str(DIAMETRE) + " V9" + str(V9))

    return [ETA,ETA_TH,ETA_PROP]


Liste = []
for i in range(20,60):
    a = cycle_entier(P0,T0,M0,i,Tt4,217)
    Liste.append(a[0])









