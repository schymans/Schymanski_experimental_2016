
# coding: utf-8

# # Equations to derive leaf energy balance components from wind tunnel measurements and compare against leaf model</h1>
# 
# 

# In[1]:

get_ipython().run_cell_magic(u'capture', u'storage', u"# The above redirects all output of the below commands to the variable 'storage' instead of displaying it.\n# It can be viewed by typing: 'storage()'\n# Setting up worksheet and importing equations from other worksheets\nload('temp/Worksheet_setup.sage')\nload_session('temp/leaf_enbalance_eqs.sobj')")


# In[2]:

# Energy balance as a function of only input variables
eqenbalTl1 = (eq_Rs_enbal - R_s).rhs().subs(eq_El, eq_Hl, eq_Rll).subs(eq_Elmol).subs(eq_Cwl, eq_gtw).subs(eq_gbw_hc, eq_Pwl).subs(eq_hc, eq_rhoa_Pwa_Ta).subs(eq_Nu_forced_all, eq_Le).subs(eq_Re, eq_Dva, eq_ka).subs(eq_alphaa, eq_nua).subs(eq_PN2, eq_PO2)
eqenbalTl1.subs(cdict).args()


# In[3]:

# Energy balance as a function of variables that are independent of leaf temperature
eqenbalTl = (eq_Rs_enbal - R_s).rhs().subs(eq_El, eq_Hl, eq_Rll).subs(eq_Elmol).subs(eq_Cwl).subs(eq_Pwl)
eqenbalTl.subs(cdict).args()


# In[4]:

vdict = cdict.copy()
vdict[a_s] = 1.0    # one sided stomata
vdict[g_sw] = 0.01    
vdict[T_a] = 273 + 25.5
vdict[T_w] = vdict[T_a] # Wall temperature equal to air temperature
vdict[P_a] = 101325
rha = 1
vdict[P_wa] = rha*eq_Pwl.rhs()(T_l = T_a).subs(vdict)
vdict[L_l] = 0.03
#vdict[L_A] = vdict[L_l]^2
vdict[Re_c] = 3000
vdict[R_s] = 600
#vdict[Q_in] = 0
vdict[v_w] = 1

vdict[C_wa] = eq_Cwl.rhs()(P_wl = P_wa, T_l = T_a).subs(vdict)
vdict[nu_a] = eq_nua.rhs().subs(vdict)
vdict[Re] = eq_Re.rhs().subs(vdict)
vdict[Nu] = eq_Nu_forced_all.rhs().subs(vdict)
vdict[k_a] = eq_ka.rhs().subs(vdict)
vdict[h_c] = eq_hc.rhs().subs(vdict)
vdict[alpha_a] = eq_alphaa.rhs().subs(vdict)
vdict[D_va] = eq_Dva.rhs().subs(vdict)
vdict[Le] =  eq_Le.rhs().subs(vdict)
vdict[P_N2] = eq_PN2.rhs().subs(vdict)
vdict[P_O2] = eq_PO2.rhs().subs(vdict)
vdict[rho_a] =  eq_rhoa_Pwa_Ta.rhs().subs(vdict)
vdict[g_bw] = eq_gbw_hc.rhs().subs(vdict)
vdict[g_tw] =  eq_gtw.rhs().subs(vdict)


print eqenbalTl.subs(vdict).args()
vdict[T_l] = find_root(eqenbalTl.subs(vdict), 273, 373)
vdict[P_wl] = eq_Pwl.rhs().subs(vdict)
vdict[C_wl] = eq_Cwl.rhs().subs(vdict)
vdict[E_lmol] = eq_Elmol.rhs().subs(vdict)
eq_El.subs(vdict)


# In[5]:

def fun_SS(vdict1):
    '''
    Steady-state T_l, R_ll, H_l and E_l under forced conditions.
    Parameters are given in a dictionary (vdict) with the following entries:
    a_s, a_sh, L_l, P_a, P_wa, R_s, Re_c, T_a, g_sw, v_w
    ''' 
    vdict = vdict1.copy()

    vdict[C_wa] = eq_Cwl.rhs()(P_wl = P_wa, T_l = T_a).subs(vdict)
    vdict[nu_a] = eq_nua.rhs().subs(vdict)
    vdict[Re] = eq_Re.rhs().subs(vdict)
    vdict[Nu] = eq_Nu_forced_all.rhs().subs(vdict)
    vdict[k_a] = eq_ka.rhs().subs(vdict)
    vdict[h_c] = eq_hc.rhs().subs(vdict)
    vdict[alpha_a] = eq_alphaa.rhs().subs(vdict)
    vdict[D_va] = eq_Dva.rhs().subs(vdict)
    vdict[Le] =  eq_Le.rhs().subs(vdict)
    vdict[P_N2] = eq_PN2.rhs().subs(vdict)
    vdict[P_O2] = eq_PO2.rhs().subs(vdict)
    vdict[rho_a] =  eq_rhoa_Pwa_Ta.rhs().subs(vdict)
    vdict[g_bw] = eq_gbw_hc.rhs().subs(vdict)
    vdict[g_tw] =  eq_gtw.rhs().subs(vdict)
    
    vdict[T_l] = find_root(eqenbalTl.subs(vdict), 273, 373)
    vdict[P_wl] = eq_Pwl.rhs().subs(vdict)
    vdict[C_wl] = eq_Cwl.rhs().subs(vdict)
    vdict[E_lmol] = eq_Elmol.rhs().subs(vdict)
    vdict[E_l] = eq_El.rhs().subs(vdict)
    vdict[R_ll] = eq_Rll.rhs().subs(vdict)
    vdict[H_l] = eq_Hl.rhs().subs(vdict)

    
    # Test for steady state
    if n((E_l + H_l + R_ll - R_s).subs(vdict))>1.:
        return 'error in energy balance: El + Hl + Rll - R_s = ' + str(n((E_l + H_l + R_ll - R_s).subs(vdict))) 
    return vdict


# In[6]:

# Test
vdict = cdict.copy()
vdict[a_s] = 1.0    # one sided stomata
vdict[g_sw] = 0.01    
vdict[T_a] = 273 + 25.5
vdict[T_w] = vdict[T_a] # Wall temperature equal to air temperature
vdict[P_a] = 101325
rha = 0.5
vdict[P_wa] = rha*eq_Pwl.rhs()(T_l = T_a).subs(vdict)
vdict[L_l] = 0.03
#vdict[L_A] = vdict[L_l]^2
vdict[Re_c] = 3000
vdict[R_s] = 0.
#vdict[Q_in] = 0
vdict[v_w] = 1

dict_print(fun_SS(vdict))


# ## Gas and energy exchange in a leaf chamber
# Calculations based on leaf_capacitance_steady_state1. However, following the LI-6400XT user manual (Eq. 17-3), we replace the air temperature by wall temperature in the calculation of the net longwave balance of the leaf, as wall temperature can be measured in the chamber. Following the same equation, we also add the leaf thermal emissivity of 0.95 (P. 17-3). 
# **Note that in order to measure sensible heat flux from the leaf, wall temperature must be equal to air temperature!**

# In[7]:

width = 0.05
height = 0.03
volume = 310/(100^3)
print 'Volume = ' + str((volume*1000).n()) + ' l'

print 'min flow rate for flushing = ' + str((volume*3/100^3/60).n()) + ' m3/s'
print 'min flow rate for flushing = ' + str((volume*3*1000).n()) + ' l/min'
print 'flow rate for 1 m/s direct wind = ' + str(width*height*1) + ' m3/s'
print 'flow rate for 1 m/s direct wind = ' + str(width*height*1*1000*60) + ' l/m'
print 'flow rate for 5 m/s direct wind = ' + str(width*height*5) + ' m3/s'
print 'flow rate for 5 m/s direct wind = ' + str(width*height*5*1000*60) + ' l/m'


# In[8]:

width = 0.05
height = 0.03
volume = 400/(100^3)
print 'Volume = ' + str((volume*1000).n()) + ' l'

print 'min flow rate for flushing = ' + str((volume*3/100^3/60).n()) + ' m3/s'
print 'min flow rate for flushing = ' + str((volume*3*1000).n()) + ' l/min'
print 'flow rate for 1 m/s direct wind = ' + str(width*height*1) + ' m3/s'
print 'flow rate for 1 m/s direct wind = ' + str(width*height*1*1000*60) + ' l/m'


# <h2>Chamber insulation material</h2>

# In[9]:

var2('c_pi', 'Heat capacity of insulation material', joule/kilogram/kelvin, latexname='c_{pi}')
var2('lambda_i', 'Heat conductivity of insulation material', joule/second/meter/kelvin)
var2('rho_i', 'Density of insulation material', kilogram/meter^3)
var2('L_i', 'Thickness of insulation material', meter)
var2('A_i', 'Conducting area of insulation material', meter^2)
var2('Q_i', 'Heat conduction through insulation material', joule/second)
var2('dT_i', 'Temperature increment of insulation material', kelvin)


# In[10]:

assumptions(c_pi)


# In[11]:

eq_Qi = Q_i == dT_i*lambda_i*A_i/L_i
units_check(eq_Qi)


# In[12]:

eq_Li = solve(eq_Qi, L_i)[0]
units_check(eq_Li)


# In[13]:

# The amount of heat absorbed by the insulation material per unit area to increase the wall temperature by the same amount as dT_i for given heat flux Q_i
units_check(c_pi*rho_i*dT_i*L_i)


# In[14]:

(c_pi*rho_i*dT_i*L_i).subs(eq_Li)


# In[15]:

# From http://www.waermedaemmstoffe.com/
vdict_PS = {}
vdict_PS[lambda_i] = 0.0375
vdict_PS[c_pi] = 1500
vdict_PS[rho_i] = 22.5

vdict_PU = {}
vdict_PU[lambda_i] = 0.025
vdict_PU[c_pi] = 1300
vdict_PU[rho_i] = 32.5


# In[16]:

(c_pi*rho_i*dT_i*L_i).subs(eq_Li).subs(vdict_PS)(A_i = 0.3, dT_i = 0.1, Q_i = 0.01)


# In[17]:

(c_pi*rho_i*dT_i*L_i).subs(eq_Li).subs(vdict_PU)(A_i = 0.3, dT_i = 0.1, Q_i = 0.01)


# In[18]:

# Best PS
vdict_PSb = {}
vdict_PSb[lambda_i] = 0.035
vdict_PSb[c_pi] = 1500
vdict_PSb[rho_i] = 10
(c_pi*rho_i*dT_i*L_i).subs(eq_Li).subs(vdict_PSb)(A_i = 0.3, dT_i = 0.1, Q_i = 0.01)


# In[19]:

# Best PU
vdict_PUb = {}
vdict_PUb[lambda_i] = 0.02
vdict_PUb[c_pi] = 1200
vdict_PUb[rho_i] = 30
(c_pi*rho_i*dT_i*L_i).subs(eq_Li).subs(vdict_PUb)(A_i = 0.3, dT_i = 0.1, Q_i = 0.01)


# In[20]:

# From http://www.sager.ch/default.aspx?navid=25, PIR
vdict = {}
vdict[lambda_i] = 0.022
vdict[c_pi] = 1400
vdict[rho_i] = 30
(c_pi*rho_i*dT_i*L_i).subs(eq_Li).subs(vdict)(A_i = 0.3, dT_i = 0.1, Q_i = 0.01)


# In[21]:

# From http://www.sager.ch/default.aspx?navid=25, PIR
vdict = {}
vdict[lambda_i] = 0.022
vdict[c_pi] = 1400
vdict[rho_i] = 30
(c_pi*rho_i*dT_i*L_i).subs(eq_Li).subs(vdict)(A_i = 0.3, dT_i = 0.1, Q_i = 0.01)


# In[22]:

# From http://www.sager.ch/default.aspx?navid=25, Sagex 15
vdict = {}
vdict[lambda_i] = 0.038
vdict[c_pi] = 1400
vdict[rho_i] = 15
(c_pi*rho_i*dT_i*L_i).subs(eq_Li).subs(vdict)(A_i = 0.3, dT_i = 0.1, Q_i = 0.01)


# In[23]:

units_check(lambda_i*A_i*dT_i*L_i)


# In[24]:

# Assuming a 30x10x5 cm chamber, how thick would the insulation have to be in order to lose less 0.01 W heat for 0.1 K dT_i?
eq_Li(A_i = 0.3*0.1*2 + 0.3*0.05*2, dT_i = 0.1, Q_i = 0.01).subs(vdict)


# In[25]:

n(10000/399)


# In[26]:

# From http://www.sager.ch/default.aspx?navid=25, Sagex 30
vdict = {}
vdict[lambda_i] = 0.033
vdict[c_pi] = 1400
vdict[rho_i] = 30
(c_pi*rho_i*dT_i*L_i).subs(eq_Li).subs(vdict)(A_i = 0.3, dT_i = 0.1, Q_i = 0.01)


# <h1>Definition of variables</h1>

# In[27]:

# Additional variables
var2('alpha_l', 'Leaf albedo, fraction of shortwave radiation reflected by the leaf', watt/watt)
var2('R_d', 'Downwelling global radiation', joule/second/meter^2)
var2('R_la', 'Longwave radiation absorbed by leaf', joule/second/meter^2, latexname='R_{la}')
var2('R_ld', 'Downwards emitted/reflected global radiation from leaf', joule/second/meter^2, latexname='R_{ld}')
var2('R_lu', 'Upwards emitted/reflected global radiation from leaf', joule/second/meter^2, latexname='R_{lu}')
var2('R_u', 'Upwelling global radiation', joule/second/meter^2)
var2('S_a', 'Radiation sensor above leaf reading', joule/second/meter^2)
var2('S_b', 'Radiation sensor below leaf reading', joule/second/meter^2)
var2('S_s', 'Radiation sensor beside leaf reading', joule/second/meter^2)


# In[28]:

variables = sorted([str(variable) for variable in udict.keys()],key=str.lower)
tabledata = [('Variable', 'Description', 'Units')] + [(eval(variable),docdict[eval(variable)],udict[eval(variable)]/eval(variable)) for variable in variables]
table(tabledata)


# In[29]:

#alpha_a = var2('alpha_a', 'Thermal diffusivity of dry air (m2/s)',meter^2/second)
#a_s = var2('a_s', 'Fraction of projected leaf area covered by stomata (1 if stomata are on one side only, 2 if they are on both sides)', meter^2/meter^2)
#a_sh = var2('a_sh', 'Fraction of projected leaf area exchanging sensible heat flux (1 if one side only, 2 if both sides)', meter^2/meter^2, value = 2)
#alpha_l = var2('alpha_l', 'Leaf albedo, fraction of shortwave radiation reflected by the leaf', watt/watt)
#c_pa = var2('c_pa', 'Specific heat of dry air (1010 J kg-1 K-1)', joule/kilogram/kelvin, latexname = 'c_{pa}', value = 1005) # http://www.engineeringtoolbox.com/air-properties-d_156.html
#c_pv = var2('c_pv', 'Specific heat of water vapour at 300 K', joule/kilogram/kelvin, latexname = 'c_{pa}', value = 1864) # source: http://www.engineeringtoolbox.com/water-vapor-d_979.html
#c_pl = var2('c_pl', 'The leaf heat capacity per leaf area (J K-1 m-2)',joule/kelvin/meter^2,'c_{pl}')
#c_pw = var2('c_pw', 'Heat capacity of liquid water at constant pressure',joule/kilogram/kelvin,'c_{pw}', value = 4187)         # J/kgK, source: http://www.engineeringtoolbox.com/water-thermal-properties-d_162.html
#C_va = var2('C_va', 'Concentration of water in the free air (mol m-3)',mole/meter^3, latexname = 'C_{va}')
#C_vl = var2('C_vl', 'Concentration of water in the leaf air space (mol m-3)',mole/meter^3, latexname = 'C_{vl}')
#D_va = var2('D_va', 'Binary diffusion coefficient of water vapour in air (m2/s)',meter^2/second, latexname = 'D_{va}')
#eps_l = var2('eps_l', 'Thermal emissivity of the leaf', value = 0.95, latexname = '\\epsilon_l')
#E_lmol = var2('E_lmol', 'Transpiration rate in molar units',mole/second/meter^2,latexname='E_{l,mol}')
#E_l = var2('E_l', 'Latent heat flux from leaf (W m-2)',joule/second/meter^2)
#E_l = var2('E_l', 'Latent heat flux from leaf (W m-2)',joule/second/meter^2)
#g = var2('g', 'Gravitational acceleration', meter/second^2, value= 9.81)
#g_bv = var2('g_bv', 'Boundary layer conductance to water vapour (m/s)',meter/second,latexname='g_{bv}')
#Gr = var2('Gr', 'Grashof number', latexname = 'N_{Gr_L}')
#g_sv = var2('g_sv', 'Stomatal conductance to water vapour (m/s)',meter/second,latexname='g_{sv}')
#g_tv = var2('g_tv', 'Total leaf conductance to water vapour (m/s)',meter/second,latexname='g_{tv}')
#g_vmol = var2('g_vmol', 'Total leaf conductance to water vapour (mol/m2/s)',mole/meter^2/second,latexname='g_{v,mol}')
#h_cl = var2('h_cl', 'free convection coefficient for lower surface', joule/(kelvin*meter^2*second))
#h_cu = var2('h_cu', 'free convection coefficient for upper surface', joule/(kelvin*meter^2*second))
#h_c = var2('h_c', 'average convective transfer coefficient', joule/(kelvin*meter^2*second))
#H_ll = var2('H_ll','sensible heat flux from lower side of a plate',joule/second/meter^2)
#H_lu = var2('H_lu','sensible heat flux from upper side of a plate',joule/second/meter^2)
#H_l = var2('H_l', 'Sensible heat flux from leaf (W m-2)',joule/second/meter^2)
#k_a = var2('k_a', 'Thermal conductivity of dry air (W m-1 K-1)',joule/second/meter/kelvin)
#k_l = var2('k_l', 'Thermal conductivity of a fresh leaf', joule/second/meter/kelvin)   # between 0.268 and 0.573 according to Hays_1975_The_thermal.pdf
#lambda_E = var2('lambda_E', 'Latent heat of evaporation (J kg-1)',joule/kilogram,value = 2.45e6)
#Le = var2('Le', 'Lewis number, Le = kappa_a/D_va', latexname = 'N_{Le}')
#L_l = var2('L_l', 'characteristic length scale for convection (size of leaf)',meter)
#M_N2 = var2('M_N2', 'Molar mass of nitrogen (kg mol-1)',kilogram/mole,value = 0.028)
#M_O2 = var2('M_O2', 'Molar mass of oxygen (kg mol-1)',kilogram/mole,value = 0.032)
#M_air = var2('M_air', 'Molar mass of air (kg mol-1)', kilogram/mole, value = 0.02897)  # http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
#M_w = var2('M_w', 'Molar mass of water (kg mol-1)',kilogram/mole,value = 0.01802)
#nu_a = var2('nu_a', 'kinematic viscosity of dry air (m2/s)',meter^2/second)
#Nu = var2('Nu', 'Nusselt number',latexname = 'N_{Nu_L}')
#phi_a = var2('phi_a', 'Relative humidity in the leaf', value = 1)
#phi_l = var2('phi_l', 'Relative humidity in the leaf', value = 1)
#P_a = var2('P_a', 'Air pressure',pascal)
#P_r = var2('P_r', 'Reference pressure',pascal, value = 101325)
#Pr = var2('Pr', 'Prandtl number, 0.71 for air',latexname = 'N_{Pr}',value = 0.71)
#P_N2 = var2('P_N2', 'Partial pressure of nitrogen in the atmosphere', pascal, latexname = 'P_{N2}')
#P_N2 = var2('P_O2', 'Partial pressure of oxygen in the atmosphere', pascal, latexname = 'P_{O2}')
#P_va = var2('P_va', 'Vapour pressure in the atmosphere', pascal, latexname = 'P_{va}')
#P_vl = var2('P_vl', 'Vapour pressure inside the leaf', pascal, latexname = 'P_{va}')
#P_v = var2('P_v', 'Vapour pressure', pascal, pascal)
#Q_in = var2('Q_in', 'Internal heat sources, such as fan', joule/second, latexname = 'Q_{in}')
#Q_l = var2('Q_l', 'Conductive heat flux from upper to lower side of leaf', joule/second/meter^2)
#Re_c = var2('Re_c', 'Critical Reynolds number for the onset of turbulence',latexname='N_{Re_c}')
#Re = var2('Re', 'Reynolds number',latexname='N_{Re_L}')
#rho_a = var2('rho_a', 'density of dry air', kilogram/meter^3)
#rho_al = var2('rho_al', 'Density of air at the leaf surface', kilogram/meter^3)
#rho_w = var2('rho_w', 'density of water (kg m-3)',kilogram/meter^3,value = 1000)
#R_d = var2('R_d', 'Downwelling global radiation', joule/second/meter^2)
#R_la = var2('R_la', 'Longwave radiation absorbed by leaf', joule/second/meter^2)
#R_ld = var2('R_ld', 'Downwards emitted/reflected global radiation from leaf', joule/second/meter^2)
#R_lu = var2('R_lu', 'Upwards emitted/reflected global radiation from leaf', joule/second/meter^2)
#R_u = var2('R_u', 'Upwelling global radiation', joule/second/meter^2)
#R_ll = var2('R_ll', 'Net longwave radiation away from leaf (W m-2)',joule/second/meter^2,'R_{ll}')
#R_mol = var2('R_mol', 'Molar gas constant (J/mol/K)',joule/mole/kelvin,latexname='R_{mol}', value = 8.314472)
#R_s = var2('R_s', 'Solar shortwave flux absorbed by leaf',joule/second/meter^2)
#S_a = var2('S_a', 'Radiation sensor above leaf reading', joule/second/meter^2)
#S_b = var2('S_b', 'Radiation sensor below leaf reading', joule/second/meter^2)
#S_s = var2('S_s', 'Radiation sensor beside leaf reading', joule/second/meter^2)
#Sh = var2('Sh', 'Sherwood number',latexname = 'N_{Sh_L}')
#sigm = var2('sigm', 'Stefan-Boltzmann constant (5.67e-8 Wm-2 K-4)',joule/second/meter^2/kelvin^4,'real','\sigma',5.67e-8)
#t = var2('t', 'Time', second, 'real')
#T0 = var2('T0', 'Freezing point in kelvin',kelvin, value = 273.15)
#T_a = var2('T_a', 'Air temperature (K)',kelvin)
#T_b = var2('T_b', 'Boundary layer temperature (K)',kelvin)
#T_l = var2('T_l', 'Leaf temperature (K)',kelvin)
#T_lu = var2('T_lu', 'Leaf surface temperature of upper side', kelvin)
#T_ll = var2('T_ll', 'Leaf surface temperature of lower side', kelvin)
#T_r = var2('T_r', 'Reference temperature', kelvin)
#T_w = var2('T_w', 'Chamber wall temperature (K)', kelvin)
#v_w = var2('v_w', 'wind velocity (m/s)',meter/second)
#z_l = var2('z_l', 'Leaf thickness (m)',meter)


# 
# 
# ![GitHub Logo](/images/Leaf_radbalance.png)
# Format: ![Alt Text](url)
# 
# 

# <h2>Leaf radiation balance</h2>
# 
# <p>The leaf is exposed to downwelling radiation ($R_d$) originating from shortwave irradiance entering through the glass window plus the longwave irradiance transmitted througha and emitted by the glass window, plus the upwelling radiation ($R_u$) emitted by the bottom glass window.</p>
# <p>The leaf itself reflects some of the radiation in both direction and emits its own black body longwave radiation. The sum of refelcted and emitted radiation away from the leaf is denoted as $R_{lu}$ and $R_{ld}$ for upward and downwards respectively. We have three net radiation sensors in place, one above the leaf measuring $S_a$, one below the leaf measureing $S_b$ and one at the same level beside the leaf measureing $S_s$. These sensor measure:</p>
# <p><img style="float: right;" src="http://localhost:8888/notebooks/Schymanski_Or_experimental_2016/images/Leaf_radbalance.png" alt="" width="400" height="300" /></p>
# <p>$S_a = R_d - R_{lu}$</p>
# <p>$S_b = R_{ld} - R_u$</p>
# <p>$S_s = R_d - R_u$</p>
# <p>This leaves us with 3 equations with 4 unknows, so we either have to assume that $R_{lu} = R_{ld}$, assuming that both sides of the leaf have the same temperature or $R_u = 0$ to solve the algebraic problem. In daylight, $R_d >> R_u$, so this assumption should not lead to a big bias, however this would imply that $S_b = R_{ld}$, which is certainly incorrect.</p>
# <p>Unfortunately, the assumption $R_{lu} = R_{ld}$ does not help solve the problem as it just implies that $S_s = S_a + S_b$:</p>

# In[30]:

eq_Rs_Rd = R_s == (1-alpha_l)*R_d
eq_Sa = S_a == R_d - R_lu
eq_Sb = S_b == R_ld - R_u
eq_Ss = S_s == R_d - R_u


# In[31]:

# Assuming R_lu = R_ld
eq_assRldRlu = R_ld == R_lu
solve([eq_Sa, eq_Sb.subs(eq_assRldRlu), eq_Ss], R_d, R_lu, R_u)


# In[32]:

# More specifically,
eq1 = solve(eq_Sa, R_d)[0]
eq2 = solve(eq_Sb, R_ld)[0].subs(eq_assRldRlu)
eq3 = solve(eq_Ss, R_u)[0] 
solve(eq1.subs(eq2).subs(eq3), S_s)


# In[33]:

# Assuming that R_u = 0
eq_assRu0 = R_u == 0
soln = solve([eq_Sa, eq_Sb, eq_Ss, eq_assRu0], R_d, R_lu, R_ld, R_u)
print soln


# In[34]:

#eq_Rd = soln[0][0]
#eq_Rlu = soln[0][1]
#eq_Rld = soln[0][2]
#eq_Rlu


# <p>However, what we can do in any case is to quantify the net radiative energy absorbed by the leaf as</p>
# <p>$\alpha_l R_s - R_{ll} = S_a - S_b$:</p>

# In[35]:

# Leaf radiation balance
eq_Rs_Rll = R_s - R_ll == R_d + R_u - R_lu - R_ld
eq_Rbalance = R_s - R_ll == S_a - S_b


# In[36]:

solve([eq_Sa, eq_Sb, eq_Ss, R_d + R_u - R_lu - R_ld == S_a - S_b], R_d, R_lu, R_ld, R_u)


# <h2>Leaf water vapour exchange and energy balace</h2>

# In[37]:

## From leaf_capacitance_steady_state1
#dTl_dt = (R_s - R_ll - H_l - E_l)/c_pl
#eq_cpl = c_pl == z_l*rho_w*c_pw
#eq_Rll = R_ll == a_sh*sigm*(T_l^4 - T_w^4)
#eq_hc = h_c == k_a*Nu/L_l
#eq_H_lu = H_lu == -(T_a - T_l)*h_cu
#eq_H_ll = H_ll == -(T_a - T_l)*h_cl
#eq_H_l = H_l == a_sh*(T_l - T_a)*h_c
#eq_Re = Re == L_l*v_w/nu_a
#eq_Tb = T_b == (T_a + T_l)/2
#eq_Nu_forced_lam = Nu == 83/125*Pr^(1/3)*sqrt(Re)
#eq_Nu_forced_mix = Nu == 1/1000*(37*Re^(4/5) - 37*Re_c^(4/5) + 664*sqrt(Re_c))*Pr^(1/3)
#eq_Nu_forced_all = eq_Nu_forced_mix(Re_c = (Re_c - (abs(Re - Re_c) - (Re - Re_c))/2))
#eq_Elmol = E_lmol == g_tv*(C_vl - C_va)
#eq_gtv = g_tv == g_sv*g_bv/(g_sv + g_bv)
#eq_Cva = C_va == P_va/(R_mol*T_a)
#eq_Pvl = P_vl == 611*exp(lambda_E*M_w/R_mol*(1/273 - 1/T_l))
#eq_El = E_l == eq_Elmol.rhs()*M_w*lambda_E
#eq_gbv = g_bv == Sh*D_va/L_l
#eq_gbv_hc = g_bv == a_s*h_c/(rho_a*c_pa*Le^(1-1/3))
#eq_Le = Le == alpha_a/D_va
#eq_rho_a = rho_a == ((M_w*P_va + M_N2*0.79*(P_a - P_va) + M_O2*0.21*(P_a - P_va))/(R_mol*T_a))


# In[38]:

## From Table A.3 in Monteith&Unsworth
#temp_tab = np.array([268.2, 273.2, 278.2, 283.2, 288.2, 293.2, 298.2, 303.2, 308.2, 313.2, 318.2])
#Dva_tab = np.array([20.5,21.2,22,22.7,23.4,24.2,24.9,25.7,26.4,27.2,28])*1e-6
#var('x m c') 
#model(x) = m*x + c
#eq_D_va = D_va == model(T_a).subs(find_fit(zip(temp_tab,Dva_tab), model,solution_dict=True))
## From Leaf_energy_balance_dry
#eq_alpha_a = alpha_a == (1.3236363633641127e-07)*T_a - 1.7281745441936981e-05
#eq_k_a = k_a == (6.836363632085321e-05)*T_a + 0.005628509103451979
#eq_nu_a = nu_a == (8.9999999844637335e-08)*T_a - 1.1297090863409663e-05


# In[39]:

#eq_hc1 = eq_hc.subs(eq_Nu_forced_all).subs(eq_Re).subs(eq_k_a,eq_nu_a).subs(eq_Tb)
#eq_El_conv = E_lmol == (P_vl - P_va)/P_a*g_vmol
#eq_Elmol1 = eq_Elmol.subs(eq_gtv,eq_Cva,C_vl == eq_Cva.rhs()(T_a = T_l, P_va = P_vl)).subs(eq_Pvl).subs(eq_gbv_hc).subs(eq_Le,eq_rho_a).subs(eq_alpha_a,eq_D_va).subs(eq_Tb)
#eq_El1 = eq_El.subs(eq_gtv,eq_Cva,C_vl == eq_Cva.rhs()(T_a = T_l, P_va = P_vl)).subs(eq_Pvl).subs(eq_gbv_hc).subs(eq_Le,eq_rho_a).subs(eq_alpha_a,eq_D_va).subs(eq_Tb)
#eq_El2 = eq_El.subs(eq_Cva,C_vl == eq_Cva.rhs()(T_a = T_l, P_va = P_vl)).subs(eq_Pvl)
#dTl_dt1 = dTl_dt.subs(H_l == eq_H_l.rhs(), eq_El1, eq_Rll, eq_cpl).subs(eq_hc1).subs(cdict)


# In[40]:

## Exact conversions between molar and length units for conductances
#soln = solve(eq_Elmol.rhs().subs(eq_Cva, C_vl == eq_Cva.rhs()(P_va = P_vl, T_a = T_l)) == eq_El_conv.rhs(), g_tv)
#eq_gtv_gvmol = soln[0]
#soln = solve(eq_Elmol.rhs().subs(eq_Cva, C_vl == eq_Cva.rhs()(P_va = P_vl, T_a = T_l)) == eq_El_conv.rhs(), g_vmol)
#eq_gvmol_gtv = soln[0]


# In[41]:

eq_rhoa_Pwa_Ta.subs(eq_PN2, eq_PO2).subs(cdict)


# In[42]:

# Comparison with http://www.engineeringtoolbox.com/density-air-d_680.html
x = 0.62198*P_wa / (P_a - P_wa) 
rhoda = (P_a) / (286.9*T_a)
rhoa = rhoda*(1 + x) / (1 + 1.609*x )
rhoa


# In[43]:

(R_mol/M_w).subs(cdict)


# In[44]:

1/286.9


# In[45]:

vdict = cdict.copy()
vdict[P_a] = 101325
vdict[P_wa] = 100
P = plot(eq_rhoa_Pwa_Ta.rhs().subs(eq_PN2, eq_PO2).subs(vdict), (T_a, 280,310), frame = True, axes = False, legend_label = 'eq_rho_a')
P += plot(rhoa.subs(vdict), (T_a, 280,310), color = 'green', legend_label = 'Engineeringtoolbox')
P.axes_labels(['$T_a$ (K)', '$\\rho_a$ (kg m$^{-3}$)'])
P


# <p>The difference in the above is due to 0.01802 kg/mol as molecular weight of water and 0.02897 kg/mol as molecular weight of air in the Engieeringtoolbox (<a href="http://www.engineeringtoolbox.com/humid-air-ideal-gas-d_677.html">http://www.engineeringtoolbox.com/humid-air-ideal-gas-d_677.html</a>)</p>

# In[46]:

# Grashof number to calculate transition between forced and free convection (note that for some Ra-correlations, L_l is the leaf area divided by perimeter
#eq_Gr = Gr == g*(rho_a - rho_al)/rho_al*L_l^3/nu_a^2


# In[47]:

vdict = cdict.copy()
vdict[L_l] = 0.05
vdict[T_a] = 300
#vdict[T_l] = 310

vdict[P_a] = 101325.
vdict[P_wa] = eq_Pwl.rhs()(T_l = T_a).subs(vdict)
vdict[P_wl] = eq_Pwl.rhs().subs(vdict)

vdict


# In[48]:

vdict[L_l]^2/(vdict[L_l]*4)


# In[49]:

n(5^2/(5*4))


# In[50]:

vdict1 = vdict.copy()
vdict1[rho_al] = eq_rhoa_Pwa_Ta.rhs().subs(eq_PN2, eq_PO2)(P_wa = P_wl, T_a = T_l).subs(vdict1)
print vdict1[rho_al](T_l = 310)
vdict1[rho_a] = eq_rhoa_Pwa_Ta.rhs().subs(eq_PN2, eq_PO2).subs(vdict1)
print vdict1[rho_a](T_l = 310)
vdict1[nu_a] = eq_nua.rhs().subs(vdict1)
print vdict1[nu_a](T_l = 310)
eq_Gr.subs(vdict1)(T_l = 310)


# In[51]:

eq_Gr.subs(eq_nua, eq_rhoa_Pwa_Ta, rho_al == eq_rhoa_Pwa_Ta.rhs()(P_wa = P_wl, T_a = T_l)).subs(vdict)(T_l = 310)


# In[52]:

eq_Gr.subs(L_l == vdict[L_l]^2/(vdict[L_l]*4),eq_nua, eq_rhoa_Pwa_Ta, rho_al == eq_rhoa_Pwa_Ta.rhs()(P_wa = P_wl, T_a = T_l)).subs(eq_PN2, eq_PO2).subs(vdict)(T_l = 310)


# In[53]:

P = plot(eq_Gr.rhs().subs(eq_nua, eq_rhoa, rho_al == eq_rhoa.rhs()(P_wa = P_wl, T_a = T_l)).subs(vdict), (T_l, 300,310), frame = True, axes = False)
P.axes_labels(['Leaf temperature (K)', 'Gr'])
P


# In[54]:

eq_Re.subs(eq_nua).subs(vdict)(T_l = 310)


# In[55]:

P = plot((eq_Re.rhs()^2).subs(eq_nua).subs(vdict)(T_l = 310), (v_w, 0.01,0.5), frame = True, axes = False)
P.axes_labels(['Wind speed (m s$^{-1}$)', 'Re$^2$'])
P.fontsize(16)
P


# In[56]:

# Blasius solution for BL thickness (http://en.wikipedia.org/wiki/Boundary-layer_thickness)
var2('B_l', 'Boundary layer thickness', meter)
vdict = cdict.copy()
Ta = 300
vdict[T_a] = Ta
vdict[T_l] = Ta
vdict[L_l] = 0.15
vdict[v_w] = 0.5
vdict[Re_c] = 3000
vdict[a_s] = 1
vdict[P_a] = 101325
vdict[P_wa] = 0
eq_Bl = B_l == 4.91*sqrt(nu_a*L_l/v_w)
print eq_Bl.subs(eq_nua).subs(vdict)
eq_gbw_hc.subs(eq_hc).subs(eq_Nu_forced_all).subs(eq_Re).subs(eq_rhoa).subs(eq_Le).subs(eq_Dva).subs(eq_alphaa, eq_nua, eq_ka).subs(vdict)


# In[57]:

vdict[L_l] = 0.03
vdict[v_w] = 8
print eq_Bl.subs(eq_nua).subs(vdict)
eq_gbw_hc.subs(eq_hc).subs(eq_Nu_forced_all).subs(eq_Re).subs(eq_rhoa).subs(eq_Le).subs(eq_Dva).subs(eq_alphaa, eq_nua, eq_ka).subs(vdict)


# In[58]:

# Does B_l scale with g_bw?
vdict[v_w] = v_w
plot((eq_gbw_hc/eq_Bl).rhs().subs(eq_hc).subs(eq_Nu_forced_all).subs(eq_Re).subs(eq_rhoa).subs(eq_Le).subs(eq_Dva).subs(eq_alphaa, eq_nua, eq_ka).subs(vdict), (v_w, 0.5,5))


# In[ ]:

# Maximum sensible heat flux of a 3x3 cm leaf irradiated by 600 W/m2
600*0.03^2


# <h1>Chamber mass and energy balance</h1>
# <p>Usually, we know the volumetric inflow into the chamber, so to convert to molar inflow (mol s$^{-1}$), we will use the ideal gas law: $P_a V_c = n R_{mol} T_{in}$, where $n$ is the amount of matter in the chamber (mol). To convert from a volume to a flow rate, we replace $V_c$ by $F_{in,v}$. Note that partial pressures of dry air and vapour are additive, such that</p>
# <p>$P_a = P_w + P_d$</p>
# <p>However, the volumes are not additive, meaning that:</p>
# <p>$P_d V_a = n_d R_{mol} T_{a}$</p>
# <p>$(P_a - P_d) V_a = n_a R_{mol} T_{a}$</p>
# <p>i.e. we use the same volume ($V_a$) for both the vapour and the dry air. This is because both the vapour and the dry air are well mixed and occupy the same total volume. Their different amounts are expressed in their partial pressures. If we wanted to calculate the partial volumes they would take up in isolation from each other, we would need to specify at which pressure this volume is taken up and if we used the same pressure for both (e.g. $P_a$), we would obtain a volume fraction for water vapour equivalent to its partial pressure fraction in the former case.</p>
# <p>Therefore, we will distinguish the molar flow rates of water vapour ($F_{in,mol,v}$) and dry air ($F_{in,mol,a}$) but they share a common volumetric flow rate ($F_{in,v}$).</p>

# In[ ]:

var2('W_c', 'Chamber width', meter)
var2('L_c', 'Chamber length', meter)
var2('H_c', 'Chamber height', meter)
var2('V_c', 'Chamber volume', meter^3)
var2('n_c', 'molar mass of gas in chamber', mole)
var2('F_in_v', 'Volumetric flow rate into chamber', meter^3/second, latexname='F_{in,v}')
#var2('F_in_va', 'Volumetric flow rate of dry air into chamber', meter^3/second, latexname='F_{in,v,a}')
#var2('F_in_vw', 'Volumetric flow rate of water vapour into chamber', meter^3/second, latexname='F_{in,v,w}')
var2('F_in_mola', 'Molar flow rate of dry air into chamber', mole/second, latexname='F_{in,mol,a}')
var2('F_in_molw', 'Molar flow rate of water vapour into chamber', mole/second, latexname='F_{in,mol,w}')
var2('F_out_mola', 'Molar flow rate of dry air out of chamber', mole/second, latexname='F_{out,mol,a}')
var2('F_out_molw', 'Molar flow rate of water vapour out of chamber', mole/second, latexname='F_{out,mol,w}')
var2('F_out_v', 'Volumetric flow rate out of chamber', meter^3/second, latexname='F_{out,v}')
#var2('F_out_va', 'Volumetric flow rate of dry air out of chamber', meter^3/second, latexname='F_{out,v,a}')
#var2('F_out_vw', 'Volumetric flow rate of water vapour out of chamber', meter^3/second, latexname='F_out,v,w}')
var2('T_d', 'Dew point temperature of incoming air', kelvin)
var2('T_in', 'Temperature of incoming air', kelvin, latexname='T_{in}')
var2('T_out', 'Temperature of outgoing air (= chamber T_a)', kelvin, latexname='T_{out}')
var2('T_room', 'Lab air temperature', kelvin, latexname='T_{room}')
var2('P_w_in', 'Vapour pressure of incoming air', pascal, latexname='P_{w,in}')
var2('P_w_out', 'Vapour pressure of outgoing air', pascal, latexname='P_{w,out}')
var2('R_H_in', 'Relative humidity of incoming air', latexname='R_{H,in}')
var2('L_A', 'Leaf area', meter^2)


# In[ ]:

eq_V_c = fun_eq(V_c == W_c*L_c*H_c)
eq_F_in_mola = fun_eq(F_in_mola == (P_a - P_w_in)*F_in_v/(R_mol*T_in))
eq_F_in_molw = fun_eq(F_in_molw == (P_w_in)*F_in_v/(R_mol*T_in))
eq_F_out_mola = fun_eq(F_out_mola == (P_a - P_w_out)*F_out_v/(R_mol*T_out))
eq_F_out_molw = fun_eq(F_out_molw == (P_w_out)*F_out_v/(R_mol*T_out))
eq_F_out_v = fun_eq(F_out_v == (F_out_mola + F_out_molw)*R_mol*T_out/P_a)


# <p>At steady state, $F_{out,mola} = F_{in,mola}$ and $F_{out,molw} = F_{in,molw} + E_{l,mol} L_A$. In the presence of evaporation, we can simply add Elmol to get F_out_v as a function of F_in_v</p>
# <p>Assuming that the pressure inside the chamber is constant and equal to the pressure outside, we compute the change in volumetric outflow due to a change in temperature and due to the input of water vapour by transpiration as:</p>

# In[ ]:

eq_Foutv_Finv_Tout = eq_F_out_v.subs(F_out_mola = F_in_mola, F_out_molw = F_in_molw + E_lmol*L_A).subs(eq_F_in_mola, eq_F_in_molw).simplify_full()
units_check(eq_Foutv_Finv_Tout)


# In[ ]:

eq_Foutv_Finv_Tout


# In[ ]:

# Other way, using molar in and outflow
eq_Foutmolw_Finmolw_Elmol = F_out_molw == (F_in_molw + E_lmol*L_A)
units_check(eq_Foutmolw_Finmolw_Elmol)


# In[ ]:




# <h2>Change in air temperature</h2>
# <p>See also <a href="http://www.engineeringtoolbox.com/mixing-humid-air-d_694.html">http://www.engineeringtoolbox.com/mixing-humid-air-d_694.html</a> and <a href="http://www.engineeringtoolbox.com/enthalpy-moist-air-d_683.html">http://www.engineeringtoolbox.com/enthalpy-moist-air-d_683.html</a> for reference.</p>
# <p>We will assume that the air entering the chamber mixes with the air inside the chamber at constant pressure, i.e. the volume of the mixed air becomes the chamber volume plus the volume of the air that entered. The temperature of the mixed air is then the sum of their enthalpies plus the heat added by the fan and by sensible heaflux, divided by the sum of their heat capacities. The addition of water vapour through evaporation by itself should not affect the air temperature, but the volume of the air.</p>
# <p> </p>
# <p>Alternatively, we could assume that a given amount of air is added to a constant volume, leading to an increase in pressure. Addition of water vapour would lead to an additional increase in pressure. In addition, addition/removal of heat by sensible heat flux and the chamber fan would affect both temperature and pressure.To calculate both temperature and pressure, we need to track the internal energy in addition to the mole number. According to Eq. 6.1.3 in Kondepudi & Prigogine (2006), the internal energy of an ideal gas is given by (see also Eq. 2.2.15 in Kondepuid & Prigogine):</p>
# <p>$U = N(U_0 + C_v T)$</p>
# <p>where</p>
# <p>$U_0 = M c^2$</p>
# <p>The relation between molar heat capacities at constant pressure and volume is given as :</p>
# <p>$C_v = C_p - R_{mol}$</p>
# <p>Any heat exchanged by sensible heat flux, across the walls and the fan can be added to total $U$, and then knowledge about total $C_v$ will let us calculate air temperature inside the chamber. After that, we can use the ideal gas law to calculate volume or pressure, depending in which of those we fixed:</p>
# <p>$P V = n R T$</p>
# <p> </p>
# <p>The difference in water vapour pressure and temperature between the incoming and outgoing air is a function of the latent and sensible heat flux, as well as the flow rate. The heat fluxes associated with the incoming and outgoing air are $T_{in} (c_{pa} F_{in,mola} M_{air} + c_{pv} F_{in,molw} M_{w})$ and $T_{out} (c_{pa} F_{out,mola} M_{air} + c_{pv} F_{out,molw} M_{w})$ respectively. The difference between the two plus any additional heat sources/sinks ($Q_{in}$) equals the sensible heat flux at constant air temperature (steady state).</p>

# In[ ]:

units_check(eq_F_out_v)


# In[ ]:

var2('M_air', 'Molar mass of air (kg mol-1)', kilogram/mole, value = 0.02897, latexname='M_{air}')  # http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
var2('c_pv', 'Specific heat of water vapour at 300 K', joule/kilogram/kelvin, latexname = 'c_{pv}', value = 1864) # source: http://www.engineeringtoolbox.com/water-vapor-d_979.html
var2('Q_in', 'Internal heat sources, such as fan', joule/second, latexname = 'Q_{in}')
eq_chamber_energy_balance = 0 == H_l*L_A + Q_in + T_in*(c_pa*M_air*F_in_mola + c_pv*M_w*F_in_molw) - (T_out*(c_pa*M_air*F_out_mola + c_pv*M_w*F_out_molw)); show(eq_chamber_energy_balance)
print latex(eq_chamber_energy_balance)
eq_Hl_enbal = solve(eq_chamber_energy_balance.subs(F_out_mola == F_in_mola, F_out_molw == F_in_molw + L_A*E_lmol), H_l)[0].expand()
print units_check(eq_Hl_enbal)
latex(eq_Hl_enbal)


# In[ ]:

soln = solve(eq_chamber_energy_balance.subs(eq_Foutmolw_Finmolw_Elmol, F_out_mola == F_in_mola), T_out)
print soln
eq_Tout_Finmol_Tin = soln[0]
units_check(eq_Tout_Finmol_Tin).simplify_full().convert().simplify_full()


# In[ ]:

eq_Tout_Finv_Tin = eq_Tout_Finmol_Tin.subs(eq_F_in_mola, eq_F_in_molw).simplify_full()
print eq_Tout_Finv_Tin
show(eq_Tout_Finv_Tin)
units_check(eq_Tout_Finv_Tin).simplify_full().convert()


# In[ ]:




# <p>The molar outflux of dry air equals the molar influx of dry air, while the molar outflux of water vapour equals the molar influx plus the evaporation rate. The sum of both can be used to obtain the volumetric outflow:</p>

# In[ ]:

# F_out_v as function of F_inv and T_in
eq1 = (eq_F_out_molw + eq_F_out_mola).simplify_full(); show(eq1)
eq2 = eq1.subs(F_out_mola == F_in_mola, eq_Foutmolw_Finmolw_Elmol).subs(eq_F_in_mola, eq_F_in_molw); show(eq2)
soln = solve(eq2,F_out_v); print soln
eq_Foutv_Finv = soln[0]; show(eq_Foutv_Finv)
units_check(eq_Foutv_Finv)


# In[ ]:

eq_Foutv_Finv


# In[ ]:

eq_Foutv_Finv_Tout


# In[ ]:

eq_Tout_Finmol_Tin


# In[ ]:

# Finding the T_in that would balance sensible heat release by the plate for given F_inv
assume(F_in_v > 0)
assume(F_out_v > 0)
assume(E_lmol >=0)
show(eq_chamber_energy_balance)
soln = solve(eq_chamber_energy_balance.subs(F_out_mola == F_in_mola).subs(eq_Foutmolw_Finmolw_Elmol).subs(eq_F_in_mola, eq_F_in_molw), T_in)
print soln
eq_T_in_ss = soln[0]
show(eq_T_in_ss)


# In[ ]:

# Calculating Q_in from T_in, F_in_v and T_out
soln = solve(eq_T_in_ss,Q_in)
print soln
eq_Qin_Tin_Tout = soln[0]


# In[ ]:

# Calculating H_l from T_in, F_in_v and T_out
soln = solve(eq_T_in_ss,H_l)
print soln
eq_Hl_Tin_Tout = soln[0]
eq_Hl_Tin_Tout.show()


# In[ ]:

# Calculating H_l from T_in, T_out and Fmol
soln = solve(eq_chamber_energy_balance.subs(F_out_mola == F_in_mola), H_l)
print soln
eq_Hl_Tin_Tout_Fmol = soln[0].simplify_full()
eq_Hl_Tin_Tout_Fmol.show()


# In[ ]:

# Check if eq_Hl_Tin_Tout_Fmol is equivalent to eq_Hl_Tin_Tout
eq_Hl_Tin_Tout_Fmol.subs(eq_Foutmolw_Finmolw_Elmol).subs(eq_F_in_mola, eq_F_in_molw).simplify_full().show()


# <p>The gas pump delivers approximately 6.4 l/min. To achieve the same flow rate from a typical compressed N2 bottle with 10m^3 volume, we could run the experiment for how long?</p>

# In[ ]:

print str(10/(6.4e-3)/60) + ' hours'


# <h2>Calculation of volumetric flow rate based on Cellkraft measurements</h2>
# <p>Cellcraft uses Arden-Buck equation to convert between vapour pressure and dew point (<a href="http://en.wikipedia.org/wiki/Arden_Buck_Equation">http://en.wikipedia.org/wiki/Arden_Buck_Equation</a>).<br />The air flow rate is given by the Cellkraft humidifier in l/min, but it refers to dry air at 0 oC and 101300 Pa (see emails from Joakim Nordlund, 8.7. - 9.7.2014.)</p>
# <p>We will use the reported dew point temperature to obtain the vapour pressure of the air coming out from the Cellkraft humidifier, then the ideal gas law to obtain the molar flow of dry air, leading to three equations with three unknowns:</p>
# <p>$F_{\mathit{in}_{\mathit{mola}}} = \frac{F_{\mathit{in}_{\mathit{va}_{n}}} P_{r}}{{R_{mol}} T_{r}}$</p>
# <p>$F_{\mathit{in}_{v}} = \frac{{\left(F_{\mathit{in}_{\mathit{mola}}} + F_{\mathit{in}_{\mathit{molw}}}\right)} {R_{mol}} T_{\mathit{in}}}{P_{a}}$</p>
# <p>$F_{\mathit{in}_{\mathit{molw}}} = \frac{F_{\mathit{in}_{v}} P_{v_{\mathit{in}}}}{{R_{mol}} T_{\mathit{in}}}$</p>

# In[ ]:

# Calculating vapour pressure in the incoming air from reported dew point by the Cellkraft humidifier
var2('T0', 'Freezing point in kelvin',kelvin, value = 273.15)
eq_Pwin_Tdew = P_w_in == 611.21*exp((18.678 - (T_d - T0)/234.5)*((T_d - T0)/(257.14-T0+T_d))) # c
list_Tdew = np.array(srange(-40,50,6))+273.25
list_Pwin = [eq_Pwin_Tdew.rhs()(T_d = dummy).subs(cdict) for dummy in list_Tdew]
print list_Pwin
P = plot(eq_Pwin_Tdew.rhs().subs(cdict), (T_d, 253,303), frame = True, axes = False, legend_label = 'Arden Buck')
P += plot(eq_Pwl.rhs().subs(cdict), (T_l, 253,303), color = 'red', linestyle = '--', legend_label = 'Clausius-Clapeyron')
P += list_plot(zip(list_Tdew,list_Pwin), legend_label = 'Cellcraft')
P.axes_labels(['Dew point (K)', 'Saturation vapour pressure (Pa)'])
P


# In[ ]:

eq_F_in_molw.show()


# In[ ]:

eq_Finv_Finmol = F_in_v == (F_in_mola + F_in_molw)*R_mol*T_in/P_a
print latex(eq_Finv_Finmol)
units_check(eq_Finv_Finmol)


# In[ ]:

var2('F_in_va_n', 'Volumetric inflow of dry air at 0oC and 101325 Pa', meter^3/second, latexname='F_{in,v,a,n}') 
var2('P_r', 'Reference pressure',pascal, value = 101325)
var2('T_r', 'Reference temperature', kelvin)

eq_Finmola_Finva_ref = fun_eq(F_in_mola == F_in_va_n * P_r/(R_mol*T_r))


# <p>To get $F_{in,v}$ and $F_{in,mol,w}$, we will consider that:</p>
# <p>$P_d = P_a - P_w$</p>
# <p>$P_w F_{in,v} = F_{in,mol,w} R_{mol} T_{in}$</p>
# <p>$(P_a - P_w) F_{in,v} = F_{in,mol,a} R_{mol} T_{a}$</p>

# In[ ]:

eq_Finmolw_Finv = F_in_molw == (P_w_in*F_in_v)/(R_mol*T_in)
print units_check(eq_Finmolw_Finv)
eq_Finv_Finmola = F_in_v == F_in_mola*R_mol*T_in/(P_a - P_w_in)
print units_check(eq_Finv_Finmola)
eq_Finmolw_Finmola_Pwa = fun_eq(eq_Finmolw_Finv.subs(eq_Finv_Finmola))
eq_Finv_Finva_ref = fun_eq(eq_Finv_Finmola.subs(eq_Finmola_Finva_ref))


# In[ ]:

print latex(eq_Finv_Finmola)
print latex(eq_Finv_Finva_ref)


# In[ ]:

vdict = cdict.copy()
vdict[F_in_va_n] = 10e-3/60   # 10 l/min reported by Cellkraft
vdict[T_d] = 273.15 + 10    # 10oC dew point
vdict[P_a] = 101325.
vdict[T_r] = 273.15
vdict[P_w_in] = eq_Pwin_Tdew.rhs().subs(vdict)
print vdict[P_w_in]
print vdict[F_in_va_n]

vdict[T_in] = 273.15+0 
inflow0 = eq_Finv_Finva_ref.rhs().subs(vdict)
print 'Volumentric flow at 0 oC: ' + str(inflow0) + ' m3/s'

vdict[T_in] = 273.15+25  
inflow25 = eq_Finv_Finva_ref.rhs().subs(vdict)
print 'Volumentric flow at 25 oC: ' + str(inflow25) + ' m3/s'


print '25oC/0oC: ' + str(inflow25/inflow0)

vdict[T_in] = 273.15+25 
vdict[P_w_in] = 0. 
inflow25 = eq_Finv_Finva_ref.rhs().subs(vdict)
print 'Volumentric flow at 25 oC without added vapour: ' + str(inflow25) + ' m3/s'


# <h2>Vapour pressure</h2>
# <p>The water fluxes associated with the incoming and the outgoing air according to the ideal gas law are $P_{v,in} F_{in,v}/(R_{mol} T_{in})$ and $P_{v,out}  F_{out,v}/(R_{mol} T_{out})$ respectively.</p>

# In[ ]:

eq_Foutv_Finv_Tout.subs(eq_F_in_molw)


# In[ ]:

eq_Foutv_Finv_Tout.show()


# In[ ]:

eq_F_out_v.show()


# In[ ]:

eq_F_in_mola.subs(eq_Finv_Finva_ref).show()


# In[ ]:

eq_Finv_Finva_ref.show()


# In[ ]:

eq_F_in_molw.show()
eq_F_out_molw.show()
eq_Foutmolw_Finmolw_Elmol.show()


# In[ ]:

eq1 = eq_F_out_molw.rhs() == eq_Foutmolw_Finmolw_Elmol.rhs().subs(eq_F_in_molw)
eq1.show()
print latex(eq1)
soln = solve(eq1, P_w_out)
print soln
eq_Pwout_Elmol = fun_eq(soln[0].subs(eq_Foutv_Finv_Tout.subs(eq_F_in_molw)).simplify_full())
print latex(eq_Pwout_Elmol)
units_check(eq_Pwout_Elmol)


# <p><span style="color: #ff0000;">It is a bit surprising that steady-state $P_{v_{out}}$ does not depend on $T_{out}$.</span></p>

# In[ ]:

soln[0].subs(eq_Foutv_Finv_Tout.subs(eq_F_in_molw)).simplify_full().show()
soln[0].subs(eq_F_out_v).subs(F_out_mola = F_in_mola, F_out_molw = F_in_molw + E_lmol*L_A).simplify_full().show()


# <p><span style="color: #ff0000;">The above are equivalent, because $F_{in,v} P_a = (F_{in,mol,a} + F_{in,mol,w}) R_{mol} T_{in}$<br /></span></p>

# In[ ]:

eq_Pwout_Elmol.subs(eq_Finv_Finva_ref).simplify_full().show()


# In[ ]:

show(soln[0])
show(eq_Foutv_Finv_Tout)
show(eq_F_in_molw)


# In[ ]:

show(soln[0].subs(eq_Foutv_Finv_Tout.subs(eq_F_in_molw)).simplify())


# In[ ]:

# T_out cancels out when the above is expanded
show(soln[0].subs(eq_Foutv_Finv_Tout.subs(eq_F_in_molw)).expand())


# <p>To convert from energetic to molar units, we need to divide $E_l$ by $\lambda_E M_w$:</p>

# In[ ]:

eq_Elmol_El = E_lmol == E_l/(lambda_E*M_w)
print units_check(eq_Elmol_El)
eq_El_Elmol = E_l == E_lmol*lambda_E*M_w


# In[ ]:

# In order to keep P_w_out = P_wa = const., we need to adjust P_w_in accordingly.
soln = solve(eq_Pwout_Elmol, P_w_in)
eq_Pwin_Elmol = soln[0].simplify_full()
show(eq_Pwin_Elmol)


# In[ ]:

vdict = cdict.copy()
vdict[F_in_va_n] = 10e-3/60   # 10 l/min reported by Cellkraft
vdict[T_d] = 273.15 + 10    # 10oC dew point
vdict[P_a] = 101325.
vdict[T_r] = 273.15
vdict[P_w_in] = eq_Pwin_Tdew.rhs().subs(vdict)
print vdict[P_w_in]
print vdict[F_in_va_n]

vdict[T_in] = 273.15+0 
inflow0 = eq_Finv_Finva_ref.rhs().subs(vdict)
print 'Volumentric flow at 0 oC: ' + str(inflow0) + ' m3/s'

vdict[T_in] = 273.15+25  
inflow25 = eq_Finv_Finva_ref.rhs().subs(vdict)
vdict[F_in_v] = eq_Finv_Finva_ref.rhs().subs(vdict)
print 'Volumentric flow at 25 oC: ' + str(inflow25) + ' m3/s'


print '25oC/0oC: ' + str(inflow25/inflow0)

vdict[T_in] = 273.15+25 
vdict[P_w_in] = 0. 
inflow25 = eq_Finv_Finva_ref.rhs().subs(vdict)
print 'Volumentric flow at 25 oC without added vapour: ' + str(inflow25) + ' m3/s'

vdict[E_l] = 100. # assuming 100 W/m2 El
vdict[L_A] = 0.03^2
vdict[E_lmol] = eq_Elmol_El.rhs().subs(vdict)
vdict[T_out] = 273+20. # Assuming 20oC T in chamber

eq_Pwout_Elmol.subs(vdict)


# <h2>Net radiation measurement</h2>
# <p>According to Incropera_fundamentals, Table 13.1, the view factor (absorbed fraction of radiation emitted by another plate) of a small plate of width $w_i$ at a distance $L$ from a parallel larger plate of width $w_j$ is calculated as:</p>

# In[ ]:

var2('L_s', 'Width of net radiation sensor', meter)
var2('L_ls', 'Distance between leaf and net radiation sensor', meter)
var2('F_s', 'Fraction of radiation emitted by leaf, absorbed by sensor', 1)
Wi = L_s/L_ls
Wj = L_l/L_ls
eq_Fs = F_s == (sqrt((Wi + Wj)^2 + 4) - sqrt((Wj - Wi)^2 + 4))/(2*Wi)
units_check(eq_Fs)
print latex(eq_Fs)


# In[ ]:

vdict = cdict.copy()
vdict[L_l] = 0.03
vdict[L_s] = 0.01
print eq_Fs.rhs().subs(vdict)(L_ls = 0.01)
plot(eq_Fs.rhs().subs(vdict), (L_ls, 0.001, 0.02))


# In[ ]:

eq_Rbalance


# In[ ]:

# Now looking at the fraction of radiation absorbed by the leaf that was emitted by the wall
vdict = cdict.copy()
vdict[L_l] = 0.1
vdict[L_s] = 0.03
print eq_Fs.rhs().subs(vdict)(L_ls = 0.01)
plot(eq_Fs.rhs().subs(vdict), (L_ls, 0.001, 0.03))


# In[ ]:

# ... and the sensor
vdict = cdict.copy()
vdict[L_l] = 0.1
vdict[L_s] = 0.01
print eq_Fs.rhs().subs(vdict)(L_ls = 0.01)
plot(eq_Fs.rhs().subs(vdict), (L_ls, 0.001, 0.02))


# In[ ]:

save_session('leaf_chamber_eqs.sobj')


# # Table of symbols

# In[ ]:

# Creating dictionary to substitute names of units with shorter forms
var('m s J Pa K kg mol')
subsdict = {meter: m, second: s, joule: J, pascal: Pa, kelvin: K, kilogram: kg, mole: mol}
variables = sorted([str(variable) for variable in udict.keys()],key=str.lower)
tabledata = [('Variable', 'Description (value)', 'Units')] + [(eval(variable),docdict[eval(variable)],(udict[eval(variable)]/eval(variable)).subs(subsdict)) for variable in variables]
table(tabledata)


# In[ ]:



