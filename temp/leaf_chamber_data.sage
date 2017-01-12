
# coding: utf-8

# #  Analysis of leaf chamber experiments
# This worksheet produces figures for:
# 
# **Technical note: An experimental setup to measure latent and sensible heat fluxes from (artificial) plant leaves.**
# 
# hess-2016-643
# 
# Author: Stan Schymanski (stan.schymanski@env.ethz.ch)

# ## Worksheet setup and importing equations and functions from other worksheets

# In[1]:

get_ipython().run_cell_magic(u'capture', u'storage', u"# The above redirects all output of the below commands to the variable 'storage' instead of displaying it.\n# It can be viewed by typing: 'storage()'\n\n# Setup of worksheet, incl. importing of packages and defining general custom functions\nload('temp/Worksheet_setup.sage')\nfrom shutil import copyfile   # for copying files between directories\nfrom matplotlib.ticker import MaxNLocator\nimport csv\nimport datetime\nimport time\nimport matplotlib.pyplot as plt\nimport matplotlib.dates as mdates")


# In[2]:

# From leaf_chamber_eqs, Leaf_enbalance2s_eqs and E_PM_eqs
load_session('temp/leaf_enbalance2s_eqs.sobj')
dict_vars1 = dict_vars.copy()

load_session('temp/E_PM_eqs.sobj')
dict_vars1.update(dict_vars)

load_session('temp/leaf_chamber_eqs.sobj')
dict_vars1.update(dict_vars)

dict_vars = dict_vars1.copy()
fun_loadvars(vardict=dict_vars)   # re-loading variable definitions


# In[3]:

T0


# In[4]:

path_figs = 'figures/'
path_data = 'data/'
path_data_orig = '/home/sschyman/Documents/STEP/Lab_data/leaf_chamber/'


# ## Functions to compute steady-state leaf energy balance components

# In[5]:

def fun_SS(vdict1):
    '''
    Steady-state T_l, R_ll, H_l and E_l under forced conditions.
    Parameters are given in a dictionary (vdict) with the following entries:
    a_s, a_sh, L_l, P_a, P_wa, R_s, Re_c, T_a, g_sw, v_w
    ''' 
    vdict = vdict1.copy()
    if not T_w in vdict1.keys():
        vdict[T_w] = vdict[T_a]


    # Nusselt number
    vdict[nu_a] = eq_nua.rhs().subs(vdict)
    vdict[Re] = eq_Re.rhs().subs(vdict)
    vdict[Nu] = eq_Nu_forced_all.rhs().subs(vdict)
    
    # h_c
    vdict[k_a] = eq_ka.rhs().subs(vdict)
    vdict[h_c] = eq_hc.rhs().subs(vdict)
 
    # gbw
    vdict[D_va] = eq_Dva.rhs().subs(vdict)
    vdict[alpha_a] = eq_alphaa.rhs().subs(vdict)
    vdict[rho_a] =  eq_rhoa.rhs().subs(vdict)
    vdict[Le] =  eq_Le.rhs().subs(vdict)
    vdict[g_bw] = eq_gbw_hc.rhs().subs(vdict)   
    
    # Hl, Rll
    vdict[R_ll] = eq_Rll.rhs().subs(vdict)
    vdict[H_l] = eq_Hl.rhs().subs(vdict)   

    # El
    vdict[g_tw] =  eq_gtw.rhs().subs(vdict)
    vdict[C_wa] = eq_Cwl.rhs()(P_wl = P_wa, T_l = T_a).subs(vdict)
    vdict[P_wl] = eq_Pwl.rhs().subs(vdict)
    vdict[C_wl] = eq_Cwl.rhs().subs(vdict)
    vdict[E_lmol] = eq_Elmol.rhs().subs(vdict)
    vdict[E_l] = eq_El.rhs().subs(eq_Elmol).subs(vdict)    

    # Tl
    try: vdict[T_l] = find_root((eq_Rs_enbal - R_s).rhs().subs(vdict), 273, 373)
    except: print 'something missing: ' + str((eq_Rs_enbal - R_s).rhs().subs(vdict))
    
    # Re-inserting T_l
    Tlss = vdict[T_l]
    for name1 in [C_wl, P_wl, R_ll, H_l, E_l, E_lmol]:
        vdict[name1] = vdict[name1].subs(T_l = Tlss)
    
    # Test for steady state
    if n((E_l + H_l + R_ll - R_s).subs(vdict))>1.:
        return 'error in energy balance: El + Hl + Rll - R_s = ' + str(n((E_l + H_l + R_ll - R_s).subs(vdict))) 
    return vdict


# In[6]:

from scipy.optimize import fsolve
def fun_SS1(vdict1):
    '''
    Steady-state T_lu, T_ll, R_lll, R_llu, H_ll, H_lu, E_lu and E_ll under forced conditions.
    Parameters are given in a dictionary (vdict) with the following entries:
    a_s, a_sh, k_l, z_l, L_l, P_a, P_wa, R_s, Re_c, T_a, g_swu, g_swl, v_w
    If g_sw is given instead of g_swu and g_swl, their values are deduced from a_s.
    ''' 
    vdict = vdict1.copy()
    if g_swl not in vdict.keys():
        if vdict[a_s] == 1:
            vdict[g_swl] = vdict[g_sw]
        if vdict[a_s] == 2:
            vdict[g_swl] = vdict[g_sw]/2

    if g_swu not in vdict.keys():
        if vdict[a_s] ==1:
            vdict[g_swu] = 0.
        if vdict[a_s] == 2:
            vdict[g_swu] = vdict[g_sw] - vdict[g_swu]
    vdict[C_wa] = eq_Cwl.rhs()(P_wl = P_wa, T_l = T_a).subs(vdict)
    vdict[nu_a] = eq_nua.rhs().subs(vdict)
    vdict[Re] = eq_Re.rhs().subs(vdict)
    vdict[Nu] = eq_Nu_forced_all.rhs().subs(vdict)
    vdict[k_a] = eq_ka.rhs().subs(vdict)
    vdict[h_cu] = eq_hcu.rhs().subs(vdict)
    vdict[h_cl] = eq_hcl.rhs().subs(vdict)
    vdict[alpha_a] = eq_alphaa.rhs().subs(vdict)
    vdict[D_va] = eq_Dva.rhs().subs(vdict)
    vdict[Le] =  eq_Le.rhs().subs(vdict)
    vdict[P_N2] = eq_PN2.rhs().subs(vdict)
    vdict[P_O2] = eq_PO2.rhs().subs(vdict)
    vdict[rho_a] =  eq_rhoa_Pwa_Ta.rhs().subs(vdict)
        
    eq1 = eq_enbalTlu.subs(eq_Elu, eq_Hlu, eq_Ql, eq_Rllu).subs(vdict)
    eq2 = eq_enbalTll.subs(eq_Ell, eq_Hll, eq_Ql, eq_Rlll).subs(vdict)
    def ff(v):
        T_lu, T_ll = map(float, v)
        return [eq1(T_lu=T_lu, T_ll=T_ll), eq2(T_lu=T_lu, T_ll=T_ll)]

    soln = fsolve(ff, [300, 300])
    
    
    vdict[T_lu] = soln[0]
    vdict[T_ll] = soln[1]
    vdict[Q_l] = eq_Ql.rhs().subs(vdict)
    vdict[E_lu] = eq_Elu.rhs().subs(eq_Pwl(T_l = T_lu)).subs(vdict)
    vdict[E_ll] = eq_Ell.rhs().subs(eq_Pwl(T_l = T_ll)).subs(vdict)
    vdict[H_lu] = eq_Hlu.rhs().subs(vdict)
    vdict[H_ll] = eq_Hll.rhs().subs(vdict)
    vdict[R_llu] = eq_Rllu.rhs().subs(vdict)
    vdict[R_lll] = eq_Rlll.rhs().subs(vdict)
    vdict[E_l] = vdict[E_lu] + vdict[E_ll]
    vdict[H_l] = vdict[H_lu] + vdict[H_ll]
    vdict[R_ll] = vdict[R_llu] + vdict[R_lll]

    
    # Test for steady state
    eq1 = eq_enbalTlu
    if n(eq1.subs(vdict))>1.:
        return 'error in energy balance: ' + str(eq1) + ' = ' + str(n(eq1.subs(vdict))) 
    eq1 = eq_enbalTll
    if n(eq1.subs(vdict))>1.:
        return 'error in energy balance: ' + str(eq1) + ' = ' + str(n(eq1.subs(vdict))) 
    
    eq1 = R_s - R_llu - R_lll - H_lu - H_ll - E_lu - E_ll
    if n(eq1.subs(vdict))>1.:
        return 'error in energy balance: ' + str(eq1) + ' = ' + str(n(eq1.subs(vdict)))     
            
    return vdict


# In[7]:

# Test using data from Fig. 8 in Ball et al. 1988
vdict = cdict.copy()
vdict[a_s] = 1
vdict[L_l] = 0.07
vdict[P_a] = 101325
vdict[P_wa] = 20/1000*101325
vdict[R_s] = 400
vdict[Re_c] = 3000
vdict[T_a] = 273+30
vdict[T_w] = vdict[T_a]
vdict[g_sw] = 0.15/40
vdict[v_w] = 1.
resdict = fun_SS(vdict)

vdict[k_l] = 0.3  # between 0.268 and 0.573 according to Hays_1975_The_thermal.pdf
vdict[z_l] = 0.0004
resdict1 = fun_SS1(vdict)

dict_print(resdict, list_names = [H_l, E_l, T_l])
print '2 sides:'
dict_print(resdict1, list_names = [H_lu, H_ll, E_lu, E_ll, T_lu, T_ll])


# ### Model by Ball et al., 1988

# In[8]:

# for model by Ball et al. 1988
var2('g_svmol', 'Stomatal condutance to vapour in molar units', mole/meter^2/second)
var2('c_pmol', 'Molar heat capacity of air', value = 29.2, units = joule/mole/kelvin)   # Units in the appendix are wrongly given as J mol/K
var2('L_E', 'Latent heat of vaporisation of water', value = 44000, units = joule/mole)
var2('r_bstar', 'Boundary layer resistance to heat transfer from unit area of leaf', meter^2*second/mole)
var2('r_b', 'Boundary layer resistnace to vapour transfer', meter^2*second/mole)
var2('r_s', 'Stomatal resistance to vapour transfer', meter^2*second/mole)

eq_Hl_Ball = H_l == c_pmol*(T_l - T_a)/r_bstar
eq_El_Ball = E_l == L_E*(P_wl - P_wa)/P_a/(r_s + r_b)
eq_rbstar = r_bstar == (3.8*L_A^(1/4)*v_w^(-1/2))
eq_rb = r_b == (1.78/a_s*r_bstar)
print units_check(eq_Hl_Ball).simplify_full()
print units_check(eq_El_Ball).simplify_full()


# In[9]:

def fun_SS_Ball(vdict1):
    '''
    Steady-state T_l, h_c, g_bv, g_tv, R_ll, H_l and E_l under forced conditions.
    h_c and g_bv are calculated using Appendix in Ball et al. (1988).
    see Ball_1988_Maintenance_of_Leaf2.pdf
    Parameters are given in a dictionary (vdict) with the following entries:
    a_s, L_l, P_a, P_wa, R_s, Re_c, T_a, g_svmol, v_w
    ''' 
    vdict = vdict1.copy()
    if not L_A in vdict.keys():
        vdict[L_A] = (pi()*L_l^2).subs(vdict)
    if not T_w in vdict.keys():
        vdict[T_w] = vdict[T_a]

    vdict[r_s] = (1/(40*g_sw)).subs(vdict)
    try:
        vdict[r_bstar] = eq_rbstar.rhs().subs(vdict).n()   #two-sided resistance for sensible heat flux
    except:
        vdict[r_bstar] = eq_rbstar.rhs().subs(vdict)
        print 'r_bstar = ' + str(vdict[r_bstar])
    vdict[r_b] = eq_rb.rhs().subs(vdict)             # one-sided resistance to latent heat flux
    
    
    
    Rll = eq_Rll.rhs().subs(vdict)
    Hl = eq_Hl_Ball.rhs().subs(vdict)
    El = eq_El_Ball.rhs().subs(eq_Pwl).subs(vdict)
    
    enbal = El + Hl + Rll - R_s == 0 
    #print enbal(R_ll = Rll, H_l = Hl, E_l = El).subs(vdict)
    Tss = find_root(enbal(R_ll = Rll, H_l = Hl, E_l = El).subs(vdict), 273, 373)
    
    Rll1 = Rll(T_l = Tss)
    Hl1= Hl(T_l = Tss)
    El1 = El(T_l = Tss)

    # Test for steady state
    if (El1 + Hl1 + Rll1 - vdict[R_s])>1.:
        print (El, Hl, Rll, vdict[R_s])
        print Tss
        return 'error in energy balance'
    Pwl = eq_Pwl.rhs()(T_l = Tss).subs(vdict)
    
    #dict_print(vdict)
    return {'Tl_ball':n(Tss), 'Rll_ball':n(Rll1), 'Hl_ball':n(Hl1), 'El_ball':n(El1), 'rs_ball': vdict[r_s], 'rbstar_ball': vdict[r_bstar], 'rb_ball': vdict[r_b]}


# In[10]:

# Fig. 8 in Ball et al. 1988
vdict = cdict.copy()
vdict[a_s] = 1
vdict[L_l] = 0.07
vdict[P_a] = 101325
vdict[P_wa] = 20/1000*101325
vdict[R_s] = 400
vdict[Re_c] = 3000
vdict[T_a] = 273+30
vdict[T_w] = vdict[T_a]
vdict[g_sw] = 0.15/40

vdict[v_w] = 1.
ssdict = fun_SS(vdict)
#print ssdict
balldict = fun_SS_Ball(vdict)
#print balldict
print 'h_c(Ball): ' + str((c_pmol/balldict['rbstar_ball']).subs(vdict))
print 'g_vmol(Ball): ' + str((1/(r_b + r_s))(r_b = balldict['rb_ball'], r_s = balldict['rs_ball']))
print 'g_vmol(SS): ' + str(eq_gtwmol_gtw(T_l = T_a).subs(eq_gtw).subs(ssdict).subs(vdict))
print 'T_l(Ball): ' + str(balldict['Tl_ball'])
print 'T_l(SS): ' + str(ssdict[T_l])
print 'T_a: ' + str(vdict[T_a])


# <p><span style="color: #ff0000;">According to Fig. 8 in Ball et al., 1988, steady-state leaf temperature should be higher by 10 K than air temperature! Here it is only 6. <br /></span></p>

# ### Analytical models from Schymanski & Or, 2016

# In[11]:

def fun_SS_PM(vdict1, twosides = False):
    '''
   Analytical equations from Worksheet E_PM_eqs, including detailed steady-state.
    '''
    vdict = vdict1.copy()
    if not L_A in vdict1.keys():
        vdict[L_A] = (pi()*L_l^2).subs(vdict1)
    if not T_w in vdict1.keys():
        vdict[T_w] = vdict[T_a]
    if not P_wa in vdict1.keys():
        print 'P_wa is missing'
    ss = fun_SS(vdict)
    if twosides:
        ss2 = fun_SS1(vdict)

    vdict[P_was] = eq_Pwl.rhs()(T_l = T_a).subs(vdict) 
    vdict[Delta_eTa] = eq_Deltaeta_T.rhs().subs(vdict)
    vdict[k_a] = eq_ka.rhs().subs(vdict)
    vdict[nu_a] = eq_nua.rhs().subs(vdict)
    vdict[Re] = eq_Re.rhs().subs(vdict)
    vdict[Nu] = eq_Nu_forced_all.rhs().subs(vdict)
    vdict[h_c] = eq_hc.rhs().subs(vdict) 
    
    vdict[P_N2] = eq_PN2.rhs().subs(vdict)
    vdict[P_O2] = eq_PO2.rhs().subs(vdict)
    vdict[alpha_a] = eq_alphaa.rhs().subs(vdict)
    vdict[k_a] = eq_ka.rhs().subs(vdict)
    vdict[D_va] = eq_Dva.rhs().subs(vdict)
    vdict[Le] = eq_Le.rhs().subs(vdict)
    vdict[rho_a] = eq_rhoa_Pwa_Ta.rhs().subs(vdict) 
    vdict[g_bw] = eq_gbw_hc.rhs().subs(vdict)
    vdict[g_tw] = eq_gtw.rhs().subs(vdict)
    vdict[g_twmol] = eq_gtwmol_gtw_iso.rhs().subs(vdict)

    # Generalized Penman, getting Rll first
    vdict_GPRll = vdict.copy()
    vdict_GPRll[R_ll] = 0.
    vdict_GPRll[T_l] = eq_Tl_Delta.rhs().subs(eq_ce_conv, eq_ch_hc).subs(vdict_GPRll)
    vdict_GPRll[R_ll] = eq_Rll.rhs().subs(vdict_GPRll)
    vdict_GPRll[E_l] = eq_El_Delta.rhs().subs(eq_ce_conv, eq_ch_hc).subs(vdict_GPRll)
    vdict_GPRll[H_l] = eq_Hl_Delta.rhs().subs(eq_ce_conv, eq_ch_hc).subs(vdict_GPRll)


    # Using linearised R_ll solution
    vdict_lin = vdict.copy()
    namesdict = [E_l, H_l, T_l, R_ll]
    vdict_lin[E_l] = eq_El_Delta_Rlllin.rhs().subs(eq_ce_conv, eq_ch_hc).subs(vdict_lin)
    vdict_lin[H_l] = eq_Hl_Delta_Rlllin.rhs().subs(eq_ce_conv, eq_ch_hc).subs(vdict_lin)
    vdict_lin[T_l] = eq_Tl_Delta_Rlllin.rhs().subs(eq_ce_conv, eq_ch_hc).subs(vdict_lin)
    vdict_lin[R_ll] = eq_Rll_tang.rhs().subs(vdict_lin)


    # 'R_N = R_s:' MU-book, P. 79, under cloudy skies, implying that R_ll = 0
    vdict2 = vdict.copy()
    vdict2[R_ll] = 0
    
    # Generalized Penman
    vdict_GP = vdict2.copy()
    vdict_GP[E_l] = eq_El_Delta.rhs().subs(eq_ce_conv, eq_ch_hc).subs(vdict_GP)
    vdict_GP[H_l] = eq_Hl_Delta.rhs().subs(eq_ce_conv, eq_ch_hc).subs(vdict_GP)
    vdict_GP[T_l] = eq_Tl_Delta.rhs().subs(eq_ce_conv, eq_ch_hc).subs(vdict_GP)


    # Penman-stomata
    vdict_PS = vdict2.copy()
    vdict_PS[S] = eq_S_gbw_gsw.rhs().subs(vdict_PS)
    vdict_PS[f_u] = eq_fu_gbw.rhs().subs(vdict_PS)
    vdict_PS[gamma_v] = eq_gammav_as.rhs().subs(vdict_PS)
    vdict_PS[E_l] = eq_El_P52.rhs().subs(vdict_PS)
    vdict_PS[H_l] = eq_Hl_P52.rhs().subs(vdict_PS)
    vdict_PS[T_l] = eq_Tl_P52.rhs().subs(vdict_PS)
    
    # PM equation
    vdict_PM = vdict2.copy() 
    vdict_PM[r_s] = 1/vdict_PM[g_sw]
    vdict_PM[r_a] = eq_ra_hc.rhs().subs(vdict_PM)
    vdict_PM[gamma_v] = eq_gammav_MU.rhs().subs(vdict_PM)
    vdict_PM[epsilon] = eq_epsilon.rhs().subs(vdict_PM)
    vdict_PM[E_l] = eq_El_PM2.rhs().subs(vdict_PM)
    vdict_PM[H_l] = (R_s - R_ll - E_l).subs(vdict_PM)
    
    # MU equation
    vdict_MU = vdict2.copy() 
    vdict_MU[n_MU] = (a_sh/a_s).subs(vdict_MU)
    vdict_MU[r_s] = 1/vdict_MU[g_sw]
    vdict_MU[r_a] = eq_ra_hc.rhs().subs(vdict_MU)
    vdict_MU[gamma_v] = eq_gammav_MU.rhs().subs(vdict_MU)
    vdict_MU[E_l] = eq_El_MU2.rhs().subs(vdict_MU)
    vdict_MU[H_l] = (R_s - R_ll - E_l).subs(vdict_MU)
    
    # 'Corrected MU-equation: '
    vdict_MUc = vdict2.copy()
    vdict_MUc[r_s] = 1/vdict_MUc[g_sw]
    vdict_MUc[r_a] = eq_ra_hc.rhs().subs(vdict_MUc)    
    vdict_MUc[gamma_v] = eq_gammav_MU.rhs().subs(vdict_MUc)
    vdict_MUc[E_l] = eq_El_MU_corr.rhs().subs(vdict_MUc)
    vdict_MUc[H_l] = (R_s - R_ll - E_l).subs(vdict_MUc)

    resdict = ss.copy()
    str1 = 'GPRll'
    namesdict = [E_l, H_l, T_l, R_ll]
    for name1 in namesdict:
        resdict[str1+'_'+str(name1)] = eval('vdict_'+str1)[name1]    
    
    str1 = 'lin'
    namesdict = [E_l, H_l, T_l, R_ll]
    for name1 in namesdict:
        resdict[str1+'_'+str(name1)] = eval('vdict_'+str1)[name1]

    str1 = 'GP'
    namesdict = [E_l, H_l, T_l]
    for name1 in namesdict:
        resdict[str1+'_'+str(name1)] = eval('vdict_'+str1)[name1]

    str1 = 'PS'
    namesdict = [E_l, H_l, T_l]
    for name1 in namesdict:
        resdict[str1+'_'+str(name1)] = eval('vdict_'+str1)[name1]    

    str1 = 'PM'
    namesdict = [E_l, H_l]
    for name1 in namesdict:
        resdict[str1+'_'+str(name1)] = eval('vdict_'+str1)[name1]

    str1 = 'MU'
    namesdict = [E_l, H_l]
    for name1 in namesdict:
        resdict[str1+'_'+str(name1)] = eval('vdict_'+str1)[name1]        

    str1 = 'MUc'
    namesdict = [E_l, H_l]
    for name1 in namesdict:
        resdict[str1+'_'+str(name1)] = eval('vdict_'+str1)[name1]

    if twosides:
        str1 = '2s'
        namesdict = ss2.keys()
        for name1 in namesdict:
            resdict[str1+'_'+str(name1)] = ss2[name1]
            
    return resdict


# # Data reading, computing and display functions

# ## Functions to read data

# In[12]:

def fun_read_csv(fname1, lc_datetime_format = "%Y-%m-%d %H:%M:%S", lc_timepos = [], split1 = ';'):
    '''
    Reads file written by leaf_chamber_read_data and saves them as lc_data
    If csv file fname1 does not exist in path_data, copies the file from path_data_orig to path_data
    before opening and reading into numpy array. Contents of data file are printed as table before returning numpy array.
    '''
    lc_timepos = lc_timepos[:]
    fname = path_data + fname1
    try:
        reader = csv.reader(open(fname, 'r'), delimiter=split1, quotechar='"')
    except:
        copyfile(path_data_orig + fname1, fname)
        reader = csv.reader(open(fname, 'r'), delimiter=split1, quotechar='"')
    htmldata = []
    lc_nameslist = reader.next()
    htmldata.append(lc_nameslist)
    htmldata[-1][-1] = htmldata[-1][-1]+ '...................................................'  # To avoid very thick rows
    #print htmldata
    lc_formatslist = ['S200' for i in srange(len(lc_nameslist))]
    lc_unitslist = reader.next()
    htmldata.append(lc_unitslist)
    #print lc_unitslist
    #print reader.next()
    ncols = len(lc_nameslist)
    lc_datetime_format = "%Y-%m-%d %H:%M:%S"
    csvdata = []
    with open(fname, 'r') as file_out:
        rows = file_out.readlines()
       
    # Determining dtypes
    for row in rows[2:]:
        row1 = row.strip('\r\n').split(split1)
        for i in srange(len(row1)-1):
            try:
                blah = float(row1[i])
                lc_formatslist[i] = 'float'
            except:
                try:
                    blah = datetime.datetime.fromtimestamp(time.mktime(time.strptime(row1[i].strip(), lc_datetime_format)))
                    lc_timepos.append(i)
                    lc_formatslist[i] = 'datetime64[us]'
                except:
                    lc_formatslist[i] = 'S200'
    
    lc_timepos = list(set(lc_timepos))     # to get unique values, i.e. drop repeated positions 

    for row in rows[2:]:
        row1 = row.strip('\r\n').split(split1)
        htmldata.append(row1[:])
        for i in lc_timepos:
            try:
                row1[i] = datetime.datetime.fromtimestamp(time.mktime(time.strptime(row1[i].strip(), lc_datetime_format)))
            except:
                print lc_timepos
                print row1
        row2 = tuple(row1) 
        csvdata.append(row2)
    
    
    try:
        lc_data = np.array(csvdata,dtype = zip(lc_nameslist,lc_formatslist))
    except: 
        print 'Error in dtype'
        return csvdata
    pretty_print(table(htmldata, header_row=True, align='center'))
    return lc_data


# In[13]:

def fun_read_campbell(fname1, nr_datetime_format = "%Y-%m-%d %H:%M:%S.%f", datelen = 21, datelast = '.0'):
    '''
    Reads Campbell data logger file and returns numpy array. 
    If file fname1 does not exist in path_data, copies the file from path_data_orig to path_data
    before opening.
    '''
    import sys, traceback
    fname = path_data + fname1
    try:
        reader = csv.reader(open(fname, 'rb'), delimiter=',', quotechar='"')
    except:
        copyfile(path_data_orig + fname1, fname)
        reader = csv.reader(open(fname, 'rb'), delimiter=',', quotechar='"')  
    
    print reader.next()
    nr_nameslist = reader.next()
    print nr_nameslist
    nr_unitslist = reader.next()
    #print nr_unitslist
    reader.next()
    ncols = len(nr_nameslist)
      
    csvdata = []
    contd = True
    while contd:
        try:
            row1 = reader.next()
            date1 = row1[0]
            # Need to add milliseconds if missing
            if len(date1) < datelen: date1 = date1+datelast
            row1[0] = datetime.datetime.fromtimestamp(time.mktime(time.strptime(date1.strip(), nr_datetime_format)))
            row2 = tuple(row1)
            csvdata.append(row2)
        except StopIteration:
            contd = False
        except:
            traceback.print_exc(file=sys.stdout)
            print row1
            #contd = False
    
    #print csvdata
    nr_formatslist = ['datetime64[us]', 'int']
    for i in srange(len(nr_nameslist) - 2):
            nr_formatslist.append('float')
    nr_data = np.array(csvdata,dtype = zip(nr_nameslist,nr_formatslist))
    return nr_data


# In[14]:

def fun_read_li(li_fname, prefix1 = ['1900-01-01_']):
    '''
    Reads Licor data file and returns as numpy array. If data contains several days, add a prefix for each day.
    If file li_fname does not exist in path_data, copies the file from path_data_orig to path_data
    before opening.
    '''
    time_format = "%H:%M:%S"
    datetime_format = "%Y-%m-%d %H:%M:%S"
    fname = path_data + li_fname
    try:
        reader = csv.reader(open(fname, 'rb'), delimiter='\t')
    except:
        copyfile(path_data_orig + li_fname, fname)
        reader = csv.reader(open(fname, 'rb'), delimiter='\t')

    
    row = ['']
    csvdata = []
    outdata = []
    nameflag = 'n'
    dataflag = 'n'
    prefpos = 0
    timeold = 0
    for row in reader:
        #print row
        if dataflag == 'n':
            if row[0][0]!='<':
                print row
        if row == ['$STARTOFDATA$']:
            print ' '
            print 'NEW data set'  
            print 'length previous = '+str(len(csvdata))
            if len(csvdata)>1:    # Converting csvdata to np.array and adding to outdata
                    li_formatslist = ['int','datetime64[us]']
                    for i in srange(len(li_nameslist) - 3):
                            li_formatslist.append('float')
                    li_formatslist.append('str')
                    li_data = np.array(csvdata,dtype = zip(li_nameslist,li_formatslist))
                    outdata.append(li_data)

            # Starting new data series
            csvdata = []
            nameflag = 'y'  
  

        if len(row)>1:
            if nameflag == 'y':
                li_nameslist = row
                nameflag = 'n'
                dataflag = 'y'
            elif dataflag == 'y':
                row1 = row[:]
                prefix = prefix1[prefpos]
                TS = prefix[0:10] + " " + row[1].strip()
                try:
                    row1[1] = datetime.datetime.fromtimestamp(time.mktime(time.strptime(TS, datetime_format)))
                    if timeold!=0:
                        if row1[1]<timeold:    # increment to next date if new time smaller old time
                            prefpos = prefpos + 1
                            TS = prefix[0:10] + " " + row[1].strip()
                            row1[1] = datetime.datetime.fromtimestamp(time.mktime(time.strptime(TS, datetime_format)))
                    timeold = row1[1]
                    row3 = tuple(row1)
                    csvdata.append(row3)        
                except:
                    print 'failed'
        elif dataflag == 'y':
            print row                    
                    
                
    li_formatslist = ['int','datetime64[us]']
    for i in srange(len(li_nameslist) - 3):
            li_formatslist.append('float')
    li_formatslist.append('str')
    li_data = np.array(csvdata,dtype = zip(li_nameslist,li_formatslist))
    outdata.append(li_data)
    print 'number of data sets: '+str(len(outdata))
    return outdata


# ## Equations to infer conductances from flux measurements

# In[15]:

eq_gsw_gtw = solve(eq_gtw, g_sw)[0]
units_check(eq_gsw_gtw).simplify_full()


# In[16]:

eq_gtw_Elmol = solve(eq_Elmol.subs(eq_Cwa, eq_Cwl), g_tw)[0]
units_check(eq_gtw_Elmol)


# In[17]:

eq_hc_Hl = solve(eq_Hl, h_c)[0]
units_check(eq_hc_Hl)


# ## Functions to automate computation of derived results and plotting

# In[18]:

def fun_results(lc_data, vdict1 = {}, poslist = [], nameRs = '', ndict1 = {}, Pwinoffset = 0, Tdewfac = 1.06,                 Tdewoffset = -2.45, twosides=False):
    '''
    Simulates steady-state conditions for every data set contained in lc_data
    and returns a numpy array with input data and results. Mapping of fields in lc_data
    is given in ndict, e.g.
     ndict = {R_s: '', T_d: 'T_dew', T_l: 'T_leaf_in', Q_in: 'fan_power_avg', S_a: 'Rn_above_leaf', S_b: 'Rn_below_leaf', S_s: 'Rn_beside_leaf', T_in: 'T_in1', F_in_v: 'Air_inflow', v_w: 'Wind', T_a: 'T_chamber', E_lmol: 'water_flow_avg'}
    If vdict1 is not given or values are missing, predefined values are assumed. Pwinoffset allows accounting for bias in Cellkraft P_w_in by negative 50 to negative 200 Pa, see LI6400_data_repaired2.sws.
    Usage:
        sage: fun_results(lc_data, vdict1 = {}, poslist = [1,2,4])
    '''

    poslist = poslist[:]
    def fun_subs(expr, vdict):
        '''
        Substitutes vdict in equation and returns result. 
        In case of a ValueError, e.g. because an entry is NaN,
        returns NaN, rather than interrupting because of exception.
        '''
        try:
            return float(expr.subs(vdict))
        except:
            #print expr.subs(vdict)
            return float('NaN')
    
    
    ndict = {R_s: '', T_d: 'T_dew', T_l: 'T_leaf_in', Q_in: 'fan_power_avg', S_a: 'Rn_above_leaf', S_b: 'Rn_below_leaf', S_s: 'Rn_beside_leaf',              T_in: 'T_in1', F_in_v: 'Air_inflow', v_w: 'Wind', T_a: 'T_chamber', E_lmol: 'water_flow_avg'}
    if len(ndict1)>0:
        ndict = fun_dict(ndict,ndict1)
    # Standard values
    vdict = cdict.copy()
    vdict[L_l] = 0.03  # 3x3 cm area
    vdict[a_s] = 1    # one sided stomata
    vdict[alpha_l] = 0  # assuming 0 albedo
    vdict[g_sw] = 999   # unlimiting, filter paper only
    vdict[P_a] = 101325
    vdict[Re_c] = 3000
    vdict[R_s] = 0
    vdict[T_r] = vdict[T0]
    vdict[P_a] = 101325
    if len(vdict1)>0:
        vdict = fun_dict(vdict,vdict1)
    if L_A not in vdict.keys():
        vdict[L_A] = vdict[L_l]^2   
    if R_d in vdict.keys():
        vdict[R_s] = eq_Rs_Rd.rhs().subs(vdict) 
    #print vdict
    results = []
    allinput = []
    
    if len(poslist) == 0:
        poslist = srange(len(lc_data))
    
    for i in poslist:
        rdict1 = {}   # For storing results that are not in vdict, e.g. if keys are strings rather than variables.        
        
        # Tdew reported by Cellkraft is biased (see worksheet Cellkraft_tests). Corrections range between y=1.09*x - 2.42 and y=1.06*x - 2.45
        Tdew = Tdewfac*lc_data[ndict[T_d]][i]+Tdewoffset+vdict[T0]
        if P_w_in in ndict1:
            vdict[P_w_in] = lc_data[ndict[P_w_in]][i]
        else:
            vdict[P_w_in] = eq_Pwl.rhs().subs(T_l = Tdew).subs(vdict) + Pwinoffset
        Qin = lc_data[ndict[Q_in]][i]
        Tlmeas = lc_data[ndict[T_l]][i]+vdict[T0]
        try:
            Rn_above = lc_data[ndict[S_a]][i]
            Rn_below = lc_data[ndict[S_b]][i]
            #Rn_beside = abs(lc_data[ndict[S_s]][i])
            Rn_leaf = eq_Rbalance.rhs()(S_a = Rn_above, S_b = Rn_below)
            rdict1['Rn_leaf'] = Rn_leaf   
        except:
            pass
        if len(ndict[R_s])>0:
            vdict[R_d] = abs(lc_data[ndict[R_s]][i]) 
            vdict[R_s] = eq_Rs_Rd.rhs().subs(vdict)    
            #print 'R_s from ' + nameRs   
        
        vdict[T_in] = lc_data[ndict[T_in]][i]+vdict[T0]
        vdict[F_in_va_n] = lc_data[ndict[F_in_v]][i]/1000/60
        vdict[F_in_v] = eq_Finv_Finva_ref.rhs().subs(vdict)
        vdict[F_in_mola] = eq_Finmola_Finva_ref.rhs().subs(vdict)
        vdict[F_in_molw] = eq_Finmolw_Finmola_Pwa.rhs().subs(vdict)

        #print 'F_in_v = ' + str(vdict[F_in_v])
        #print vdict[P_v_in]
    
        vdict[v_w] = lc_data[ndict[v_w]][i]
        vdict[T_a] = lc_data[ndict[T_a]][i]+ 273.15
        if T_w not in vdict1.keys():
            vdict[T_w] = vdict[T_a]
        
        Elmolmeas = lc_data[ndict[E_lmol]][i]*1e-6/60/vdict[L_A]/vdict[M_w]
        Elmeas = eq_El_Elmol.rhs()(E_lmol = Elmolmeas).subs(vdict)
        vdict[P_wa] = eq_Pwout_Elmol.rhs().subs(E_lmol = Elmolmeas).subs(vdict)
        vdict[F_out_molw] = eq_Foutmolw_Finmolw_Elmol.rhs().subs(E_lmol = Elmolmeas).subs(vdict)
        #print 'P_w_out = ' + str(vdict[P_wa])
        #Hlmeas = eq_Hl_Tin_Tout.rhs().subs(E_lmol = Elmolmeas, P_a = 101325, T_out = T_a, Q_in = Qin).subs(vdict)  
        Hlmeas = eq_Hl_Tin_Tout_Fmol.rhs().subs(E_lmol = Elmolmeas, P_a = 101325, T_out = T_a, Q_in = Qin).subs(vdict)  
        
        

        Balldict = fun_SS_Ball(vdict)
        PMdict = fun_SS_PM(vdict, twosides=twosides)
        
        # Inferring g_tw and g_sw from Elmolmeas, Tlmeas, P_wa and g_bw
        vdict1 = vdict.copy()
        vdict1[E_lmol] = Elmolmeas
        vdict1[T_l] = Tlmeas
        vdict1[P_wl] = eq_Pwl.rhs().subs(vdict1)
        vdict1[g_bw] = PMdict[g_bw]
        vdict1[g_tw] = fun_subs(eq_gtw_Elmol.rhs(), vdict1)  
        vdict1[g_sw] = fun_subs(eq_gsw_gtw.rhs(),vdict1)
 
        
        # Inferring h_c from Hlmeas, Tlmeas and Tameas
        vdict1[H_l] = Hlmeas
        vdict1[h_c] = fun_subs(eq_hc_Hl.rhs(),vdict1)
        
           
        #resdict = dict(vdict.items() + SSdict.items() + rdict1.items() + Balldict.items() + PMdict.items())
        resdict = dict(vdict.items() + rdict1.items() + Balldict.items() + PMdict.items())
        resdict['Tlmeas'] = Tlmeas
        resdict['Elmeas'] = Elmeas
        resdict['Elmolmeas'] = Elmolmeas
        resdict['Hlmeas'] = Hlmeas
        resdict['g_twmeas'] = vdict1[g_tw]
        resdict['g_swmeas'] = vdict1[g_sw]
        resdict['h_cmeas'] = vdict1[h_c]
        results.append(tuple(list(lc_data[i]) + resdict.values()))
        allinput.append(vdict)
    #print resdict    
    names = [blah[0] for blah in lc_data.dtype.descr] + resdict.keys()
    nameslist = [str(names[i]) for i in srange(len(names))]
    formatslist = [blah[1] for blah in lc_data.dtype.descr] + ['f8' for i in srange(len(nameslist))]
    #print results
    #print zip(nameslist,formatslist)
    try:
        results = np.array(results,dtype = zip(nameslist,formatslist))
    except:
        print results
        results1 = np.nan_to_num(results)  # Converts any NaN to something that can be expressed as a float
        print nameslist
        print formatslist
        print vdict
        print results
        return results
   
    return results


# In[19]:

def fun_plot_TN(results1, varname1 = 'v_w', axeslabel1 = 'Wind speed (m s$^{-1}$)', EHTmods = None,                    Emods = [('E_l', '(mod.)', '-')],                    Hmods = [('H_l', '(mod.)', '-')],                    Tmods = [('T_l', '(bulk)', '-')],                    energy=True, esums = True, leaftemp = True,                    alltemps = [('T_leaf1', 'TC1', 'v'), ('T_leaf2', 'TC2', '^'), ('T_leaf_IR', 'TCIR', '.')],                    Hobs = True, gtwcomp = True, hccomp = True, hccomp1 = False, font_size = 25, fsize = 14, lfsize = 28,                   axfsize = 22, psize = 12, dpi1 = 300, leglength = 2, lwidth = 2, fname = False, xfac = 1,                   fext = '.png', width = 12, height = 9, scaling = 0.6, hspace1 = 0., pltitle = False, lrot = 60, ylaboffs = -0.17):
    '''
    Sorts results1 for variable varname1 and plots diagnostic plots 
    for energy, esums, leaftemp, alltemps (set either of them to False
    if not desired). 
    '''
    
    list_labels = ['0', '(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
    #### Choice of data to be plotted
    if EHTmods == 'Ball':
        Emods = [('E_l', '(mod.)', '-'), ('El_ball', '(Ball)', '--')]
        Hmods = [('H_l', '(mod.)', '-'), ('Hl_ball', '(Ball)', '--')]
        Tmods = [('T_l', '(bulk)', '-'), ('Tl_ball', '(Ball)', '--')]
    if EHTmods == '2s':
        Emods = [('E_l', '(mod.)', '-'), ('2s_E_l', '(2s)', '--')]
        Hmods = [('H_l', '(mod.)', '-'), ('2s_H_l', '(2s)', '--')]
        Tmods = [('T_l', '(bulk)', '-'), ('2s_T_lu', '(upper side)', '--'), ('2s_T_ll', '(lower side)', ':')]
    
    #### Computations
    # Sorting array along v_w
    results2 = results1.copy()
    results2 = np.sort(results2, order = varname1)
    pos_vw = srange(len(results2))
    xdata = xfac*results2[varname1][pos_vw]  # unit conversions can be defined in xfac
    Talist = results2['T_a'][pos_vw]
    
    #### Setup plot
    ylabelxposmin = 0
    nrows = 0
    if energy: nrows = nrows+1
    if esums: nrows = nrows+1
    if leaftemp: nrows = nrows+1
    if gtwcomp: nrows = nrows+1
    if hccomp: nrows = nrows+1


    height = height*nrows
    ncols = 1
    fsize = [scaling*width,scaling*height]

    pylab.close('all')
    fig = plt.figure(figsize=fsize, dpi = dpi1)
    pylab.subplots_adjust(hspace=hspace1)  # set to very small value >0 if you want plots to merge

    ### Plotting
    nplot = 0
    if energy:  
        nplot = nplot+1
        ax2 = plt.subplot(nrows,ncols,nplot)
        ax2.plot(xdata,results2['Elmeas'][pos_vw], linestyle='', marker='o', markersize=psize, label = '$E_l$ (obs.)')

        for i in srange(len(Emods)):
            tup1 = Emods[i]
            if len(tup1)<4: 
                tup1 = tuple(list(tup1) + ['blue'])
            ax2.plot(xdata,results2[tup1[0]][pos_vw], color = tup1[3], linewidth = lwidth, linestyle=tup1[2], marker=None, label = '$E_l$ ' + tup1[1])

        if Hobs:
            ax2.plot(xdata,results2['Hlmeas'][pos_vw], linestyle='', marker='o', markerfacecolor = 'white', markeredgecolor = 'black', markersize=psize, label = '$H_l$ (obs.)')
        for i in srange(len(Hmods)):
            tup1 = Hmods[i]
            if len(tup1)<4: 
                tup1 = tuple(list(tup1) + ['black'])
            ax2.plot(xdata,results2[tup1[0]][pos_vw], color = tup1[3], linewidth = lwidth, linestyle=tup1[2], marker=None, label = '$H_l$ ' + tup1[1])
        #ax2.tick_params(pad=axfsize/2)
        ax2.margins(y=0.1)
        ylab = ax2.set_ylabel('Energy flux from leaf \n (W m$^{-2}$)') 
        ax2.get_yaxis().set_label_coords(ylaboffs,0.5)

        lgd = ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':lfsize})  # aim is to put legend on the right side of figure without shrinking the figure. Need to add bbox_extra_artists=(lgd,), bbox_inches='tight' to savefig
        # Setting axes font sizes
        ax2.yaxis.label.set_fontsize(lfsize)
        for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
            item.set_fontsize(axfsize)
        if pltitle: ax2.set_title(list_labels[nplot], fontsize = font_size)
        if nplot==1: 
            ax1 = ax2
        else:
            plt.draw()   # to generate tick labels
            list_tickpos = [blah.label.get_position()[-1] for blah in ax2.yaxis.get_major_ticks()]
            tickmax = max(list_tickpos)
            tick1 = ax2.yaxis.get_major_ticks()[1].label.get_position()[-1]    
            tick2 = ax2.yaxis.get_major_ticks()[2].label.get_position()[-1]  
            maxval = tickmax + 0.9*abs(tick2 - tick1)
            ax2.set_ylim(top=maxval)   # to avoid overlapping tick marks with plot above
        


    if esums: 
        nplot = nplot+1
        if hspace1<0.1:
            try: 
                pylab.setp(ax2.get_xticklabels(), visible=False)   # Remove x-axis labels from previous plot
            except:
                pass
        ax2 = plt.subplot(nrows,ncols,nplot, sharex = ax1)
        if Hobs:
            ax2.plot(xdata,results2['Elmeas'][pos_vw]+results2['Hlmeas'][pos_vw], linestyle='', marker='o', markersize=psize, label = '$E_l + H_l$ (obs.)')
        if 'Rn_leaf' in results2.dtype.names:
            ax2.plot(xdata,results2['Rn_leaf'][pos_vw], linestyle='', marker='o', markersize=psize, color = 'red', label = '$R_{nleaf}$ (obs.)')

        for i in srange(len(Hmods)):
            ax2.plot(xdata,results2[Emods[i][0]][pos_vw]+results2[Hmods[i][0]][pos_vw], color = 'blue', linewidth = lwidth, linestyle=Emods[i][2], marker=None, label = '$E_l + H_l$ ' + Emods[i][1])

        ax2.plot(xdata,results2['R_s'][pos_vw] - results2['R_ll'][pos_vw], linestyle=":", color = 'red', linewidth = lwidth, marker=None, label = '$R_s - R_{ll}$ (mod.)')
        #ax2.tick_params(pad=axfsize/2)
        ax2.margins(y=0.1)
        ax2.set_ylim(bottom=0)
        ax2.set_ylabel('Net energy flux \n (W m$^{-2}$)') 
        ax2.get_yaxis().set_label_coords(ylaboffs,0.5)
        
        lgd = ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':lfsize})
        ax2.yaxis.label.set_fontsize(lfsize)
        for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
            item.set_fontsize(axfsize)
        if pltitle: ax2.set_title(list_labels[nplot], fontsize = font_size)
        if nplot==1: 
            ax1 = ax2
        else:
            plt.draw()   # to generate tick labels
            list_tickpos1 = [blah.label.get_position()[-1] for blah in ax2.yaxis.get_major_ticks()]
            # filtering out all numbers that are not floats
            list_tickpos = []
            for val1 in list_tickpos1:
                if isinstance(val1, np.float64):
                    list_tickpos.append(val1)
            # calculating maxval that is as far as possible above highest tick without causing an additional tick
            tickmax = max(list_tickpos)
            tick1 = ax2.yaxis.get_major_ticks()[1].label.get_position()[-1]    
            tick2 = ax2.yaxis.get_major_ticks()[2].label.get_position()[-1]  
            maxval = tickmax + 0.9*abs(tick2 - tick1)
            ax2.set_ylim(top=maxval)   # to avoid overlapping tick marks with plot above


    if leaftemp:
        nplot = nplot+1
        if hspace1<0.1:
            try: 
                pylab.setp(ax2.get_xticklabels(), visible=False)   # Remove x-axis labels from previous plot
            except:
                pass        
        ax2 = plt.subplot(nrows,ncols,nplot, sharex = ax1)
        ax2.plot(xdata,results2['Tlmeas'][pos_vw]-Talist, linestyle='', marker='v', markersize=psize, label = '$T_l - T_a$ (obs.)')
        for i in srange(len(Tmods)):
            ax2.plot(xdata,results2[Tmods[i][0]][pos_vw]-Talist, color = 'gray', linewidth = lwidth, linestyle=Tmods[i][2], marker=None, label = '$T_l - T_a$ ' + Tmods[i][1])
        if alltemps:  
                for temp1 in alltemps: 
                    ax2.plot(xdata,results2[temp1[0]][pos_vw]+cdict[T0]-Talist, linestyle='', color='red', marker=temp1[2], markersize=psize, label = '$T_l - T_a$ ' + temp1[1])

        #ax2.tick_params(pad=axfsize/2)
        ax2.margins(y=0.1)
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':lfsize})
        ax2.set_ylabel('Leaf - air temp. \n (K)')  
        ax2.get_yaxis().set_label_coords(ylaboffs,0.5)
        ylim1 = ax2.get_ylim()
        if ylim1[0] > 0:
            ax2.set_ylim(bottom=0.)
        if ylim1[1] <0:
            ax2.set_ylim(top=0.)
        lgd = ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':lfsize})
        ax2.yaxis.label.set_fontsize(lfsize)
        for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
            item.set_fontsize(axfsize)
        if pltitle: ax2.set_title(list_labels[nplot], fontsize = font_size)
        if nplot==1: 
            ax1 = ax2
        else:
            plt.draw()   # to generate tick labels
            list_tickpos1 = [blah.label.get_position()[-1] for blah in ax2.yaxis.get_major_ticks()]
            # filtering out all numbers that are not floats
            list_tickpos = []
            for val1 in list_tickpos1:
                if isinstance(val1, np.float64):
                    list_tickpos.append(val1)
            # calculating maxval that is as far as possible above highest tick without causing an additional tick
            tickmax = max(list_tickpos)
            tick1 = ax2.yaxis.get_major_ticks()[1].label.get_position()[-1]    
            tick2 = ax2.yaxis.get_major_ticks()[2].label.get_position()[-1]  
            maxval = tickmax + 0.9*abs(tick2 - tick1)
            ax2.set_ylim(top=maxval)   # to avoid overlapping tick marks with plot above



    if gtwcomp:
        nplot = nplot+1
        if hspace1<0.1:
            try: 
                pylab.setp(ax2.get_xticklabels(), visible=False)   # Remove x-axis labels from previous plot
            except:
                pass
        ax2 = plt.subplot(nrows,ncols,nplot, sharex = ax1)
        ax2.plot(xdata,results2['g_twmeas'][pos_vw], linestyle='', color = 'green', marker='o', markersize=psize, label = '$g_{tw}$ (obs.)')
        ax2.plot(xdata,results2['g_tw'][pos_vw], color = 'green', linewidth = lwidth, linestyle='-', marker=None, label = '$g_{tw}$ (bulk)')

        ax2.margins(y=0.3, tight=False)
        ax2.set_ylim(bottom=0)
        #ax2.tick_params(pad=axfsize/2)
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':lfsize})
        ax2.set_ylabel('Leaf conductance \n (m s$^{-1}$)') 
        ax2.get_yaxis().set_label_coords(ylaboffs,0.5)
        ax2.yaxis.label.set_fontsize(lfsize)
        for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
            item.set_fontsize(axfsize)
        if pltitle: ax2.set_title(list_labels[nplot], fontsize = font_size)
        if nplot==1: 
            ax1 = ax2
        else:
            plt.draw()   # to generate tick labels
            list_tickpos1 = [blah.label.get_position()[-1] for blah in ax2.yaxis.get_major_ticks()]
            # filtering out all numbers that are not floats
            list_tickpos = []
            for val1 in list_tickpos1:
                if isinstance(val1, np.float64):
                    list_tickpos.append(val1)
            # calculating maxval that is as far as possible above highest tick without causing an additional tick
            tickmax = max(list_tickpos)
            tick1 = ax2.yaxis.get_major_ticks()[1].label.get_position()[-1]    
            tick2 = ax2.yaxis.get_major_ticks()[2].label.get_position()[-1]  
            maxval = tickmax + 0.9*abs(tick2 - tick1)
            ax2.set_ylim(top=maxval)   # to avoid overlapping tick marks with plot above

    if hccomp:
        nplot = nplot+1
        if hspace1<0.1:
            try: 
                pylab.setp(ax2.get_xticklabels(), visible=False)   # Remove x-axis labels from previous plot
            except:
                pass
        ax2 = plt.subplot(nrows,ncols,nplot, sharex = ax1)
        ax2.plot(xdata,results2['h_cmeas'][pos_vw], linestyle='', color = 'green', marker='o', markersize=psize, label = '$h_c$ (obs.)')
        ax2.plot(xdata,results2['h_c'][pos_vw], color = 'green', linewidth = lwidth, linestyle='-', marker=None, label = '$h_c$ (bulk)')

        ax2.margins(y=0.1, tight=False)
        ax2.set_ylim(bottom=0)
        #ax2.tick_params(pad=axfsize/2)
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':lfsize})
        ax2.set_ylabel('Heat transf. coeff. \n (W m$^{-2}$ K$^{-1}$)')  
        ax2.get_yaxis().set_label_coords(ylaboffs,0.5)
        ax2.yaxis.label.set_fontsize(lfsize)
        for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
            item.set_fontsize(axfsize)
        if pltitle: ax2.set_title(list_labels[nplot], fontsize = font_size)
        if nplot==1: 
            ax1 = ax2
        else:
            plt.draw()   # to generate tick labels
            list_tickpos1 = [blah.label.get_position()[-1] for blah in ax2.yaxis.get_major_ticks()]
            # filtering out all numbers that are not floats
            list_tickpos = []
            for val1 in list_tickpos1:
                if isinstance(val1, np.float64):
                    list_tickpos.append(val1)
            # calculating maxval that is as far as possible above highest tick without causing an additional tick
            tickmax = max(list_tickpos)
            tick1 = ax2.yaxis.get_major_ticks()[1].label.get_position()[-1]    
            tick2 = ax2.yaxis.get_major_ticks()[2].label.get_position()[-1]  
            maxval = tickmax + 0.9*abs(tick2 - tick1)
            ax2.set_ylim(top=maxval)   # to avoid overlapping tick marks with plot above


    if not pltitle: 
        ax1.xaxis.set_tick_params(labeltop='on', labelsize=axfsize)
        xticklabels1 = ax1.get_xticklabels()
        if max(xdata)>1000:
            ax1.set_xticklabels(xticklabels1, rotation=lrot)
        ax1.set_xlabel(axeslabel1, fontsize=lfsize)    
        ax1.xaxis.set_label_position('top') 

        

    # Setting x-axis label
    ax2.set_xlabel(axeslabel1)
    ax2.xaxis.label.set_fontsize(lfsize)  
    ax2.margins(x=0.05)
    if max(xdata)>1000:
        pylab.xticks(rotation=lrot)

    if fname:
        pylab.savefig(fname, bbox_extra_artists=(lgd,), bbox_inches='tight')
    pylab.savefig('foo.png', bbox_extra_artists=(lgd,), bbox_inches='tight')


# In[20]:

def fun_plottime_explicit(xdata1=None, ydata1=None, ylabel1=None, list_series = None, ynames1 = [''], yformats = ['b-'], fname = 'foo.png', dpi=50, facecolor='w',                           edgecolor='w',  orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches='tight', pad_inches=0.1,                           fontsize=12,legfontsize = 10, labfontsize = 10, frameon=None, dateformat = "%H:%M", ylim = (None, None), xlim = (None, None)):
    '''
    Plots time series of arrays or lists saved in ydata1, with time information in xdata1. y-axis label given in label1, series names in ynames1 and formatting information in yformats. Input can also be in the form of list_series = [(xdata1, ydata1, ynames1, yformats1), (xdata2, ydata2, ynames2, yformats2)].
    Other possible options:
        fname = 'foo.png', dpi=50, facecolor='w', edgecolor='w',  orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches='tight', pad_inches=0.1, fontsize=24,legfontsize = 20, labfontsize = 20, frameon=None, dateformat = "%H:%M", ylim = (None, None)
        xlim = (datetime.datetime(2014, 12, 11, 9, 57, 25), datetime.datetime(2014, 12, 11, 10, 57, 25))
        Examle:
            fun_plottime(xdata1,ydata1,'Amount')
        
    '''
    if list_series:
        xdata1 = [blah[0] for blah in list_series]
        ydata1 = [blah[1] for blah in list_series]
        ynames1 = [blah[2] for blah in list_series]
        yformats = [blah[3] for blah in list_series]
    
    pylab.close('all')
    #pylab.clf()
    pylab.xticks(rotation=70)
    for i in srange(len(ydata1)):
        xdata = xdata1[i].astype(datetime.datetime) 
        pylab.plot(xdata,ydata1[i], yformats[i], label = ynames1[i])
    pylab.gca().xaxis.set_major_formatter(mdates.DateFormatter(dateformat))   # Formats time axis (could add :%S to add seconds and .%f for decimals, see http://stackoverflow.com/questions/11107748/showing-milliseconds-in-matplotlib)
    pylab.tick_params(axis='x',which='minor')
    pylab.grid(color='b', which='both', axis = 'x')
    pylab.ylabel(ylabel1, fontsize = fontsize) # label the axes
    pylab.xlabel("Local time (hours:minutes)", fontsize = fontsize)
    if len(ydata1)>1:
        pylab.legend(loc=0, fontsize=legfontsize)
    pylab.tick_params(labelsize=labfontsize)
    pylab.xlim(xlim)
    pylab.ylim(ylim)
    plt.tight_layout()   # This makes automatically space for the vertical axes labels!
    pylab.savefig(fname, dpi=dpi, facecolor=facecolor, edgecolor=edgecolor,
        orientation=orientation, papertype=papertype, format=format,
        transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches,
        frameon=frameon)


# In[ ]:




# # Comparison of computed chamber vapour pressure with that measured by LI-6400XT

# In[21]:

str1 = ''
local_path = '2016_03_17_wet_leaf_Licor/'
# Create local_path if it does not exist yet
try:
    os.mkdir(path_data+local_path)
except:
    pass
nr_fname = local_path+'CR1000_net_radiometer_NetRadiometer'+str1+'.dat'
lc_fname = local_path+'CR1000_Leaf_chamber_LeaveChamber_PT100'+str1+'.dat'
csv_fname = local_path+'current_sage'+str1+'.csv'

li_fname = local_path+str1+'2016-03-17_stan_windtunnel'
print li_fname

pow_nameslist = ['Date', 'Time', 'Flow']
pow_formatslist = ['datetime64[us]','datetime64[us]', 'float']
sens_date_format = "%d.%m.%Y"
sens_time_format = "%H:%M:%S.%f"


li_data_all = fun_read_li(li_fname, prefix1 = ['2016-03-17_'])  
nr_data = fun_read_campbell(nr_fname)
lc_data = fun_read_campbell(lc_fname)
csv_data = fun_read_csv(csv_fname, lc_timepos = [])
csv_data['SensirionSLI'] = -csv_data['SensirionSLI'] # Making flow positive


# In[22]:

li_data = li_data_all[0]


# From the comments listed in the Licor file (see above), we can deduce the time intervals when the LI-6400 recorded the vapour concentration of the air entering the wind tunnel (whenever comments says "changed to humidifier". At all other times, the LI-6400XT was recording the vapour concentration in the air coming out of the wind tunnel.
# 
# ```
# ['11:42:42 changing to humidifier']
# ['11:45:40 increasinf tdew']
# ['11:47:12 back to chamber air']
# ['13:35:38 changing to humidifier']
# ['13:37:35 changing dew point']
# ['13:39:49 changing back to chamber air']
# ['15:27:42 changing to humidifier']
# ['15:31:45 back to chamber air and increasing Tdew']
# ['17:00:41 changing to humidifier']
# ['17:02:21 back to chamber and increasing dew point']
# ['18:08:42 changing to humidifier']
# ```

# In[23]:

# According to the comments listed in the Licor file (see above), the LI-6400 recorded the vapour concentration of the air entering the wind tunnel at the following times:

list_times_Pwin = [(datetime.datetime(2016, 3, 17, 11, 44, 0), datetime.timedelta(minutes=np.int(2))),                    (datetime.datetime(2016, 3, 17, 13, 36, 0), datetime.timedelta(minutes=np.int(2))),                    (datetime.datetime(2016, 3, 17, 15, 30, 0), datetime.timedelta(minutes=np.int(2))),                    (datetime.datetime(2016, 3, 17, 17, 1, 0), datetime.timedelta(minutes=np.int(1.5)))]
for (t1, dt1) in list_times_Pwin:
    posli1 = np.where(li_data['HHMMSS'] <= t1)[0][-1]
    posli2 = np.where(li_data['HHMMSS'] <= t1+dt1)[0][-1]
    print [posli1, posli2]


# In[24]:

# Solving saturation vapour pressure equation for temperature to calculate dewpoint from vapour pressure
eq_Td_Pwa = solve(eq_Pwl(P_wl = P_wa, T_l = T_d), T_d)[0]
show(eq_Td_Pwa)


# In[25]:

# Computing dew point and vapour pressure from H2OR in LI-6400XT and adding to li_data1
import numpy.lib.recfunctions as rfn
list_Td = []
list_Pwa = []
for i in srange(len(li_data)):
    Pwa = li_data['H2OR'][i]/1000*li_data['Press'][i]*1000.
    Td = eq_Td_Pwa.rhs()(P_wa = Pwa).subs(cdict)
    list_Td.append(Td)
    list_Pwa.append(Pwa)
li_data1 = rfn.rec_append_fields(li_data,['Td'], [list_Td], dtypes = ['float'])
li_data1 = rfn.rec_append_fields(li_data1,['Pwa'], [list_Pwa], dtypes = ['float'])
print li_data1['Pwa']


# In[26]:

# Adding vapour pressure from LI-6400XT to csv_data
import numpy.lib.recfunctions as rfn
list_Pwin = []
list_pos = []
for i in srange(len(csv_data)):
    (t1, dt1) = list_times_Pwin[i]
    posli1 = np.where(li_data1['HHMMSS'] >= t1)[0][0]
    posli2 = np.where(li_data1['HHMMSS'] <= t1+dt1)[0][-1]
    Pwin = np.average(li_data1['Pwa'][posli1:posli2])
    list_Pwin.append(Pwin)
    list_pos.append((posli1, posli2))
csv_data1 = rfn.rec_append_fields(csv_data,['Pwin'], [list_Pwin], dtypes = ['float'])
print csv_data1['Pwin']


# We will now compute the wind tunnel results using the vapour pressure of the incoming air measured by the LI-6400XT, and then compare the computed outgoing vapour pressure with that measured by the LI-6400XT:

# In[27]:

vdict = cdict.copy()
vdict[a_s] = 1
vdict[L_l] = 0.03
vdict[P_a] = 101325.
vdict[R_s] = 0
vdict[Re_c] = 3000.
vdict[g_sw] = 999
vdict[k_l] = 0.27  # between 0.268 and 0.573 according to Hays_1975_The_thermal.pdf
vdict[z_l] = 0.0005
ndict = {P_w_in: 'Pwin', T_l: 'T_leafTC', E_lmol: 'SensirionSLI'}
results_orig = fun_results(csv_data1, vdict1 = vdict, ndict1 = ndict, twosides=True, Pwinoffset=0)

#results1 = fun_results(csv_data1, ndict1 = {P_w_in: 'Pwin', T_l: 'T_leafTC', E_lmol: 'SensirionSLI', })
results1 = results_orig.copy()
list_vars = ['R_s', 'P_wa', 'T_a', 'v_w', 'g_sw']
for var1 in list_vars:
    print var1 + ' = (' + str(min(results1[var1])) + ', ' + str(max(results1[var1])) + ') ' + str(udict[eval(var1)]/eval(var1))
fun_plot_TN(results1, varname1 = 'P_wa', axeslabel1 = 'Vapour pressure (kPa)', xfac = 1/1000., fsize = 14,  alltemps = [('T_leafTC', 'Inside', 'o')], Emods = [('E_l', '(bulk)', '-'), ('2s_E_l', '(2s)', '--')], Hmods = [('H_l', '(bulk)', '-'), ('2s_H_l', '(2s)', '--')], Tmods = [('T_l', '(bulk)', '-'), ('2s_T_lu', '(upper side)', '--'), ('2s_T_ll', '(lower side)', ':')])


# In[28]:

xdata = [li_data1['HHMMSS'], results1['lc_time'], results1['lc_time']]
ydata = [li_data1['Pwa'], results1['P_w_in'], results1['P_wa']]
ylabel1 = 'Vapour pressure (Pa)'
ynames1 = ['Chamber (obs.)', 'Incoming (obs.)', 'Chamber (sim.)']
yformats = ['b-', 'ko', 'ro']
fname1 = path_figs + 'Pwa_tseries.pdf'
fun_plottime_explicit(xdata, ydata, ylabel1, ynames1 = ynames1, yformats = yformats, fname = fname1, dpi = 100)


# In the above plot, the blue line represents the continous measurements of the LI-6400XT, which were periodically switched between incoming and outgoing air. The blue dots represent the values for the incoming air, as extracted from the blue line, while the red dots represent the calculated values for the outgoing air. Correspondence between the red dots and the blue line suggests that the vapour pressure of the outgoing air was correctly simulated.

# # Experiments in the dark

# ## Only filter paper

# In[29]:

fname = 'exp_maxgs1.csv'
lc_data_orig = fun_read_csv(fname)
lc_data = lc_data_orig.copy()
lc_data['water_flow_avg'] = -lc_data_orig['water_flow_avg'] # Making flow positive
lc_data_orig = lc_data.copy()


# In[30]:

pos_Tdew = [0..7]
pos_vw = [8..13]
pos_vw_orig = pos_vw[:]
pos_Tdew_orig = pos_Tdew[:]


# In[31]:

vdict = cdict.copy()
vdict[a_s] = 1
vdict[L_l] = 0.03
vdict[P_a] = 101325.
vdict[R_s] = 0
vdict[Re_c] = 3000.
vdict[g_sw] = 999  # Non-limiting
vdict[k_l] = 0.27  # between 0.268 and 0.573 according to Hays_1975_The_thermal.pdf
vdict[z_l] = 0.0005
results_orig = fun_results(lc_data, vdict1 = vdict, twosides=True)


# <p><span style="color: #ff0000;">Note that the filter paper had dried out on the edges during the first part experiment. Only the last few data points of the Tdew variation were better, when the water reservoir was elevated. The wind variation should be fine.<br /></span></p>

# ### Varying wind speed

# In[32]:

results1 = results_orig[pos_vw_orig]
list_vars = ['R_s', 'P_wa', 'T_a', 'v_w', 'g_sw']
for var1 in list_vars:
    print var1 + ' = (' + str(min(results1[var1])) + ', ' + str(max(results1[var1])) + ') ' + str(udict[eval(var1)]/eval(var1))
fun_plot_TN(results1, Emods = [('E_l', '(mod.)', '-'), ('El_ball', '(Ball)', '--')], Hmods = [('H_l', '(mod.)', '-'), ('Hl_ball', '(Ball)', '--')], Tmods = [('T_l', '(mod.)', '-'), ('Tl_ball', '(Ball)', '--')], alltemps = False, fname = '/home/sschyman/Documents/papers/leaf-windtunnel/figures/maxgs_vw_Ball')


# ### Calibrated thermal conductivity, Fig. 9a in Schymanski et al. (2016, HESSD)

# In[33]:

vdict = cdict.copy()
vdict[a_s] = 1
vdict[L_l] = 0.03
vdict[P_a] = 101325.
vdict[R_s] = 0
vdict[Re_c] = 3000.
vdict[g_sw] = 999
vdict[k_l] = 0.1  # between 0.268 and 0.573 according to Hays_1975_The_thermal.pdf
vdict[z_l] = 0.0005

results_orig = fun_results(lc_data, vdict1 = vdict, twosides=True)

results1 = results_orig[pos_vw_orig]
list_vars = ['R_s', 'P_wa', 'T_a', 'v_w', 'g_sw']
for var1 in list_vars:
    print var1 + ' = (' + str(min(results1[var1])) + ', ' + str(max(results1[var1])) + ') ' + str(udict[eval(var1)]/eval(var1))
fname = path_figs + 'maxgs_vw_kl1_2s_all.eps'
fun_plot_TN(results1, varname1 = 'v_w', axeslabel1 = 'Wind speed (m s$^{-1}$)', fsize = 14,  alltemps = False,         Emods = [('E_l', '(bulk)', '-'), ('2s_E_l', '(2s)', '--')], Hmods = [('H_l', '(bulk)', '-'), ('2s_H_l', '(2s)', '--')],         Tmods = [('T_l', '(bulk)', '-'), ('2s_T_lu', '(upper side)', '--'), ('2s_T_ll', '(lower side)', ':')], dpi1=50, fname=fname)


# ## Black foil, maxgs, Fig. 9b in Schymanski et al. (2016, HESSD)

# In[34]:

fname = 'exp_maxgs_black_Tdew.csv'
lc_data_orig = fun_read_csv(fname, lc_timepos = [0, 15, 24, 25, 29, 30])
lc_data = lc_data_orig.copy()
lc_data['water_flow_avg'] = -lc_data_orig['water_flow_avg'] # Making flow positive
lc_data_orig = lc_data.copy()


# In[35]:

# Calibrated thermal conductivity
lc_data = lc_data_orig.copy()
vdict = cdict.copy()
vdict[a_s] = 1
vdict[L_l] = 0.03
vdict[P_a] = 101325.
vdict[R_s] = 0
vdict[Re_c] = 3000.
vdict[g_sw] = 999
vdict[k_l] = 0.1  # between 0.268 and 0.573 according to Hays_1975_The_thermal.pdf
vdict[z_l] = 0.0005

results_orig = fun_results(lc_data, vdict1 = vdict, twosides=True)

results1 = results_orig.copy()
list_vars = ['R_s', 'P_wa', 'T_a', 'v_w', 'g_sw']
for var1 in list_vars:
    print var1 + ' = (' + str(min(results1[var1])) + ', ' + str(max(results1[var1])) + ') ' + str(udict[eval(var1)]/eval(var1))
fname = path_figs + 'maxgs_Pwa_kl1_2s_all.eps'
fun_plot_TN(results1, varname1 = 'P_wa', axeslabel1 ='Vapour pressure (kPa)', xfac = 1/1000., fsize = 14,  alltemps = False,         Emods = [('E_l', '(bulk)', '-'), ('2s_E_l', '(2s)', '--')], Hmods = [('H_l', '(bulk)', '-'), ('2s_H_l', '(2s)', '--')],         Tmods = [('T_l', '(bulk)', '-'), ('2s_T_lu', '(upper side)', '--'), ('2s_T_ll', '(lower side)', ':')], dpi1=50, fname=fname)


# <p><span style="color: #ff0000;">If the effective leaf thermal conductivity of the artificial leaf is only 30% of that of real leaves, the upper side of the leaf could have similar temperatures as the TC readings, without significant change to the fluxes.<br /></span></p>

# ## Perforated foil 63_1, Figure A4a in Schymanski et al. (2016, HESSD)

# In[36]:

fname = 'exp_63_1_thin.csv'
lc_data_orig = fun_read_csv(fname)
lc_data = lc_data_orig.copy()
lc_data['water_flow_avg'] = -lc_data_orig['water_flow_avg'] # Making flow positive
lc_data_orig = lc_data.copy()

# Calibrated k_l
vdict = cdict.copy()
vdict[a_s] = 1
vdict[L_l] = 0.03
vdict[P_a] = 101325.
vdict[R_s] = 0
vdict[Re_c] = 3000.
vdict[k_l] = 0.03  # between 0.268 and 0.573 according to Hays_1975_The_thermal.pdf
vdict[z_l] = 0.00035  # Thin leaf
vdict[g_sw] = 0.05   # Tuned
results_orig = fun_results(lc_data, vdict1 = vdict, twosides=True)
results1 = results_orig.copy()
list_vars = ['R_s', 'P_wa', 'T_a', 'v_w', 'g_sw']
for var1 in list_vars:
    print var1 + ' = (' + str(min(results1[var1])) + ', ' + str(max(results1[var1])) + ') ' + str(udict[eval(var1)]/eval(var1))
fname = path_figs + '63_vw_kl03_2s_all.eps'
fun_plot_TN(results1, varname1 = 'v_w', axeslabel1 = 'Wind speed (m s$^{-1}$)', fsize = 14,  alltemps = False,         Emods = [('E_l', '(bulk)', '-'), ('2s_E_l', '(2s)', '--')], Hmods = [('H_l', '(bulk)', '-'), ('2s_H_l', '(2s)', '--')],         Tmods = [('T_l', '(bulk)', '-'), ('2s_T_lu', '(upper side)', '--'), ('2s_T_ll', '(lower side)', ':')], dpi1=50, fname=fname)


# In[37]:

fun_plot_TN(results1, fsize = 14,  alltemps = False, Emods = [('E_l', '(S-mod.)', '-'), ('El_ball', '(B-mod.)', '--')], Hmods = [('H_l', '(S-mod.)', '-'), ('Hl_ball', '(B-mod.)', '--')], Tmods = [('T_l', '(S-mod.)', '-'), ('Tl_ball', '(B-mod.)', '--')])


# ## Black foil, 35.4 holes/mm2, Figure A5a in Schymanski et al. (2016, HESSD)
# In the below experiments, we also measured chamber air temperature just after the leaf, in order to assess the uncertainty in $H_l$ due to leaks in this area, given that the measured outflow is only half of the measured inflow.

# In[38]:

fname = 'exp_35_4_Tdew.csv'
lc_data_orig = fun_read_csv(fname, lc_timepos = [0, 16, 25, 26, 30, 31])
lc_data = lc_data_orig.copy()
lc_data['water_flow_avg'] = -lc_data_orig['water_flow_avg'] # Making flow positive
lc_data_orig = lc_data.copy()
lc_data_35_Pwa = lc_data.copy()


# In[39]:

vdict = cdict.copy()
vdict[a_s] = 1
vdict[L_l] = 0.03
vdict[P_a] = 101325.
vdict[R_s] = 0
vdict[Re_c] = 3000.
vdict[k_l] = 0.1  # between 0.268 and 0.573 according to Hays_1975_The_thermal.pdf
vdict[z_l] = 0.0005
vdict[g_sw] = 0.035   # According to perforated_foils_LO: range 0.023-0.051
results_orig = fun_results(lc_data, vdict1 = vdict, twosides=True)
results1 = results_orig.copy()
results_35 = results_orig.copy()
list_vars = ['R_s', 'P_wa', 'T_a', 'v_w', 'g_sw']
for var1 in list_vars:
    print var1 + ' = (' + str(min(results1[var1])) + ', ' + str(max(results1[var1])) + ') ' + str(udict[eval(var1)]/eval(var1))
fname = path_figs + '35_Pwa_klcal_2s_all.eps'
fun_plot_TN(results1, varname1 = 'P_wa', axeslabel1 = 'Vapour pressure (kPa)', xfac = 1/1000., axfsize = 20,  alltemps = False,         Emods = [('E_l', '(bulk)', '-'), ('2s_E_l', '(2s)', '--')], Hmods = [('H_l', '(bulk)', '-'), ('2s_H_l', '(2s)', '--')],         Tmods = [('T_l', '(bulk)', '-'), ('2s_T_lu', '(upper side)', '--'), ('2s_T_ll', '(lower side)', ':')], fname=fname)


# ## 35.4 holes/mm2, varying wind speed,  Figure A4b in Schymanski et al. (2016, HESSD)
# <p>$g_{sv}$-range according to perforated_foils_LO:</p>
# <pre class="shrunk">min/max:
# gsp: (0.043695677, 0.070854492)
# gsw: (0.019652234, 0.027908623)
# gsw_r0: (0.027980447, 0.04199113)
# gsw_r0_S: (0.033073541, 0.051419344)
# gsw_AO: (0.015633401, 0.02144536)
# 
# BUMP:
# gsp: (0.030809233, 0.042032648)
# gsw: (0.014939534, 0.019531313)
# gsw_r0: (0.020606836, 0.027298812)
# gsw_r0_S: (0.02390163, 0.031910311)
# gsw_AO: (0.011958141, 0.015504949)</pre>

# In[40]:

fname = 'exp_35_4_vw.csv'
lc_data_orig = fun_read_csv(fname, lc_timepos = [])
lc_data = lc_data_orig.copy()
lc_data['water_flow_avg'] = -lc_data_orig['water_flow_avg'] # Making flow positive
lc_data_orig = lc_data.copy()


# In[41]:

vdict = cdict.copy()
vdict[a_s] = 1
vdict[L_l] = 0.03
vdict[P_a] = 101325.
vdict[R_s] = 0
vdict[Re_c] = 3000.
vdict[k_l] = 0.1  # between 0.268 and 0.573 according to Hays_1975_The_thermal.pdf
vdict[z_l] = 0.0005
#vdict[g_sw] = (0.028+0.051)/2   # According to perforated_foils_LO: range 0.027--0.042
vdict[g_sw] = 0.042
print vdict[g_sw]
results_orig = fun_results(lc_data, vdict, twosides=True)
results1 = results_orig.copy()
list_vars = ['R_s', 'P_wa', 'T_a', 'v_w', 'g_sw']
for var1 in list_vars:
    print var1 + ' = (' + str(min(results1[var1])) + ', ' + str(max(results1[var1])) + ') ' + str(udict[eval(var1)]/eval(var1))

fname = path_figs + '35_vw_klcal_2s_all.eps'
fun_plot_TN(results1, fsize = 14,  alltemps = False,         Emods = [('E_l', '(bulk)', '-'), ('2s_E_l', '(2s)', '--')], Hmods = [('H_l', '(bulk)', '-'), ('2s_H_l', '(2s)', '--')],         Tmods = [('T_l', '(bulk)', '-'), ('2s_T_lu', '(upper side)', '--'), ('2s_T_ll', '(lower side)', ':')], dpi1=50, fname=fname)


# ## Leaf 7_1,  Figure A5b in Schymanski et al. (2016, HESSD)
# 

# In[42]:

fname1 = 'new_tunnel_chamber_2Tin_leaf.csv'
fname = path_data + fname1
try:
    reader = csv.reader(open(fname, 'rb'), delimiter=',', quotechar='"')
except:
    copyfile(path_data_orig + fname1, fname)
    reader = csv.reader(open(fname, 'rb'), delimiter=',', quotechar='"')

nameslist = reader.next()
unitslist = reader.next()
ncols = len(nameslist)
print ncols
print nameslist
print unitslist
csvdata = []
for row in reader :
    row1 = np.array(row)
    # replacing empty fields by NaN
    row1[row1==''] = 'NaN'
    row = tuple(row1)
    csvdata.append(row)

formatslist = []
for i in srange(len(unitslist)):
    if unitslist[i] == '':
        formatslist.append('S100')
    else:
        formatslist.append('float')
data = np.array(csvdata,dtype = zip(nameslist,formatslist))
data_orig = data.copy()

tabledata = []
tabledata.append(list(nameslist))
tabledata.append(list(unitslist))
for i in srange(len(data)):
    line1 = data[i]
    tabledata.append(list(line1) )
#print tabledata
table(tabledata)


# In[43]:

data.dtype


# In[44]:

pos_vlow = [3,4,5,6,7,9]
pos_vhigh = [10,11]


# In[45]:

vdict = cdict.copy()
vdict[a_s] = 1
vdict[L_l] = 0.03
vdict[P_a] = 101325.
vdict[R_s] = 0
vdict[Re_c] = 3000.
vdict[k_l] = 0.05  # between 0.268 and 0.573 according to Hays_1975_The_thermal.pdf
vdict[z_l] = 0.0005
vdict[g_sw] = 0.007   # Range: 0.005 to 0.01
ndict = {R_s: '', T_d: 'Tdew humidifier', T_l: 'Tlin', Q_in: 'Fan power', S_a: 'Rn_above_leaf', S_b: 'Rn_below_leaf', S_s: 'Rn_beside_leaf', T_in: 'Incoming2 Temp_C(5)', F_in_v: 'Inflow rate', v_w: 'Wind speed', T_a: 'chamber air Temp_C(1) ', E_lmol: 'Sensirion'}

results_orig = fun_results(data, vdict1=vdict, ndict1=ndict, twosides=True)
results1 = results_orig[pos_vlow]
list_vars = ['R_s', 'P_wa', 'T_a', 'v_w', 'g_sw']
for var1 in list_vars:
    print var1 + ' = (' + str(min(results1[var1])) + ', ' + str(max(results1[var1])) + ') ' + str(udict[eval(var1)]/eval(var1))
fname = path_figs + '7_Pwa_klcal_2s_all.eps'
fun_plot_TN(results1, varname1 = 'P_wa', axeslabel1 = 'Vapour pressure (kPa)', xfac = 1/1000., fsize = 14,  alltemps = False,         Emods = [('E_l', '(bulk)', '-'), ('2s_E_l', '(2s)', '--')], Hmods = [('H_l', '(bulk)', '-'), ('2s_H_l', '(2s)', '--')],         Tmods = [('T_l', '(bulk)', '-'), ('2s_T_lu', '(upper side)', '--'), ('2s_T_ll', '(lower side)', ':')], dpi1=50, fname=fname)


# <p><span style="color: #ff0000;">Even if the observed leaf temperature is 0.7 K above the simulated, observed leaf conductance and $E_l$ are close to simulated, as sensitivity of $E_l$ to $g_{tw}$ is much greater than to $T_l$ at low $g_{sw}$ (see below).Only observed $h_c$ is above the simulated, as expected from the difference in $T_l$.</span></p>

# ### Comparison of 35 and 7 pores/mm2

# In[46]:

leglength=2
lfsize = 12
axfsize = 18
figsize1 = [6,5]
psize = 24
varname1 = 'P_wa'
axeslabel1 = 'Vapour pressure (Pa)'
Emods = [('E_l', '(mod., 35)', ':')]
Hmods = [('H_l', '(mod., 35)', ':')]
lwidth = 3

# 35 holes data
results2 = results_35.copy()

results2 = np.sort(results2, order = varname1)
pos_vw = srange(len(results2))
xdata = results2[varname1][pos_vw]
Talist = results2['T_a'][pos_vw]

P = list_plot(zip(xdata,results2['Elmeas'][pos_vw]), frame = True, axes = False, plotjoined=False, marker='s', size=psize, legend_label = '$E_l$ (obs., 35)')
for i in srange(len(Emods)):
    tup1 = Emods[i]
    if len(tup1)<4: 
        tup1 = tuple(list(tup1) + ['blue'])
    P += list_plot(zip(xdata,results2[tup1[0]][pos_vw]), color = tup1[3], plotjoined=True, thickness = lwidth, linestyle=tup1[2], marker=None, legend_label = '$E_l$ ' + tup1[1])

    P += list_plot(zip(xdata,results2['Hlmeas'][pos_vw]), marker='s', faceted = True, color = 'white', markeredgecolor = 'black', size=psize, plotjoined=False, legend_label = '$H_l$ (obs., 35)')
for i in srange(len(Hmods)):
    tup1 = Hmods[i]
    if len(tup1)<4: 
        tup1 = tuple(list(tup1) + ['black'])
    P += list_plot(zip(xdata,results2[tup1[0]][pos_vw]), color = tup1[3], plotjoined=True, thickness = lwidth, linestyle=tup1[2], marker=None, legend_label = '$H_l$ ' + tup1[1])


# Now adding 7 pores / mm2
#Emods = [('E_l', '(S-mod.)', '-'), ('El_ball', '(B-mod.)', '--')]
#Hmods = [('H_l', '(S-mod.)', '-'), ('Hl_ball', '(B-mod.)', '--')]
Emods = [('E_l', '(mod., 7)', '--')]
Hmods = [('H_l', '(mod., 7)', '--')]
lwidth = 2

# Sorting array along v_w
results2 = results1.copy()

results2 = np.sort(results2, order = varname1)
pos_vw = srange(len(results2))
xdata = results2[varname1][pos_vw]
Talist = results2['T_a'][pos_vw]

P += list_plot(zip(xdata,results2['Elmeas'][pos_vw]), frame = True, axes = False, plotjoined=False, marker='o', size=psize, legend_label = '$E_l$ (obs., 7)')
for i in srange(len(Emods)):
    tup1 = Emods[i]
    if len(tup1)<4: 
        tup1 = tuple(list(tup1) + ['blue'])
    P += list_plot(zip(xdata,results2[tup1[0]][pos_vw]), color = tup1[3], plotjoined=True, thickness = lwidth, linestyle=tup1[2], marker=None, legend_label = '$E_l$ ' + tup1[1])

    P += list_plot(zip(xdata,results2['Hlmeas'][pos_vw]), marker='o', faceted = True, color = 'white', markeredgecolor = 'black', size=psize, plotjoined=False, legend_label = '$H_l$ (obs., 7)')
for i in srange(len(Hmods)):
    tup1 = Hmods[i]
    if len(tup1)<4: 
        tup1 = tuple(list(tup1) + ['black'])
    P += list_plot(zip(xdata,results2[tup1[0]][pos_vw]), color = tup1[3], plotjoined=True, thickness = lwidth, linestyle=tup1[2], marker=None, legend_label = '$H_l$ ' + tup1[1])



P.axes_labels([axeslabel1, 'Energy flux from leaf (W m$^{-2}$)'])
#P.save(fontsize = fsize, axes_labels_size = axfsize, fig_tight = True, figsize=figsize1, aspect_ratio = 'automatic', legend_font_size = lfsize, legend_loc=(1.01,0), legend_handlelength=leglength, filename = '/home/sschyman/Documents/papers/windtunnel_PM/figures/35_7_Pwa_anal.png') 
P.show(fontsize = lfsize, fig_tight = True, figsize=figsize1, aspect_ratio = 'automatic', legend_font_size = lfsize, legend_handlelength=leglength, legend_loc=(1.01,0.))
#P.show(fontsize = lfsize, axes_labels_size = axfsize, fig_tight = True, figsize=figsize1, aspect_ratio = 'automatic', legend_font_size = lfsize, legend_handlelength=leglength, legend_loc=(1.01,0.))


# # Experiments with light

# In[47]:

eq_Rs_Rll.show()


# ## max_gs, 400 W/m2 irradiance,  Figure 10a in Schymanski et al. (2016, HESSD)

# In[48]:

fname = 'exp_maxgs_Rs400_vw.csv'
lc_data_orig = fun_read_csv(fname, lc_timepos = [0, 16, 25, 26, 30, 31])
lc_data = lc_data_orig.copy()
lc_data['water_flow_avg'] = -lc_data_orig['water_flow_avg'] # Making flow positive
lc_data_orig = lc_data.copy()


# In[49]:

vdict = cdict.copy()
vdict[a_s] = 1
vdict[L_l] = 0.03
vdict[P_a] = 101325.
vdict[Re_c] = 3000.
vdict[k_l] = 0.27  # between 0.268 and 0.573 according to Hays_1975_The_thermal.pdf
vdict[z_l] = 0.0005

vdict[g_sw] = 999
#vdict[R_d] = 350
ndict = {R_s: 'Rn_above_leaf'}
results_orig = fun_results(lc_data, vdict, ndict1=ndict, twosides=True)
results1 = results_orig.copy()
list_vars = ['R_s', 'P_wa', 'T_a', 'v_w', 'g_sw']
for var1 in list_vars:
    print var1 + ' = (' + str(min(results1[var1])) + ', ' + str(max(results1[var1])) + ') ' + str(udict[eval(var1)]/eval(var1))
fname = path_figs + 'maxgs_vw_Rs350_klcal_2s_all.eps'
fun_plot_TN(results1, fsize = 14,  alltemps = False, hccomp=False,        Emods = [('E_l', '(bulk)', '-'), ('2s_E_l', '(2s)', '--')], Hmods = [('H_l', '(bulk)', '-'), ('2s_H_l', '(2s)', '--')],         Tmods = [('T_l', '(bulk)', '-'), ('2s_T_lu', '(upper side)', '--'), ('2s_T_ll', '(lower side)', ':')], dpi1=50, fname=fname)   


# **Both observed $R_{nleaf}$ and modelled $R_s - R_{ll}$ are well above the measured $R_s$ of 380-380 W m$^{-2}$. This is because the leaf is colder than the surrounding air, despite the irradiance, and absorbs radiative heat. As before, the emitted longwave is lower in the measurements than in the calculations by a factor of 2:**

# In[50]:

print results1['R_ll']
print 2*results1['Rn_below_leaf']


# ## 7 holes/mm2, 500 W/m2 light,  Figure 10b in Schymanski et al. (2016, HESSD)
# <p> </p>
# <pre class="shrunk">min/max:
# gsp: (0.012254157, 0.014780151)
# gsw: (0.00519383, 0.0059423419)
# gsw_r0: (0.0074208532, 0.0086463923)
# gsw_r0_S: (0.0087905126, 0.010359369)
# gsw_AO: (0.0040818779, 0.0046200179)
# 
# BUMP:
# gsp: (0.0079768673, 0.0099597964)
# gsw: (0.0036078144, 0.0044655325)
# gsw_r0: (0.0050318656, 0.0062556397)
# gsw_r0_S: (0.0058744419, 0.0073222672)
# gsw_AO: (0.0029048966, 0.0035467658)</pre>
# <p>$g_{sv}$ in the range of 0.005-0.0073 (see perforated_foils_LO.sws, 7 holes/mm2)</p>

# In[51]:

fname = 'exp_7_2_Rs500_vw.csv'
lc_data_orig = fun_read_csv(fname, lc_timepos = [0, 16, 25, 26, 30, 31])
lc_data = lc_data_orig.copy()
lc_data['water_flow_avg'] = -lc_data_orig['water_flow_avg'] # Making flow positive
lc_data_orig = lc_data.copy()


# In[52]:

vdict = cdict.copy()
vdict[a_s] = 1
vdict[L_l] = 0.03
vdict[P_a] = 101325.
vdict[R_s] = 0
vdict[Re_c] = 3000.
vdict[k_l] = 0.3 # between 0.268 and 0.573 according to Hays_1975_The_thermal.pdf
vdict[z_l] = 0.0005
vdict[g_sw] = 0.008
#vdict1[R_d] = 550
ndict = {R_s: 'Rn_above_leaf'}
results_orig = fun_results(lc_data, vdict, ndict1=ndict, twosides=True)

results1 = results_orig.copy()
list_vars = ['R_s', 'P_wa', 'T_a', 'v_w', 'g_sw']
for var1 in list_vars:
    print var1 + ' = (' + str(min(results1[var1])) + ', ' + str(max(results1[var1])) + ') ' + str(udict[eval(var1)]/eval(var1))
fname = path_figs + '7_vw_Rs550_klcal_2s_all.eps'
fun_plot_TN(results1, fsize = 14,  alltemps = False, hccomp=False, Hobs=False,        Emods = [('E_l', '(bulk)', '-'), ('2s_E_l', '(2s)', '--')], Hmods = [('H_l', '(bulk)', '-'), ('2s_H_l', '(2s)', '--')],         Tmods = [('T_l', '(bulk)', '-'), ('2s_T_lu', '(upper side)', '--'), ('2s_T_ll', '(lower side)', ':')], dpi1=50, fname=fname)       


# <p><span style="color: #ff0000;">The two low outliers were the last point measured, so it could be that the filter paper started drying out towards the end.</span></p>
