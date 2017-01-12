
# coding: utf-8

# # Analysis of images of perforated foils and calculation of diffusive conductances
# Author: stan.schymanski@env.ethz.ch

# # Functions to read data files and compute conductances

# In[1]:

get_ipython().run_cell_magic(u'capture', u'storage', u"# The above redirects all output of the below commands to the variable 'storage' instead of displaying it.\n# It can be viewed by typing: 'storage()'\n\n# Setup of worksheet, incl. importing of packages and defining general custom functions\nload('temp/Worksheet_setup.sage')\n# From stomatal_cond_eqs\nload_session('temp/stomatal_cond_eqs.sobj')\nimport csv")


# In[2]:

path_data = 'data/perforated_foils/'


# In[3]:

# Adding foil thickness to dictionary of general constants, cdict:
cdict[d_p] = 25e-6   # Foil thickness was 25 um.
dict_print(cdict)


# ## Functions to read bump data and compute conductances

# In[4]:

fname =path_data+'63_7/63_7_1a_bump_a.csv'
fh = open(fname)
for i in srange(10):
    row = fh.readline()
    print row


# In[5]:

def stat_bump(fname, delimiter1=','):
    '''
    Returns statistics of holes in image, given name of csv file
    '''

    # Determining number of rows to skip
    with open(fname, 'U') as fh:
        data = []
        # Counting rows to be skipped
        skiplines = 0
        dataline = fh.readline()
        #print dataline
        while(len(dataline)<100):
            skiplines = skiplines+1
            dataline = fh.readline()
            if len(dataline) == 0:
                return 'empty line encountered'

    #print skiplines

    # Reading whole file
    with open(fname, 'rb') as fh:
        reader = csv.reader(fh, delimiter=delimiter1, quotechar='"')
        for i in srange(skiplines):
            reader.next()  


        # Creating list of column names
        nameslist = reader.next()
        #print nameslist[0]
        ncols = len(nameslist)
        nameslist1 = []
        name = ""
        for i in srange(ncols):
            if nameslist[i]!= "":
                name = nameslist[i]
            else:
                name = "unit_"+name
            nameslist1.append(name)
        formatslist = []
        for name in nameslist1:
            fmt = 'f8'
            if (name[0:4] == 'unit') or (name[0:4] == 'Comm'):
                fmt = 'S100'
            formatslist.append(fmt)  

        # Saving data in array
        data = []
        for row in reader:
            if len(row) != len(nameslist1):
                break
            data.append(tuple(row))  

    try:    data = np.array(data,dtype = zip(nameslist1,formatslist))
    except:
        print data[0]
        print nameslist1
        print formatslist
    
    # Removing artefacts (areas smaller than 1% of median)
    med_area = np.median(data['Surface area'])
    pos = np.where(data['Surface area'] < med_area/100.)
    data1 = np.delete(data,pos)
    data = data1.copy()

    porosity = sum(data['Area ratio'])
    area_av = np.average(data['C.S. area'])
    density_mm = (porosity/100.)/(area_av*(1e-6^2))*0.001^2   # number of pores per mm2

    area_std = np.std(data['C.S. area'])
    perimeter_av = np.average(data['Perimeter'])
    perimeter_std = np.std(data['Perimeter'])
    area_perimeter_av = np.average(data['C.S. area']/data['Perimeter'])
    area_perimeter_std = np.std(data['C.S. area']/data['Perimeter'])
    diameter_av = np.average(data['Circle equivalent dia'])
    diameter_std = np.std(data['Circle equivalent dia'])
    #print 'area based on diameter_av: ' + str(((diameter_av/2)^2*pi).n())
    #print 'C.S. area: ' + str(area_av)

    resdict = {}
    list_props = ['porosity', 'area_av', 'density_mm', 'area_std', 'perimeter_av', 'perimeter_std', 'area_perimeter_av', 'area_perimeter_std', 'diameter_av', 'diameter_std']
    for dummy in list_props:
        resdict[dummy] = eval(dummy)
    return resdict


# In[6]:

def gsv_bump(vdict_in, fname, delimiter1=','):
    # Standard values
    vdict = cdict.copy()
    vdict[T_a] = 300.     # Assumption
    vdict[P_a] = 101325
    vdict[k_dv] = 1e-3
    vdict[D_va] = eq_Dva.rhs().subs(vdict)
    vdict.update(vdict_in)  # updating all entries occurring in vdict_in
    example1 = stat_bump(fname, delimiter1 = delimiter1)
    
    try:
        vdict[F_p] = example1['porosity']/100.
        vdict[A_p] = example1['area_av']*(1e-6^2)
        vdict[n_p] = (F_p/A_p).subs(vdict)  # For regular pore distribution we can derive the pore density from porosity and pore area
        vdict[s_p] = eq_sp_np.rhs().subs(vdict)
        rp_keyence = example1['diameter_av']/2*1e-6  
        vdict[r_p] = eq_rp_Ap.rhs().subs(vdict).n()
        rperror = abs(rp_keyence - vdict[r_p])/rp_keyence
        if rperror > 0.05:
            print 'mismatch by ' + str(rperror*100)+ '%! rp_keyence: ' + str(rp_keyence) + '; r_p:' + str(vdict[r_p])
        
    except Exception, e:
        print e
        vdict[g_sw] = NaN
        return vdict

    gsp = eq_gsw_gswmol_iso.rhs().subs(g_swmol == 1/eq_rsp_rp.rhs()).subs(vdict).n()
    gsw = eq_gsw_gswmol_iso.rhs().subs(eq_gswmol.subs(eq_rend_rp, eq_rsp_rp, eq_rvs_B)).subs(vdict).n()
    gsw_r0 = eq_gsw_gswmol_iso.rhs().subs(eq_gswmol.subs(r_end==0, eq_rsp_rp, eq_rvs_B)).subs(vdict).n()
    #gsw_r0_S = eq_gsw_gswmol_iso.rhs().subs(eq_gswmol.subs(r_end==0, eq_rsp_rp, eq_rvs_S)).subs(vdict).n()

    # Now computing gsw assuming 50% bigger stomatal area
    vdict1 = vdict.copy()
    vdict1[A_p] = 1.5*vdict[A_p]
    vdict1[r_p] = eq_rp_Ap.rhs().subs(vdict1)
    gsp50 = eq_gsw_gswmol_iso.rhs().subs(g_swmol == 1/eq_rsp_rp.rhs()).subs(vdict1).n()
    gsw50 = eq_gsw_gswmol_iso.rhs().subs(eq_gswmol.subs(eq_rend_rp, eq_rsp_rp, eq_rvs_B)).subs(vdict1).n()
    gsw_r050 = eq_gsw_gswmol_iso.rhs().subs(eq_gswmol.subs(r_end==0, eq_rsp_rp, eq_rvs_B)).subs(vdict1).n()    
    
    #vdict[g_sw] = gsw
    
    resdict = {}
    list_props = ['fname', 'gsp', 'gsw', 'gsw_r0', 'gsp50', 'gsw50', 'gsw_r050']
    
    for dummy in list_props:
        resdict[dummy] = eval(dummy)
    for dummy in example1.keys():
        resdict[str(dummy)] = example1[dummy]
    resdict['vdict'] = vdict
    return resdict


# In[7]:

# Test
vdict1 = {}
vdict1[T_l] = 300
dict_print(gsv_bump(vdict1, path_data + '35_4/new_analysis/35_2c_matte_10x_bump.csv', delimiter1=';'))


# <h1>Results</h1>

# In[8]:

# Creating named array for storing data of bump analyses
pathname = path_data + '63_7/'
fname = '63_7_1a_bump_a'
result_dict = gsv_bump(vdict1, pathname + fname + '.csv')

# Creating named array for storing data
names = [str(name) for name in result_dict.keys()]
formats = ['f4' for i in srange(len(names))]
for i in srange(len(names)):
    if names[i] == 'fname':
        formats[i] = 'S200'
    if names[i] == 'vdict':
        formats[i] = 'object'
           
results_array_b = np.array([],dtype = zip(names, formats))

# Arrays for computations with d_p=0
results_array_b0 = np.array([],dtype = zip(names, formats))


# In[9]:

def fun_append(results_array, result_dict, print_names = ['density_mm', 'area_av', 'diameter_av', 'radius_area_av', 'gsp', 'gsp50', 'gsw', 'gsw50', 'gsw_r0', 'gsw_r050']):
    '''
    Append content of results_dict to results_array in the order of names in results_array.
    results_array.dtype.names must be equal to result_dict.keys()
    Example:
        sage: fun_append(results_array, result_dict)
    '''
    if not list(results_array.dtype.names) == result_dict.keys():
        print sorted(results_array.dtype.names)
        print sorted(result_dict.keys())
        return 'dictionary does not have the same entries as array'
    dtype1 = results_array.dtype
    try: 
        line1 = np.array(tuple([result_dict[name1] for name1 in dtype1.names]), dtype = dtype1)
        outarray = np.append(results_array, line1)
    except Exception, e:
        print e
        print 'dtype: ' + str(dtype1)
        print 'result_dict: ' + str(result_dict)
    for name1 in print_names:
        try:    print name1 + ': ' + str(outarray[name1][-1])
        except: pass
    return outarray


# In[10]:

# Test
pathname = path_data + '63_7/'
fname = '63_7_1a_bump_a'
vdict1 = cdict.copy()
result_dict = gsv_bump(vdict1, pathname + fname + '.csv')
print 'number of holes per mm2: ' + str(result_dict['vdict'][n_p]*0.001^2)
blah = fun_append(results_array_b, result_dict)


# In[11]:

def fun_printminmax(results_array, varname, pathname, fnameconds = [], list_files = None, print1 = True):
    '''
    Finds all entries in pathname that contain all elements of
    fnameconds and returns min and max of varname. Prepend 'vdict_' to varname
    if it is stored in vdict.
    '''

    poslist = []
    for i in srange(len(results_array)):
        if list_files: 
            if any(s in results_array['fname'][i] for s in list_files):
                poslist.append(i)            
        else:
            fname = results_array['fname'][i]
            if fname[0:len(pathname)] == pathname:
                fnameconds1 = [True]
                for cond in fnameconds:
                    fnameconds1.append(cond in fname)
                    #print fname
                    #print fnameconds1
                if False not in fnameconds1:
                    #print 'true!'
                    poslist.append(i)

    if varname[0:6] == 'vdict_':
        listval = [results_array['vdict'][pos1][var(varname[6:])] for pos1 in poslist]
        minval = min(listval)
        maxval = max(listval)      
    else:   
        minval = min(results_array[varname][poslist])
        maxval = max(results_array[varname][poslist])
    if print1 == True:
        print varname + ': ' + str((minval, maxval))
    return (minval, maxval)


# In[12]:

def fun_calc_all(pathname = 'data/',                  list_fnames_b = None,                  list_printnames = ['density_mm', 'area_av', 'diameter_av', 'radius_area_av', 'gsp', 'gsp50', 'gsw', 'gsw50','gsw_r0', 'gsw_r050'],                  append1 = True, delimiter1 = ",", vdict1 = False):
    '''
    Runs all calculations, appends to arrays and prints summaries. 
    Boundary conditions are taken from cdict, unless provided as vdict1.
    '''
    global results_array_b, results_array_b0
   
    out_results_array_b = results_array_b.copy()
    out_results_array_b0 = results_array_b0.copy()
    
    if list_fnames_b is None:  # finding all csv files with "bump" in their names
        list_files = os.listdir(pathname)
        list_fnames_b1 =[s for s in list_files if "csv" in s and "bump" in s]
        list_fnames_b = sorted(list_fnames_b1)
    print list_fnames_b

    for fname in list_fnames_b:
        if vdict1 is False:
            vdict2 = cdict.copy()
        else:
            vdict2 = vdict1.copy()
            
        if fname[-4:] != '.csv':
            fullpath = pathname + fname + '.csv'
        else:
            fullpath = pathname + fname
        result_dict = gsv_bump(vdict2, fullpath, delimiter1 = delimiter1)
        #print 'number of holes per mm2: ' + str(result_dict['vdict'][n_p]*0.001^2)
        if NaN not in result_dict.values():
            out_results_array_b = fun_append(out_results_array_b, result_dict, print_names = [])

        # Same for d_p = 0
        vdict2[d_p] = 1e-12
        result_dict = gsv_bump(vdict2, fullpath, delimiter1 = delimiter1)
        #print 'number of holes per mm2: ' + str(result_dict['vdict'][n_p]*0.001^2)
        if NaN not in result_dict.values():
            out_results_array_b0 = fun_append(out_results_array_b0, result_dict, print_names = [])    
            
    print 'number of holes per mm2: ' + str(result_dict['vdict'][n_p]*0.001^2)
    print 'min/max BUMP:'
    for varname in list_printnames:
        try: blah = fun_printminmax(out_results_array_b, varname, pathname, list_files = list_fnames_b);
        except: pass

    print ''
    print 'min/max BUMP at d_p = 0:'
    for varname in list_printnames[4:]:
        try: blah = fun_printminmax(out_results_array_b0, varname, pathname, list_files = list_fnames_b);
        except: pass

    if append1:
        results_array_b = out_results_array_b
        results_array_b0 = out_results_array_b0
        
pathname1 = path_data + '63_7/'
fun_calc_all(pathname = pathname1, append1 = False)


# <h2>7 holes per mm2</h2>

# In[13]:

pathname1 = path_data + '7/'
fun_calc_all(pathname = pathname1, delimiter1 = ",")


# In[14]:

pathname1 = path_data + '7/new_analysis/'
fun_calc_all(pathname = pathname1, delimiter1 = ";")


# <h2>35.4 holes/mm2</h2>

# In[15]:

pathname1 = path_data + '35_4/'
fun_calc_all(pathname = pathname1, delimiter1 = ",")


# In[16]:

pathname1 = path_data + '35_4/new_analysis/'
fun_calc_all(pathname = pathname1, delimiter1 = ";")


# <h2>63.7 holes per mm1</h2>

# In[17]:

pathname1 = path_data + '63_7/'
fun_calc_all(pathname = pathname1, delimiter1 = ",")


# In[18]:

pathname1 = path_data + '63_7/new_analysis/'
fun_calc_all(pathname = pathname1, delimiter1 = ";")


# # Summary

# In[19]:

results_array_b['fname']


# In[20]:

list_names = ['63_7', '35_4']
list_varnames1 = ['density_mm', 'area_av', 'diameter_av']
list_varnames2 = ['gsw_r0', 'gsw_r050']

results_array = results_array_b.copy()

for name1 in list_names:
    print ''
    print name1
    print 'orig. foils: '
    pathname = 'data/perforated_foils/' + name1 + '/' + name1
    for varname in list_varnames1:
        blah = fun_printminmax(results_array, varname, pathname)
    for varname in list_varnames2:
        blah = fun_printminmax(results_array, varname, pathname, print1 = False)
        print varname + '(d_p = ' + str(int(cdict[d_p]*1e6)) + ' um) = ' + str(blah)
        blah = fun_printminmax(results_array_b0, varname, pathname, print1 = False)
        print varname + '(d_p = 0 um) = ' + str(blah)
    
    print ''    
    print 'shiny 10x: '
    pathname = 'data/perforated_foils/' + name1 + '/new_analysis/'
    for varname in list_varnames1:
        blah = fun_printminmax(results_array, varname, pathname, fnameconds = ['10x', 'shiny'])
    for varname in list_varnames2:
        blah = fun_printminmax(results_array, varname, pathname, fnameconds = ['10x', 'shiny'], print1 = False)
        print varname + '(d_p = ' + str(int(cdict[d_p]*1e6)) + ' um) = ' + str(blah)

    print ''    
    print 'matte 10x: '
    pathname = 'data/perforated_foils/' + name1 + '/new_analysis/'
    for varname in list_varnames1:
        blah = fun_printminmax(results_array, varname, pathname, fnameconds = ['10x', 'matte'])
    for varname in list_varnames2:
        blah = fun_printminmax(results_array, varname, pathname, fnameconds = ['10x', 'matte'], print1 = False)
        print varname + '(d_p = ' + str(int(cdict[d_p]*1e6)) + ' um) = ' + str(blah)

    print ''    
    print 'shiny 20x: '
    pathname = 'data/perforated_foils/' + name1 + '/new_analysis/'
    for varname in list_varnames1:
        blah = fun_printminmax(results_array, varname, pathname, fnameconds = ['20x', 'shiny'])
    for varname in list_varnames2:
        blah = fun_printminmax(results_array, varname, pathname, fnameconds = ['20x', 'shiny'], print1 = False)
        print varname + '(d_p = ' + str(int(cdict[d_p]*1e6)) + ' um) = ' + str(blah)
    
    print ''
    print 'matte 20x: '
    pathname = 'data/perforated_foils/' + name1 + '/new_analysis/'
    for varname in list_varnames1:
        blah = fun_printminmax(results_array, varname, pathname, fnameconds = ['20x', 'matte'])
    for varname in list_varnames2:
        blah = fun_printminmax(results_array, varname, pathname, fnameconds = ['20x', 'matte'], print1 = False)
        print varname + '(d_p = ' + str(int(cdict[d_p]*1e6)) + ' um) = ' + str(blah)
        
name1 = '7'
print ''
print name1
print 'orig. foils: '
pathname = 'data/perforated_foils/' + name1 + '/' + name1
for varname in list_varnames1:
    blah = fun_printminmax(results_array, varname, pathname)
for varname in list_varnames2:
    blah = fun_printminmax(results_array, varname, pathname, print1 = False)
    print varname + '(d_p = ' + str(int(cdict[d_p]*1e6)) + ' um) = ' + str(blah)
    blah = fun_printminmax(results_array_b0, varname, pathname, print1 = False)
    print varname + '(d_p = 0 um) = ' + str(blah)

print ''    
print 'shiny: '
pathname = 'data/perforated_foils/' + name1 + '/new_analysis/'
for varname in list_varnames1:
    blah = fun_printminmax(results_array, varname, pathname, fnameconds = ['shiny'])
for varname in list_varnames2:
    blah = fun_printminmax(results_array, varname, pathname, fnameconds = ['shiny'], print1 = False)
    print varname + '(d_p = ' + str(int(cdict[d_p]*1e6)) + ' um) = ' + str(blah)

print ''    
print 'matte: '
pathname = 'data/perforated_foils/' + name1 + '/new_analysis/'
for varname in list_varnames1:
    blah = fun_printminmax(results_array, varname, pathname, fnameconds = ['matte'])
for varname in list_varnames2:
    blah = fun_printminmax(results_array, varname, pathname, fnameconds = ['matte'], print1 = False)
    print varname + '(d_p = ' + str(int(cdict[d_p]*1e6)) + ' um) = ' + str(blah)


# ## Summary for Table 1 in Schymanski et al. (2016, HESSD)

# In[21]:

list_names = ['63_7', '35_4', '7']
list_varnames1 = ['density_mm', 'area_av', 'diameter_av']
list_varnames2 = ['gsw_r0', 'gsw_r050']

array_summary = []
for name1 in list_names:
    print ''
    print name1
    pathname = 'data/perforated_foils/' + name1 + '/'
    for varname in list_varnames1:
        blah = fun_printminmax(results_array_b, varname, pathname)
    for varname in list_varnames2:
        blah = fun_printminmax(results_array_b, varname, pathname, print1 = False)
        print varname + '(d_p = ' + str(int(cdict[d_p]*1e6)) + ' um) = ' + str(blah)
        #blah = fun_printminmax(results_array_b0, varname, pathname, print1 = False)
        #print varname + '(d_p = 0 um) = ' + str(blah)


# In[ ]:



