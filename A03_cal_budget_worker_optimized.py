""" This file does the calculatoin, it saves the files and deletes a file from the do list.

This scripts loads ERA5 hourly atmospheric data and tries to close the AM budget.

This is pyhton 3 !!, use cyc_ph3 environment or similar
requires files to be in processed/active!

"""
import os
exec(open(os.environ['PYTHONSTARTUP']).read())
exec(open(STARTUP_2019_surfacemom).read())

import m_general_ph3 as M
import glob
import xarray as xr
import time
import shutil
import datetime
from m_tools import write_log as logger
#import scipy.ndimage.filters as smoother

#%matplotlib inline
print('python version sys.version')
print('xr modules loaded')
print(xr.__version__)


# import varibables
def init_from_input(arguments):
    if (len(arguments) <= 1) | ('-f' in set(arguments) ) :

        date_str='1980-01-01'
        #date='2016-02'
        slavenumber=250
        config_file = 'processing_config'

        print('use standard values')
        print(date_str)
        print(config_file)
    else:
        date_str    =  arguments[1]
        slavenumber =  arguments[2]
        config_file =  arguments[3]
        print('read vars from file: ' +str(arguments[1]) )

    return date_str, slavenumber, config_file

date_str, slavenumber, config_file = init_from_input(sys.argv)

#mconfig['server']='casper' # for testing!!  mconfig['server']

""" config_file is not used jet """
# load config file
processor_path=mconfig['paths']['analysis']+'AM_budget/automated_ch/'
conf_path=processor_path+ '/budget_cal_configs/'+mconfig['server']+ '/'
Pconf =MT.json_load( config_file, conf_path)

# load date.json and move to active folder
try:
    data_conf =MT.json_load(date_str, conf_path+'../processed/active/')
    #shutil.move(conf_path+'../downloaded/'+date_str+'.json', conf_path+'../active/'+date_str+'.json')
except:
    print("could not find file, break ")
    raise ValueError('could not find file')

ttotstart=time.time()
# %%
load_path_plev = Pconf['save_path_level'] #mconfig['paths']['downloads'] + 'ERA5/temp_data_AM_budget/plevels/' # replace
load_path_srf  = Pconf['save_path_surface'] #mconfig['paths']['downloads']  + 'ERA5/temp_data_AM_budget/surface/' # replace
#load_path_srf  = mconfig['paths']['downloads'] + 'ERA5/temp_data_AM_budget/surface/'
plot_path = Pconf['plot_path'] #mconfig['paths']['plot']      + 'AM_budget/budget_control/'     # replace
save_path = Pconf['processed_save_path_base'] #mconfig['paths']['processed'] + 'AM_budget/data/'               # replace
log_path = Pconf['log_path'] +'/budget_workers/'#mconfig['paths']['plot']      + 'AM_budget/logs/budget_workers/' # replace

#MT.mkdirs_r(plot_path)
#MT.mkdirs_r(log_path)
#MT.mkdirs_r(save_path)

log_name='worker_log_num'+str(slavenumber)+'_'+date_str
key_name=Pconf['key_name'] #'ERA5_AM_eq_'              # replace
plot_cont= False
seg_dict={ 'africa_europa':slice(-17.9, 60), 'asia':slice(60.1, 180), 'americas':slice(-180.1, -18)  }

# define chunk sizes:
time_chunk_size=3
lat_chunk_size=60

#flist =os.listdir(load_path)
date_str_short=date_str.replace('-', '')
D= xr.open_dataset(load_path_plev+'ERA5_subdaily_plev_.25deg_'+date_str_short+'.nc', chunks={'time':time_chunk_size, 'latitude':lat_chunk_size})#.isel(level=range(-12, 0))
print(D.chunks)
len(D.chunks)
# help(xr.open_dataset)
# D= xr.open_dataset(load_path_plev+'ERA5_subdaily_plev_.25deg_'+date_str_short+'.nc', chunks=20)#.isel(level=range(-12, 0))
D['level'] = D.level*100.0 # to get Pascal

Dsurface= xr.open_dataset(load_path_srf+'ERA5_subdaily_plev_.25deg_gph_'+date_str_short+'.nc', chunks={'time':time_chunk_size, 'latitude':lat_chunk_size})#.sel(latitude=slice(0, -90))
# define varibables
Phi      = D.z
data     = D.drop('z')
ps       = Dsurface.sp
fluxes   = Dsurface['iews']
gravity_drag = Dsurface['megwss'] # These are the mean eastward gravity wave stresses,
#eventhough the ERA5 documentation claims that this should have units of Pa = N m**-2,
# they likely miss a factor of

# """ fluxes have units of N s / m^2 and are accumulations since the last time step. To get to N/m^2 devide by 60 sec """
# fluxes = Dsurface['ewss']
#fluxes.data = fluxes.data#/ 3600.0
#fluxes.attrs['units'] ='N m**-2'

# %%
def rho_beta(levels , ps, check=False):
    """
    This function returns the Boer beta field, given the levels in G and the ps (hPa)
    uses xarray

    returns:
    xr.DataArray in the same dimension as ps levels
    """
    aa= (levels < ps)
    if check is True:
        print('ratio ' + str( aa.sum()/float(aa.size)) )
    return aa*1 * ps#/ ps.mean()

def rho_beta_heaviside(levels , ps, check=False):
    """
    This function returns the Boer beta field, given the levels in G and the ps (hPa)
    uses xarray

    returns:
    xr.DataArray in the same dimension as ps levels
    """
    aa= (levels < ps)
    if check is True:
        print('ratio ' + str( aa.sum()/float(aa.size)) )
    return aa*1 #/ ps.mean()

def ddphi_spherical_zm (dd, ps_zm, r_e, lat, time_chunk=None  ):
    """
    This function calculates the gradient in meridional direction in a spherical system
    It takes and returns xarray.DataArrays

    inputs:
    dd          data xarray.DataArray with (latitude, time, level) or (latitude, time), or combinations there off
    ps_zm       xr.DataArray, Surfare pressure in the dimensions acoording dd, no copying to addional dimensions needed.
                2nd dimension should be latitude, if more then 2 dims.
    r_e         earth radius used in the spherical gradient
    lat         np.array, latitude values in degree, same size as dd.latitude

    returns:
    xr.DataArray same dimensions as dd

    """
    import xarray as xr

    # ensure correct chunks
    rechunk_dic=dict()
    for k in dd.dims:
        rechunk_dic[k]= dd[k].size

    if time_chunk is not None:
        rechunk_dic['time']= time_chunk
    dd= dd.chunk(rechunk_dic)

    #plt.hist(np.diff(lat_radiens))

    lat_radiens =lat *np.pi/180.0
    cos_phi= np.cos(lat_radiens)

    if ps_zm is None:
        print('no ps weight lat gradient')
        ps_dummy = dd.isel(level=1)*0+1
        grad_matrix = ps_dummy* r_e *cos_phi**2 * dd
    else:

        print('ps weight lat gradient')

        rechunk_dic=dict()
        for k in ps_zm.dims:
            rechunk_dic[k]= uzm_vzm_rep[k].size

        if time_chunk is not None:
            rechunk_dic['time']= time_chunk

        ps_zm=ps_zm.chunk(rechunk_dic)

        grad_matrix =ps_zm* r_e *cos_phi**2 * dd

    if lat.size != grad_matrix.shape[1]:
        grad_matrix= grad_matrix.T

    if lat.size != grad_matrix.shape[1]:
        raise ValueError('the 2nd dimension it not the same size as the latitude. make sure the input arrays as the cooriantes like (time, latitude, level) or (time, latitude)')

    grad_matrix_dphi = - grad_matrix.differentiate('latitude', edge_order=2)/(4.0*lat_radiens.diff('latitude').mean())
    #grad_matrix_dphi_np =np.gradient(grad_matrix, lat_radiens , axis=1)

    # ensure same order of diemnsions when data is returned
    # only for non-xarray fileds
    # trans_list=list()
    # for k in list(dd.shape):
    #     for i in [i for i,x in enumerate(list(grad_matrix_dphi.shape)) if x == k]:
    #         trans_list.append(i)

    #print(np.shape(r_e**2 *cos_phi**2))
    #print(np.shape(ps_zm * r_e**2 *cos_phi**2))


    if ps_zm is None:
        factor = r_e**2 *cos_phi**2
    else:
        factor = ps_zm * r_e**2 *cos_phi**2

    # non xarray version
    #dd_return = xr.DataArray(data=np.transpose(grad_matrix_dphi, trans_list), dims=dd.dims, coords=dd.coords ) /factor

    # xarray version
    dd_return = grad_matrix_dphi/factor

    return dd_return

def ps_weight_timemean(field, ps):
    """
    This takes the surface pressure time mean of atmos_fields
    input:
    field       xr.DataArray or xr.Dataset
    ps          surface pressure field with the same dimensions are field, it does not need the verical coordinates

    return
    same structure are field but time averaged

    """

    return (field * ps).mean('time') /ps.mean('time')

def ddlambda_spherical(dd, ps, r_e, lon, lat ):
    """
    This function calculates the gradient in meridional direction in a spherical system
    It takes and returns xarray.DataArrays
    It wraps around the longitudes, to ensure continous gradients at the boundary


    inputs:
    dd          data xarray.DataArray with (latitude, time, level) or (latitude, time), or combinations there off
    ps          !!! this variable seam to be not needed !!! xr.DataArray, Surfare pressure in the dimensions acoording dd, no copying to addional dimensions needed.
                2nd dimension should be latitude, if more then 2 dims.
    r_e         earth radius used in the spherical gradient
    lon         np.array, latitude values in degree, same size as dd.longitude
    lat         np.array, latitude values in degree, same size as dd.latitude

    returns:
    xr.DataArray same dimensions as dd

    """
    import xarray as xr


    lon_radiens =lon * np.pi/180.0
    lat_radiens =lat * np.pi/180.0
    cos_phi= np.cos(lat_radiens)

    #dd.shape
    if lon.size != dd.shape[2]:
        dd= dd.T

    # if lon.size != dd.shape[2]:
    #     raise ValueError('the 3rd dimension is not the same size as the latitude. make sure the input arrays as the cooriantes like (time, latitude, level) or (time, latitude)')

    # wrap longs
    if dd.name is None:
        dd.name='temp_name_ddy'

    dlon=lon.diff('longitude').mean()

    da = dd.isel(longitude=slice(0, 2) )
    da['longitude'] = np.array([dd.longitude[-1].data + dlon, dd.longitude[-1].data + 2* dlon])

    de = dd.isel(longitude=[-2, -1] )
    de['longitude'] = np.array([dd.longitude[0].data - 2 * dlon, dd.longitude[0].data - 1* dlon])
    dd2 =xr.merge([da, dd, de]).to_array()


    if ps is None:
        sp_dx = dd2.differentiate('longitude', edge_order=2)/(4*lon_radiens.diff('longitude').mean())
        #sp_dx = np.gradient(dd, lon_radiens, axis=dd.shape.index(dd.longitude.size))
        print('np pressure weighted lon gradient')
    else:
        print('pressure weighted lon gradient')
        sp_dx = (dd2 * ps).differentiate('longitude', edge_order=2)/(4*lon_radiens.diff('longitude').mean())
        #sp_dx = np.gradient(dd * ps, lon_radiens, axis=dd.shape.index(dd.longitude.size))


    #dd2_new = dd_new.differentiate('longitude', edge_order=2)/(4*lon_radiens.diff('longitude').mean())

    sp_dx= sp_dx.isel( longitude= (sp_dx.longitude >= lon[0]) & (sp_dx.longitude <= lon[-1]) )


    # ensure same order of diemnsions when data is returned
    # this is only need if sp_dx is not a n xarray
    # trans_list=list()
    # for k in list(dd.shape):
    #     for i in [i for i,x in enumerate(list(sp_dx.shape)) if x == k]:
    #         trans_list.append(i)

    if ps is None:
        factor = xr.DataArray(r_e*np.cos(lat*np.pi/180.0),  dims='latitude', coords=[dd.coords['latitude']])
    else:
        factor = ps*xr.DataArray(r_e*np.cos(lat*np.pi/180.0),  dims='latitude', coords=[dd.coords['latitude']])

    # non xarray version
    #sp_dx_adjust = xr.DataArray(data=np.transpose(sp_dx, trans_list), dims=dd.dims, coords=dd.coords) / factor

    # xarray version
    sp_dx_adjust = sp_dx/factor

    return sp_dx_adjust

def vertical_integal(dset):
    g =9.81
    dset_int = dset.integrate('level')/ g
    for k in dset_int.keys():
        if 'level' not in dset[k].coords.keys():
            dset_int[k] = dset_int[k] *g
            print(k)
    return dset_int

def vertical_integal_Hbeta(dset, Hb):
    g =9.81
    dset_int = (dset* Hb).integrate('level')/ g
    for k in dset_int.keys():
        if 'level' not in dset[k].coords.keys():
            dset_int[k] = dset[k]
            print(k)
    return dset_int

def plot_continent_seperation(all_dict, Gbudget):
    F = M.figure_axis_xy(6, 8)
    plt.suptitle('budget closure for continental seperation' , y=1.025)
    plt.subplot(3,1, 1)
    #all_CD_rep.keys()
    key='F_gwave_zm'
    plt.title('1 hour exmpl | '+ key)
    plt.plot(lat, all_dict['africa_europa'][key].isel(time=1), 'r-', label='africa_europa')
    plt.plot(lat, all_dict['asia'][key].isel(time=1), 'g-', label='asia')
    plt.plot(lat, all_dict['americas'][key].isel(time=1), 'b-', label='americas')

    plt.plot(lat, (all_dict['africa_europa'][key]+ all_dict['asia'][key]+ all_dict['americas'][key]).isel(time=1), 'k+' , label='sum')

    plt.plot(lat,Gbudget[key].isel(time=1), '-k', label='budget')
    plt.xlim(-90, -70)

    plt.legend()
    plt.subplot(3,1, 2)

    key='torque_lev_zm'

    plt.title('1 day mean | '+ key)
    plt.plot(lat, all_dict['africa_europa'][key].mean('time'), 'r-')
    plt.plot(lat, all_dict['asia'][key].mean('time'), 'g-')
    plt.plot(lat, all_dict['americas'][key].mean('time'), 'b-')

    plt.plot(lat, (all_dict['africa_europa'][key]+ all_dict['asia'][key]+ all_dict['americas'][key]).mean('time'), 'k-+')

    plt.plot(lat,Gbudget[key].mean('time'), '-k')
    plt.xlim(-78, -0)
    plt.ylim(-0.7, .7)
    #plt.ylim(-1*1e5, 1*1e5)


    plt.subplot(3,1, 3)
    key='torque_srf_zm'
    plt.title('1 hour exmpl | '+ key)

    plt.plot(lat, all_dict['africa_europa'][key].isel(time=1), 'r-')
    plt.plot(lat, all_dict['asia'][key].isel(time=1), 'g-')
    plt.plot(lat, all_dict['americas'][key].isel(time=1), 'b-')

    plt.plot(lat, (all_dict['africa_europa'][key]+ all_dict['asia'][key]+ all_dict['americas'][key]).isel(time=1), 'k+')

    plt.plot(lat,Gbudget[key].isel(time=1), '-k')
    plt.xlim(-78, 0)

    return F

# %%

# %%
hist = 'Start Processing bata'

tstart=tstart_start=time.time()

BATA = rho_beta(data.level, ps, check=False)
BATA.name='bata'
BATA_zm=BATA.mean('longitude').compute()

BATA_01 = rho_beta_heaviside(data.level, ps, check=False)
BATA_01.name='bata_01'
BATA_01_zm=BATA_01.mean('longitude').compute()

#data['bata']=bata
repres=dict()
budget=dict()
BATAD=dict()

BATAD['BATA_zm']       = BATA_zm
BATAD['BATA_zm_01']    = BATA_01_zm

tend=time.time()
hist = logger(hist, 'bata time: '  + str(tend-tstart) )
tstart=time.time()
print('Start Processing')
# 1. zonal mean terms
hist = logger(hist, '1. zonal mean terms', verbose = False)
print('1. zonal mean terms')





# %% 2. a) mountain torque
ps_zm= ps.mean('longitude').compute()

r_e = float(mconfig['constants']['radius'])
omega = float(mconfig['constants']['omega'])
lon=BATA.longitude
lat=BATA.latitude
f = 2 * omega * np.sin(lat *np.pi/180)


#  2. a) 1. surface tourque
sp_dlambda = ddlambda_spherical(Dsurface.sp, None, r_e, lon, lat ).sel(variable= 'sp').drop('variable')
gph_sp_dlambda= (Dsurface.z * sp_dlambda / 9.81).compute()

# take zonal mean, at this point treat it as a surface variable
gph_sp_dlambda_zm_rep = (gph_sp_dlambda).mean('longitude')
gph_sp_dlambda_zm_rep.name='zonal mean surface mountain torque'
gph_sp_dlambda_zm_rep.attrs['units']='N m**-2'
gph_sp_dlambda_zm_budget = gph_sp_dlambda_zm_rep * ps_zm
gph_sp_dlambda_zm_budget.attrs['units']='N**2 m**-4'

#gph_sp_dlambda_zm_rep.plot()
#(gph_sp_dlambda_zm_rep.mean('time')).plot()

# %%  2. a) 2. gph level tourque
gph_bata_div                    = ddlambda_spherical(BATA_01 * Phi , None, r_e, lon, lat ).compute().sel(variable='temp_name_ddy').drop('variable')
gph_bata_div_zm_rep             = (gph_bata_div * BATA).mean('longitude') / BATA_zm
gph_bata_div_zm_conventional    = gph_bata_div.mean('longitude')

# apply devinition of representative mean
gph_bata_div_zm_rep = gph_bata_div_zm_rep.where(BATA_zm != 0, gph_bata_div_zm_conventional)
gph_bata_div_zm_budget = gph_bata_div_zm_rep * BATA_zm #(data * BATA).mean('longitude')

tend=time.time()
hist = logger(hist, 'defining mountain tourque part I: '  + str(tend-tstart) )
tstart=time.time()

# %% 2. b) continenal excurse
""" start of continental excurse """
#D.longitude
# %% split out tourque terms per continent seperatly and save them


Nlon=dict()
N_tot= Nlon['N_total'] =float(D.longitude.size)
all_CD_rep=dict()
all_CD_bud=dict()

for kk in seg_dict.keys():

    CD_storage_rep=dict()
    CD_storage_bud=dict()
    lon_seg=seg_dict[kk]
    print(kk)

    Dsurface_segment    = Dsurface.sel(longitude=lon_seg)
    #ps_zm_seg           = ps.sel(longitude=lon_seg).mean('longitude').compute()
    BATA_seg            = BATA.sel(longitude=lon_seg)
    #BATA_zm_seg         = BATA_seg.mean('longitude').compute()
    BATA_01_zm_seg      = BATA_01.sel(longitude=lon_seg).mean('longitude').compute()
    Nlon[kk]            = float(Dsurface_segment.longitude.size)
    N_weight            = Nlon[kk]/N_tot

    # a) surface tourque
    #sp_dlambda_seg = ddlambda_spherical(Dsurface_segment.sp, None, r_e, lon, lat )
    #gph_sp_dlambda_seg= Dsurface_segment.z * sp_dlambda_seg / 9.81
    #better copy from above and select
    gph_sp_dlambda_seg= gph_sp_dlambda.sel(longitude=lon_seg)

    # take zonal mean, at this point treat it as a surface variable
    gph_sp_dlambda_zm_rep_seg = (gph_sp_dlambda_seg).mean('longitude')
    gph_sp_dlambda_zm_rep_seg.name='zonal mean surface mountain torque'
    gph_sp_dlambda_zm_rep_seg.attrs['units']='N m**-2'
    gph_sp_dlambda_zm_budget_seg = gph_sp_dlambda_zm_rep_seg * ps_zm
    gph_sp_dlambda_zm_budget_seg.attrs['units']='N**2 m**-4'

    CD_storage_rep['torque_srf_zm']= gph_sp_dlambda_zm_rep_seg * N_weight
    CD_storage_bud['torque_srf_zm']= gph_sp_dlambda_zm_budget_seg * N_weight

    # b) GPh level
    a =gph_bata_div.sel(longitude=lon_seg)
    a_zm_rep             = (a * BATA_seg).mean('longitude') / BATA_zm
    a_zm_conventional    = a.mean('longitude')

    # apply devinition of representative mean
    a_zm_rep = a_zm_rep.where(BATA_zm != 0, a_zm_conventional)
    a_zm_budget = a_zm_rep * BATA_zm #(data * BATA).mean('longitude')

    CD_storage_rep['torque_lev_zm']= (a_zm_rep * N_weight)
    CD_storage_bud['torque_lev_zm']= (a_zm_budget * N_weight)


    # c) gravity wave drag
    F_gravity_zm_data_seg =(Dsurface_segment['megwss']* Dsurface_segment.sp).mean('longitude')/ ps_zm
    F_gravity_zm_rep_seg = xr.DataArray(data=F_gravity_zm_data_seg, name='Zonal mean zonal gravity wave stress', attrs= gravity_drag.attrs)
    F_gravity_zm_rep_seg.attrs['units']='N m**-2'
    F_gravity_zm_budget_seg = F_gravity_zm_rep_seg * ps_zm
    F_gravity_zm_budget_seg.attrs['units']='N**2 m**-4'

    CD_storage_rep['F_gwave_zm']= (F_gravity_zm_rep_seg  * N_weight)
    CD_storage_bud['F_gwave_zm']= (F_gravity_zm_budget_seg * N_weight)

    # d) save

    CD_storage_rep = xr.Dataset(CD_storage_rep)
    G_CD_int =vertical_integal_Hbeta(CD_storage_rep,BATA_01_zm_seg ).compute()
    G_CD_int.attrs['long_name'], G_CD_int.attrs['units'] = 'Continental Surface Drag as Representetive Mean', 'Pa'

    CD_storage_bud = xr.Dataset(CD_storage_bud)
    GB_CD_int =vertical_integal_Hbeta(CD_storage_bud, BATA_01_zm_seg).compute()
    GB_CD_int.attrs['long_name'], GB_CD_int.attrs['units'] = 'Continental Surface Drag as Budget Mean', 'Pa**2'


    save_path_local = save_path + '/drag_on_continents_zm/repres/'
    MT.mkdirs_r(save_path_local)
    G_CD_int.to_netcdf(save_path_local + key_name + 'repres_'+ kk+ '_'+ date_str + '.nc')
    save_path_local = save_path + '/drag_on_continents_zm/budget/'
    MT.mkdirs_r(save_path_local)
    GB_CD_int.to_netcdf(save_path_local + key_name + 'budget_'+ kk+ '_'+ date_str + '.nc')


    all_CD_rep[kk]=G_CD_int
    #all_CD_bud[kk]=GB_CD_int

if plot_cont is False:
    del all_CD_rep
    del all_CD_bud
    del GB_CD_int
    del G_CD_int

del Dsurface_segment
del BATA_seg
del BATA_01_zm_seg
del CD_storage_rep
del CD_storage_bud
del a_zm_rep
del a_zm_budget
del a_zm_conventional


tend=time.time()
hist = logger(hist, 'continental sepeparatio time: '  + str(tend-tstart) )
tstart=time.time()

""" end of continental excurse """


# %%  2. c) finish up global mean of surface vars. compute and store them

# 2. c)  1 surface tourque

repres['torque_srf_zm']= gph_sp_dlambda_zm_rep.compute()
budget['torque_srf_zm']= gph_sp_dlambda_zm_budget.compute()

del gph_sp_dlambda_zm_rep
del gph_sp_dlambda_zm_budget
del sp_dlambda
del gph_sp_dlambda

# 2. c) 2. gph level
repres['torque_lev_zm']= gph_bata_div_zm_rep.compute()
budget['torque_lev_zm']= gph_bata_div_zm_budget.compute()

del gph_bata_div
del gph_bata_div_zm_conventional
del gph_bata_div_zm_rep
del gph_bata_div_zm_budget

tend=time.time()
hist = logger(hist, 'store and compute global surface tourque part II : '  + str(tend-tstart) )
tstart=time.time()


# %% 3 . global surface drag terms
# %% a) surface var: 1. turbulent stress
F_srf_zm_data =(fluxes* ps).mean('longitude')/ ps_zm
F_srf_zm_rep = xr.DataArray(data=F_srf_zm_data, name='Zonal mean zonal Surface Stress', attrs= fluxes.attrs)
F_srf_zm_rep.attrs['units']='N m**-2'
F_srf_zm_budget = F_srf_zm_rep * ps_zm
F_srf_zm_budget.attrs['units']='N**2 m**-4'

repres['F_tur_zm']  = F_srf_zm_rep
budget['F_tur_zm']  = F_srf_zm_budget
BATAD['ps_zm']      = ps_zm


# %% b) surface var: 2. gravity drag
F_gravity_zm_data =(gravity_drag* ps).mean('longitude')/ ps_zm
F_gravity_zm_rep = xr.DataArray(data=F_gravity_zm_data, name='Zonal mean zonal gravity wave stress', attrs= gravity_drag.attrs)
F_gravity_zm_rep.attrs['units']='N m**-2'
F_gravity_zm_budget = F_gravity_zm_rep * ps_zm
F_gravity_zm_budget.attrs['units']='N**2 m**-4'

repres['F_gwave_zm']  = F_gravity_zm_rep
budget['F_gwave_zm']  = F_gravity_zm_budget


# %% 4 . representative average
data_zm_rep          = (data * BATA).mean('longitude') / BATA_zm
data_zm_conventional = data.mean('longitude')
#levmask =data.level < 0.1e5
#data_zm_conventional = data.sel(level=~levmask).mean('longitude')
#data.sel(level=levmask)*np.nan

# apply devinition of representative mean
data_zm_rep     = data_zm_rep.where(BATA_zm != 0, data_zm_conventional).compute()
data_zm_budget  = (data_zm_rep * BATA_zm).compute() #(data * BATA).mean('longitude')

# store in dict
repres['data_zm'] = data_zm_rep
budget['data_zm'] = data_zm_budget



tend=time.time()
hist = logger(hist, 'Surface stresses and representative means: '  + str(tend-tstart) )
tstart=time.time()

# %% eddy terms
# a) mean flow:
uzm_vzm_rep    =  data_zm_rep.u * data_zm_rep.v
uzm_vzm_budget = BATA_zm * uzm_vzm_rep

# b) eddies
data_p =data - data_zm_rep
# test if primes are 0 , see Boer eq. sec. 4b.
#(data_p * BATA).mean()

upvp_zm_budget = (data_p.u * data_p.v * BATA_zm).mean('longitude').compute()
upvp_zm_conventional = (data_p.u * data_p.v).mean('longitude')

# apply devinition of representative mean
upvp_zm_rep = (upvp_zm_budget / BATA_zm)
upvp_zm_rep =upvp_zm_rep.where(BATA_zm != 0 ,  upvp_zm_conventional).compute()

#upvp_zm_rep
# store in dict
repres['uzm_vzm']             = uzm_vzm_rep
repres['uprime_vprime_zm']   = upvp_zm_rep

budget['uzm_vzm']             = uzm_vzm_budget
budget['uprime_vprime_zm']  = upvp_zm_budget

tend=time.time()
hist = logger(hist, 'eddy terms compute: '  + str(tend-tstart) )
tstart=time.time()

# %%
# 2. zonal derivative
# define constants

repres['uzm_vzm_div']            = ddphi_spherical_zm(uzm_vzm_rep ,ps_zm, r_e, lat ).where(lat !=-90, 0).where(lat !=90, 0).compute()
repres['uprime_vprime_zm_div']   = ddphi_spherical_zm(upvp_zm_rep ,ps_zm, r_e, lat ).where(lat !=-90, 0).where(lat !=90, 0).compute()

budget['uzm_vzm_div']            = ddphi_spherical_zm(uzm_vzm_budget   ,ps_zm, r_e, lat ).where(lat !=-90, 0).where(lat !=90, 0).compute()
budget['uprime_vprime_zm_div']   = ddphi_spherical_zm(upvp_zm_budget ,ps_zm, r_e, lat ).where(lat !=-90, 0).where(lat !=90, 0).compute()

tend=time.time()
hist = logger(hist, 'phi gradients compute: '  + str(tend-tstart) )
tstart=time.time()

# 3. tendency term
repres['dudt']  = data_zm_rep.u.differentiate('time', edge_order=2, datetime_unit='s').compute()
budget['dudt']  = data_zm_budget.u.differentiate('time', edge_order=2, datetime_unit='s').compute()


del data_zm_rep
del data_zm_budget
del data

tend=time.time()
hist = logger(hist, 'tendency term compute: '  + str(tend-tstart) )
tstart=time.time()


# %% 4. also process gph
Phi_zm_rep          = (Phi * BATA).mean('longitude') / BATA_zm
repres['phi_zm'] = Phi_zm_rep.where(BATA_zm != 0, Phi.mean('longitude'))
repres['phi_zm'].attrs =Phi.attrs
budget['phi_zm']  = Phi_zm_rep * BATA_zm #(data * BATA).mean('longitude')
budget['phi_zm'].attrs =Phi.attrs


tend=time.time()
hist = logger(hist, 'Phi single var compute time: '  + str(tend-tstart) )
tstart=time.time()
# %% 4. merge data to xr.DataSets
print('Repack data and Cal data')
# a) representetive means
G = repres['dudt']
G.name , G.attrs['long_name'], G.attrs['units'] = 'dudt' , 'Zonal Mean Tendency', '(m s**-2)'
G =repres['dudt'].to_dataset()
G.attrs['long_name'], G.attrs['units'] = 'Terms for Representative Zonal Mean Momentum Budget', 'm s**-2 (fields var) or N m**-2 (surface var)'

key= 'uzm_vzm'
repres[key].name , repres[key].attrs['long_name'], repres[key].attrs['units'] = 'uzm_vzm' , 'Zonal Mean mean-momentum flux', '(m**2 s**-2)'
G['uzm_vzm']=repres[key]

key= 'uprime_vprime_zm'
repres[key].name , repres[key].attrs['long_name'], repres[key].attrs['units'] = 'uprime_vprime_zm' , 'Zonal Mean eddy-moemntum flux', '(m**2 s**-2)'
G['uprime_vprime_zm']=repres[key]

key= 'uzm_vzm_div'
repres[key].name , repres[key].attrs['long_name'], repres[key].attrs['units'] = 'uzm_vzm_div' , 'Zonal Mean mean-flux divergence', '(m s**-2)'
G['uzm_vzm_div']=repres[key]

key= 'uprime_vprime_zm_div'
repres[key].name , repres[key].attrs['long_name'], repres[key].attrs['units'] = 'uprime_vprime_zm_div' , 'Zonal Mean eddy-flux divergence', '(m s**-2)'
G['uprime_vprime_zm_div']=repres[key]

key= 'torque_lev_zm'
repres[key].name , repres[key].attrs['long_name'], repres[key].attrs['units'] = 'torque_lev_zm' , 'Zonal Mean Geopotential Height Torque', '(m s**-2)'
G['torque_lev_zm']=repres[key].T

key= 'torque_srf_zm'
repres[key].name , repres[key].attrs['long_name'], repres[key].attrs['units'] = 'torque_srf_zm' , 'Zonal Mean Surface Torque', '(m s**-2)'
G['torque_srf_zm']=repres[key]

key= 'F_tur_zm'
repres[key].name , repres[key].attrs['long_name'] = 'F_tur_zm' , 'Zonal Mean Turbulent Surface Stress'
G['F_tur_zm']=repres[key]

key= 'F_gwave_zm'
repres[key].name , repres[key].attrs['long_name'] = 'F_gwave_zm' , 'Zonal Mean Zonal Gravity Wave Stress'
G['F_gwave_zm']=repres[key]


key= 'data_zm'
tempv =repres[key].v * f
tempv.name , tempv.attrs['long_name'], tempv.attrs['units'] = 'Zonal Mean Advection of Planetary Momentum' , 'Zonal Mean Advection of Planetary Momentum', '(m s**-2)'
G['v_f_zm']=tempv

# save also zonal mean winds and GPH
G_others = repres['data_zm'].u
G_others.name , G_others.attrs['long_name'], G_others.attrs['units'] = 'u_repres' , 'Representative Zonal Mean Zonal Wind' , '(m s**-1)'
G_others = G_others.to_dataset()

key='v_repres'
G_others[key]=repres['data_zm'].v
G_others[key].name , G_others[key].attrs['long_name'], G_others[key].attrs['units'] = key , 'Representative Zonal Mean Meridional Wind' , '(m s**-1)'

key= 'phi_repres'
G_others[key]=repres['phi_zm'].compute()
G_others[key].name , G_others[key].attrs['long_name'], G_others[key].attrs['units'] = key , 'Representative Zonal Mean Geopotential Height' , '(m**2 s**-2)'


# b) budget means
GB = budget['dudt']
GB.name , GB.attrs['long_name'], GB.attrs['units'] = 'dudt' , 'Zonal Mean Tendency', '(Pa m * s**-2)'
GB =budget['dudt'].to_dataset()
GB.attrs['long_name'], GB.attrs['units'] = 'Terms for Zonal Mean Momentum Budget', '(Pa m * s**-2)'

key= 'uzm_vzm'
budget[key].name , budget[key].attrs['long_name'], budget[key].attrs['units'] = 'uzm_vzm' , 'Zonal Mean mean-flux', '(Pa m**2 s**-2)'
GB['uzm_vzm']=budget[key]

key= 'uprime_vprime_zm'
budget[key].name , budget[key].attrs['long_name'], budget[key].attrs['units'] = 'uprime_vprime_zm' , 'Zonal Mean eddy-flux', '(Pa m**2 s**-2)'
GB['uprime_vprime_zm']=budget[key]

key= 'uzm_vzm_div'
budget[key].name , budget[key].attrs['long_name'], budget[key].attrs['units'] = 'uzm_vzm_div' , 'Zonal Mean mean-flux divergence', '(Pa m s**-2)'
GB['uzm_vzm_div']=budget[key]

key= 'uprime_vprime_zm_div'
budget[key].name , budget[key].attrs['long_name'], budget[key].attrs['units'] = 'uprime_vprime_zm_div' , 'Zonal Mean Eddy-flux divergence', '(Pa m s**-2)'
GB['uprime_vprime_zm_div']=budget[key]

key= 'torque_lev_zm'
budget[key].name , budget[key].attrs['long_name'], budget[key].attrs['units'] = 'torque_lev_zm' , 'Zonal Mean Geopotential Height Torque', '(Pa m s**-2)'
GB['torque_lev_zm']=budget[key].T

key= 'torque_srf_zm'
budget[key].name , budget[key].attrs['long_name'], budget[key].attrs['units'] = 'torque_srf_zm' , 'Zonal Mean Surface Torque', '(Pa m s**-2)'
GB['torque_srf_zm']=budget[key]

key= 'F_tur_zm'
budget[key].name , budget[key].attrs['long_name'] = 'F_tur_zm' , 'Zonal Mean Turbulent Surface Stress'
GB['F_tur_zm']=budget[key]

key= 'F_gwave_zm'
budget[key].name , budget[key].attrs['long_name'] = 'F_gwave_zm' , 'Zonal Mean Zonal Gravity Wave Stress'
GB['F_gwave_zm']=budget[key]

key= 'data_zm'
tempv =budget[key].v * f
tempv.name , tempv.attrs['long_name'], tempv.attrs['units'] = 'Zonal Mean Advection of Planetary Momentum' , 'Zonal Mean Advection of Planetary Momentum', '(Pa m s**-2)'
GB['v_f_zm']=tempv


# save also zonal mean winds and GPH
key='u_budget'
G_others[key]=budget['data_zm'].u
G_others[key].name , G_others[key].attrs['long_name'], G_others[key].attrs['units'] = key , 'Budget Zonal Mean Zonal Wind' , '(Pa m s**-1)'

key='v_budget'
G_others[key]=budget['data_zm'].v
G_others[key].name , G_others[key].attrs['long_name'], G_others[key].attrs['units'] = key , 'Budget Zonal Mean Meridional Wind' , '(Pa m s**-1)'

key= 'phi_budget'
G_others[key]=budget['phi_zm'].compute()
G_others[key].name , G_others[key].attrs['long_name'], G_others[key].attrs['units'] = key , 'Budget Zonal Mean Geopotential Height' , '(Pa m**2 s**-2)'


# close original files

Dsurface.close()
D.close()


tend=time.time()
print('define and zm process time :' + str(tend- tstart))
hist = logger(hist, 'define and zm process time :'+ str(tend- tstart))
tstart=time.time()


# %% test difference between Bdget and representative mean
# for k in G.keys():
#     try:
#         M.figure_axis_xy()
#         ilo=-1
#         arep= G[k].isel(time=2, level=ilo)
#         abud= GB[k].isel(time=2, level=ilo)
#
#         abata =BATA_zm.isel(time=2, level=ilo)
#         plt.plot(arep.latitude, arep*abata, '-+b', label='rep')
#         plt.plot(arep.latitude, abud, '-r', label='budget')
#         plt.title(k)
#         plt.legend()
#     except:
#         pass

# %% 5. Vertical integrals
print('5. Vertical integrals')

G_int =vertical_integal_Hbeta(G, BATA_01_zm) # bata not computed jet.
G_int.attrs['long_name'], G_int.attrs['units'] = 'Momentum Budget in the Representetive Mean', 'Pa'

GB_int =vertical_integal_Hbeta(GB, BATA_01_zm)
GB_int.attrs['long_name'], GB_int.attrs['units'] = 'Momentum Budget in the Budget Mean', 'Pa**2'

level_vars = list(G.keys())
flux_vars=['uprime_vprime_zm','uzm_vzm' ]

#level_vars.remove('F_tur_zm')
#level_vars.remove('F_gwave_zm')
for k in flux_vars:
    level_vars.remove(k)

for k in level_vars:
    G_int[k].attrs = G[k].attrs
    G_int[k].attrs['units'] = 'Pa'

    GB_int[k].attrs = GB[k].attrs
    GB_int[k].attrs['units'] = 'Pa**2'


for k in flux_vars:
    G_int[k].attrs = G[k].attrs
    G_int[k].attrs['units'] = 'Pa m'

    GB_int[k].attrs = GB[k].attrs
    GB_int[k].attrs['units'] = 'Pa**2 m'


# %% 5.b optional plotting of the continenal seperation
if plot_cont:
    F =plot_continent_seperation(all_CD_rep, G_int)
    F.save(name='exmpl_repres_'+ date_str, path=plot_path+'contitent_separ/')

    del all_CD_rep
    del all_CD_bud
    del GB_CD_int
    del G_CD_int


# %% 5.b save zonal mean data
date_str

save_path_local = save_path + '/instantanious_zm/repres/'
os.makedirs(save_path_local, exist_ok = True)
G.to_netcdf(save_path_local + key_name +'repres_zm_'+  date_str + '.nc')

save_path_local = save_path + '/instantanious_zm/budget/'
os.makedirs(save_path_local, exist_ok = True)
GB.to_netcdf(save_path_local + key_name +'budget_zm_'+  date_str + '.nc')

save_path_local = save_path + '/instantanious_eulerian_zm/'
os.makedirs(save_path_local, exist_ok = True)
G_others.to_netcdf(save_path_local + key_name +'zm_others_'+  date_str + '.nc')


G_bata = xr.Dataset(BATAD)
key='BATA_zm'
G_bata[key].name , G_bata[key].attrs['long_name'], G_bata[key].attrs['units'] = 'BATA_zm' , 'Zonal mean rho_beta', 'Pa'

key='BATA_zm_01'
G_bata[key].name , G_bata[key].attrs['long_name'], G_bata[key].attrs['units'] = 'BATA_zm_01' , 'Zonal mean H_beta', 'binary'

key='ps_zm'
G_bata[key].name , G_bata[key].attrs['long_name'], G_bata[key].attrs['units'] = 'ps_zm' , 'Zonal mean surface pressure', 'Pa'

save_path_local = save_path + '/instantanious_bata/'
os.makedirs(save_path_local, exist_ok = True)
G_bata.to_netcdf(save_path_local + key_name + 'bata_'+  date_str + '.nc')

# %%

#G_int['LHS'] = G_int['dudt'] - G_int['v_f_zm'] + G_int['uprime_vprime_zm_div'] + G_int['uzm_vzm_div']
#G_int['LHS'].name , G_int['LHS'].attrs['long_name'], G_int['LHS'].attrs['units'] = 'LHS' , 'Left Hand Side of the Butget', 'Pa'
#GB_int['LHS'] = GB_int['dudt'] - GB_int['v_f_zm'] + GB_int['uprime_vprime_zm_div'] + GB_int['uzm_vzm_div']
#GB_int['LHS'].name , GB_int['LHS'].attrs['long_name'], GB_int['LHS'].attrs['units'] = 'LHS' , 'Left Hand Side of the Butget', 'Pa**2'

#G_int['torque_lev_zm'].plot()
#G_int['torque_srf_zm'].plot()
#(G_int['torque_lev_zm'] - G_int['torque_srf_zm']).plot()

G_int['LHS'] = G_int['dudt'] - G_int['v_f_zm'] + G_int['uprime_vprime_zm_div'] + G_int['uzm_vzm_div'] + G_int['torque_lev_zm'] - G_int['torque_srf_zm'] - G_int['F_gwave_zm']
G_int['LHS'].name , G_int['LHS'].attrs['long_name'], G_int['LHS'].attrs['units'] = 'LHS' , 'Left Hand Side of the Butget', 'Pa'

GB_int['LHS'] = GB_int['dudt'] - GB_int['v_f_zm'] + GB_int['uprime_vprime_zm_div'] + GB_int['uzm_vzm_div'] + GB_int['torque_lev_zm'] - GB_int['torque_srf_zm'] - G_int['F_gwave_zm']
GB_int['LHS'].name , GB_int['LHS'].attrs['long_name'], GB_int['LHS'].attrs['units'] = 'LHS' , 'Left Hand Side of the Butget', 'Pa**2'


# %% delete old variables
del BATAD, repres, budget
del G, GB
del G_others, G_bata
del tempv
del uzm_vzm_rep, upvp_zm_rep, uzm_vzm_budget , upvp_zm_budget


# %%

# plotting
# G_int_nobeta= vertical_integal(G)
#
# for k in G.keys():
#     M.figure_axis_xy(15, 3)
#     plt.subplot(1,3, 1)
#     G_int_nobeta[k].plot()
#     plt.subplot(1,3, 2)
#
#     G_int[k].plot()
#     plt.subplot(1,3, 3)
#     (G_int_nobeta[k] - G_int[k]).plot()
#
# # plot all
# for k in G_int.keys():
#     M.figure_axis_xy()
#     G_int[k].plot()

# %% 6. time average
print('6. time average')
G_int_tmean =ps_weight_timemean(G_int, ps_zm )
G_int_tmean.attrs['long_name'], G_int_tmean.attrs['units'] = 'Daily Mean Momentum Budget of Representetive Means', 'Pa'
G_int_tmean.coords['time'] = G_int.time[12].data.astype('M8[h]') #G_int.time.mean().data.astype('M8[D]') # set to mid-da time
G_int_tmean =G_int_tmean.expand_dims('time')

GB_int_tmean =ps_weight_timemean(GB_int, ps_zm )
GB_int_tmean.attrs['long_name'], GB_int_tmean.attrs['units'] = 'Daily Mean Momentum Budget of Budget Means', 'Pa**2'
GB_int_tmean.coords['time'] = GB_int.time[12].data.astype('M8[h]')# set to mid-da time
GB_int_tmean =GB_int_tmean.expand_dims('time')

tend=time.time()
print('vert int and budget time :' + str(tend- tstart))
tstart=time.time()


# %% 6. b) save budgets
save_path_local = save_path + '/instantanious_vert_int/repres/'
os.makedirs(save_path_local, exist_ok = True)
G_int.to_netcdf(save_path_local + key_name +'repres_int_'+  date_str + '.nc')

save_path_local = save_path + '/instantanious_vert_int/budget/'
os.makedirs(save_path_local, exist_ok = True)
GB_int.to_netcdf(save_path_local + key_name +'budget_int_'+  date_str + '.nc')

save_path_local = save_path + '/daily/repres/'
os.makedirs(save_path_local, exist_ok = True)
G_int_tmean['time']=G_int_tmean.time.data.astype('M8[D]')
G_int_tmean.to_netcdf(save_path_local + key_name +'repres_tmean_int_'+  date_str + '.nc')

save_path_local = save_path + '/daily/budget/'
os.makedirs(save_path_local, exist_ok = True)
GB_int_tmean['time']=GB_int_tmean.time.data.astype('M8[D]')
GB_int_tmean.to_netcdf(save_path_local + key_name +'budget_tmean_int_'+  date_str + '.nc')


tend=time.time()
print('save time :' + str(tend- tstart))
hist = logger(hist, 'save time :'+ str(tend- tstart))
tstart=time.time()
# %%
# M.figure_axis_xy()
# for k in level_vars:
#     GB_int_tmean[k].plot(label=k)
#
# plt.legend()

# simple test
# %%


# %%  vertical integrals
#LHS = G_int['dudt'] - G_int['v_f_zm'] + G_int['uprime_vprime_zm_div'] + G_int['uzm_vzm_div']
F = M.figure_axis_xy()
plt.suptitle('Budget Mean | Veritcal integral | weighted time mean over 1 Day', y= 1.03)
plt.subplot(2,1,1)
pdata= GB_int_tmean.isel(time=0)
xlims=(pdata['LHS'].latitude.min().data , pdata['LHS'].latitude.max().data)

plt.plot(lat, pdata['LHS'],c='k',  label=' LHS')

plt.plot(lat, pdata['dudt'], label=' dudt')
plt.plot(lat, - pdata['v_f_zm'], label=' -f v')
plt.plot(lat, + pdata['uzm_vzm_div'], label=' div mean')
plt.plot(lat, + pdata['uprime_vprime_zm_div'], label=' div eddy')

plt.plot(lat, + pdata['torque_lev_zm'] - pdata['torque_srf_zm'], label='Torque Lev - Srf')

plt.plot(lat, -1*pdata['F_tur_zm'],'k--', linewidth= 2,   label=' F')
plt.legend(ncol=2)
plt.ylabel(pdata.units)

plt.title('LHS')
plt.grid()
plt.xlim(xlims)

plt.subplot(2,1,2)
plt.title('1 1hour time step')
tt=10
pdata= GB_int
plt.plot(lat, pdata['LHS'].isel(time=tt),c='k',  label=' LHS')

plt.plot(lat, pdata['dudt'].isel(time=tt), label=' dudt')
plt.plot(lat, - pdata['v_f_zm'].isel(time=tt), label=' -f v')
plt.plot(lat, + pdata['uzm_vzm_div'].isel(time=tt), label=' div mean')
plt.plot(lat, + pdata['uprime_vprime_zm_div'].isel(time=tt), label=' div eddy')

plt.plot(lat, + pdata['torque_lev_zm'].isel(time=tt) - pdata['torque_srf_zm'].isel(time=tt), label='Torque Lev - Srf')

plt.plot(lat, -1*pdata['F_tur_zm'].isel(time=tt),'k--', linewidth= 2,   label=' F')
plt.ylabel(pdata.units)

plt.title('F')
plt.grid()
plt.xlim(xlims)

#date_str =str(pdata.time[0].data.astype('M8[D]'))
F.save(name='exmpl_budget_dmean_'+ date_str, path=plot_path)

# %%
lat_radiens =lat *np.pi/180.0
cos_phi= np.cos(lat_radiens)

F = M.figure_axis_xy(6.5, 11)
plt.suptitle('Veritcal integral', y= 1.03)
plt.subplot(5,1,1)

plt.title('Representative Mean | Vertical Integral | weighed 1 Day Time mean')

pdata= G_int_tmean.isel(time=0)
plt.plot(lat, pdata['LHS'],c='k',  label=' LHS')

plt.plot(lat, pdata['dudt'], label=' dudt')
plt.plot(lat, - pdata['v_f_zm'], label=' -f v')
plt.plot(lat, + pdata['uzm_vzm_div'], label=' div mean')
plt.plot(lat, + pdata['uprime_vprime_zm_div'], label=' div eddy')

plt.plot(lat, + pdata['torque_lev_zm'] - pdata['torque_srf_zm'], label='Torque Lev - Srf')

plt.plot(lat, -1*pdata['F_tur_zm'],'k--', linewidth= 2,   label=' F')

plt.ylabel(pdata.units)
plt.legend(ncol=2)#loc='upper right',

plt.xlim(xlims)
plt.grid()


plt.subplot(5,1,2)

plt.title('Grouped | Vertical Integral | weighed 1 Day Time mean')

pdata= G_int_tmean.isel(time=0)
plt.plot(lat, pdata['LHS'],c='k',  label=' LHS')

plt.plot(lat, pdata['dudt'] + pdata['uprime_vprime_zm_div'] + pdata['uzm_vzm_div'], label='dudt - eddy div ')
#plt.plot(lat, + pdata['uzm_vzm_div'], label=' div mean')
#plt.plot(lat, , label=' div eddy')

plt.plot(lat, - pdata['v_f_zm'] + pdata['torque_lev_zm'] - pdata['torque_srf_zm'], label='- fv + Torque Lev - Srf')

plt.plot(lat, -1*pdata['F_tur_zm'],'k--', linewidth= 2,   label=' F')

plt.ylabel(pdata.units)
plt.legend(ncol=2)#loc='upper right',

plt.xlim(xlims)
plt.grid()


plt.subplot(5,1,3)

plt.title('single 1 hour time step ')

pdata= G_int
tt=1
plt.plot(lat, pdata['LHS'].isel(time=tt),c='k',  label=' LHS')

plt.plot(lat, pdata['dudt'].isel(time=tt), label=' dudt')
plt.plot(lat, - pdata['v_f_zm'].isel(time=tt), label=' -f v')
plt.plot(lat, + pdata['uzm_vzm_div'].isel(time=tt), label=' div mean')
plt.plot(lat, + pdata['uprime_vprime_zm_div'].isel(time=tt), label=' div eddy')

plt.plot(lat, -1*pdata['F_tur_zm'].isel(time=tt),'k--', linewidth= 2,   label=' F')

plt.plot(lat, + pdata['torque_lev_zm'].isel(time=tt) - pdata['torque_srf_zm'].isel(time=tt), label='Torque (Lev - Srf)')

plt.ylabel(pdata.units)
plt.xlim(xlims)
plt.grid()

plt.subplot(5,1,4)

plt.title('Balance for 1 day mean, expressed as AM/r ')
# plt.plot(lat, pdata['LHS'].isel(time=tt)* cos_phi,c='k',  label=' LHS')
# plt.plot(lat, - pdata['F_tur_zm'].isel(time=tt)* cos_phi ,'k--', linewidth= 2,   label='F')

plt.plot(lat, pdata['LHS'].mean('time')* cos_phi,c='k',  label='LHS')
plt.plot(lat, - pdata['F_tur_zm'].mean('time')* cos_phi ,'k--', linewidth= 2,   label='F')
plt.plot(lat, (pdata['LHS'] + pdata['F_tur_zm'] ).mean('time')* cos_phi,'r-', linewidth= 0.8,   label='residual 1d mean')

# plt.plot(lat, + torge_zm* cos_phi,'g-', linewidth= 1,   label='M tourque')
# plt.plot(lat, (pdata['LHS'].isel(time=tt) + pdata['F_tur_zm'].isel(time=tt) - torge_zm )* cos_phi,'g--', linewidth= 1.3,   label='residual of 1 timestep')
#
#plt.plot(lat, (pdata['LHS'].isel(time=tt) + pdata['F_tur_zm'].isel(time=tt) )* cos_phi,'r--', linewidth= 1.3,   label='residual of 1 timestep')

plt.legend(loc='upper right', ncol=2)
plt.ylabel('$N /m^2$')
plt.xlim(xlims)
plt.grid()

plt.subplot(5,1,5)

plt.title('Eddy Fluxes')
# plt.plot(lat, pdata['LHS'].isel(time=tt)* cos_phi,c='k',  label=' LHS')
# plt.plot(lat, - pdata['F_tur_zm'].isel(time=tt)* cos_phi ,'k--', linewidth= 2,   label='F')
#plt.plot(lat, (pdata['LHS'].isel(time=tt) + pdata['F_tur_zm'].isel(time=tt) )* cos_phi,'r', linewidth= 1.3,   label='residual')

plt.plot(lat, pdata['uzm_vzm'].mean('time')* cos_phi,c='g',  label='Northward Fluxes by mean')
plt.plot(lat, pdata['uprime_vprime_zm'].mean('time')* cos_phi ,'b', linewidth= 2,   label='Northward Eddy-Flux')

plt.legend(loc='upper right', ncol=2)
plt.ylabel(pdata['uprime_vprime_zm'].units)

plt.xlim(xlims)
plt.grid()

#date_str =str(pdata.time[0].data.astype('M8[D]'))
F.save(name='exmpl_repres_dmean_ps_iews_'+ date_str, path=plot_path)


print(' all done and saved')
tend=time.time()
print('plot time :' + str(tend- tstart))
hist = logger(hist, 'plot time :' + str(tend- tstart))
hist = logger(hist, 'total time :' + str(tend-tstart_start) )

MT.save_log_txt(log_name ,log_path, hist , verbose=True )

#shutil.move(conf_path+'../active/'+date_str+'.json' , conf_path+'../processed/'+date_str+'.json', )
data_conf['processed_tstamp']= datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
data_conf['processed_time']= 'total time :' + str(tend-tstart_start)

try:
    del data_conf['conf_surface']
    del data_conf['conf_level']
except:
    pass

MT.json_save(name= date_str, path=conf_path+'../processed/' , data=data_conf)
exit()
