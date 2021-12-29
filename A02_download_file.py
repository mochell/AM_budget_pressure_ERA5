"""
python 3
this script downlaods the file given by the date in the argument.
"""

import os
import cdsapi # configure cdsapi accordingly

import glob
import xarray as xr
import time
import datetime as datetime
import shutil
import sys

root_path= '/Users/Shared/Projects/codes/AM_budget_ERA5'
sys.path.append(root_path)

import tools_AM_budget as M
from tools_AM_budget import write_log as logger

save_path_base      = root_path+'/data/'
save_path_base_gph  = root_path+'/data/'
conf_path           = root_path
key_str='ERA5_subdaily_plev_.25deg_'

data_conf = dict()
data_conf['conf_surface'] = {
        'product_type':'reanalysis',
        'date':date_str,
        'time':[
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
            ],
        "param": "134/129/229/235045",
        "grid": "0.25/0.25",
        "area": "90/-180/-90/180",
        "format" : "netcdf" }

data_conf['conf_level'] = {
    'product_type':'reanalysis',
    'variable':[
        'u_component_of_wind','v_component_of_wind', 'geopotential'
    ],
    'pressure_level':[
        '1','2','3',
        '5','7','10',
        '20','30','50',
        '70','100','125',
        '150','175','200',
        '225','250','300',
        '350','400','450',
        '500','550','600',
        '650','700','750',
        '775','800','825',
        '850','875','900',
        '925','950','975',
        '1000'
    ],
    'date':date_str,
    "area": "90/-180/-90/180", #"-20/-180/-75/180",
    'time':[
        '00:00','01:00','02:00',
        '03:00','04:00','05:00',
        '06:00','07:00','08:00',
        '09:00','10:00','11:00',
        '12:00','13:00','14:00',
        '15:00','16:00','17:00',
        '18:00','19:00','20:00',
        '21:00','22:00','23:00'
    ],
"format" : "netcdf" }


print( 'download ' + date_str) #date_list[0]
date_str_save   = date_str.replace('-','')
load_filename   = date_str_save
#M.save_log_txt('worker_'+load_filename,save_log,  '',  verbose=True)
# with open(save_log+'worker_'+load_filename+'.hist.txt','a') as f:
#         f.write('\n start '+ date_str)

save_path       = save_path_base+key_str+date_str_save+'.nc'
save_path_srf   = save_path_base_gph+key_str+'gph_'+date_str_save+'.nc'

try:

    c2 = cdsapi.Client()
    c2.retrieve( 'reanalysis-era5-single-levels',    data_conf['conf_surface'],save_path_srf  )

    c = cdsapi.Client()
    c.retrieve(
        'reanalysis-era5-pressure-levels', data_conf['conf_level'],save_path )

    dtime= time.clock() - start
    print('saved at '+ save_path)
    print('time taken ' + str(dtime))
    #time.sleep(5)

    # with open(save_log+'worker_'+load_filename+'.hist.txt','a') as f:
    #         f.write('\n saved at '+ save_path)
    #         f.write('\n time taken ' + str(dtime))

    # save requested file to the stock
    data_conf['download_tstamp']= datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    data_conf['save_path_level']= save_path
    data_conf['save_path_surface']  = save_path_srf
    M.json_save(name= date_str, path=conf_path+'/downloaded/' , data=data_conf)
    #os.remove(conf_path+'../downloaded/active/'+date_str+'.json')

except:
    print('failed submiting')
    #shutil.move(conf_path+'../downloaded/active/'+date_str+'.json' , conf_path+'../requested/'+date_str+'.json')
    # with open(save_log+'worker_'+load_filename+'.hist.txt','a') as f:
    #         f.write('\n failed at ' + date_str)
