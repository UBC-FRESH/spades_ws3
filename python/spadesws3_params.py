import os
from os import listdir
from os.path import isfile, join, dirname
import sys

########################################################################################################
# import ws3 from local copy
#if use_local_ws3:
#    ws3_path = os.path.abspath('./ws3')
#    _path = os.path.abspath(join(dirname(__file__), ws3_path))

#ws3_path = os.path.abspath(join(r.getPaths()['modulePath'], r.currentModule()))
#ws3_path = os.path.abspath(join(r.getPaths()['modulePath'], r.currentModule()))
#try:
#_path = os.path.abspath(join(dirname(__file__), ws3_path))
#except:
#    _path = os.path.abspath(join(dirname('__file__'), ws3_path))
#if not _path in os.sys.path: os.sys.path.insert(0, _path)
import ws3
from ws3.forest import ForestModel, Action
from ws3.spatial import ForestRaster
from ws3.common import clean_vector_data, reproject_vector_data, rasterize_stands, hash_dt, warp_raster
########################################################################################################

print(os.getcwd())
print(os.sys.path)
#os.sys.path.insert(0, os.path.abspath('./src'))

import numpy as np
import pandas as pd
from pandas import DataFrame as DF
import geopandas as gpd
import matplotlib.pyplot as plt
import gurobipy as grb
import rasterio
import rasterio.plot
import seaborn as sns
from functools import partial
import fiona
import folium
import pathlib
try:
   import cPickle as pickle
except:
   import pickle
import copy
from pathlib import Path
import rasterio
import shutil

from spadesws3 import clean_shapefiles, rasterize_inventory, read_basenames, compile_basecodes, bootstrap_forestmodel, bootstrap_areas, schedule_harvest, sda

# configure paths and global variables
scenario_name = 'base'

sda_mode = 'randblk' # 'randpxl'
obj_mode = 'min_harea' # 'max_harea'
horizon = 2
period_length = 1
yields_period_length = 10
yields_x_unit = 'years'
time_step = 1
max_age = 1000
try:
    dat_path = os.path.abspath(r.getPaths()['inputPath'])
except:
    dat_path = '../../../input'
target_path = join(dat_path, 'targets.csv')
yld_path = '%s/yld.csv' % dat_path
tolerance = 10.
clean_inv = False
rasterize_inv = False
gdb_path = lambda bn: '%s/gis/gdb/%s.gdb' % (dat_path, bn)
shp_path = lambda bn: '%s/gis/shp/%s.shp' % (dat_path, bn)
tif_path = lambda bn: '%s/tif/%s' % (dat_path, bn)
shp_name = 'stands'
age_col = 'age'
theme_cols = ['theme0', 'theme1', 'theme2', 'theme3']
compress = 'lzw'
dtype = rasterio.uint8
base_year = 2015
prop_names = [u'THLB', u'AU', u'LdSpp', u'Age2015', u'Shape_Area']
prop_types = [(u'theme0', 'str:10'),
              (u'theme1', 'str:1'),
              (u'theme2', 'str:5'), 
              (u'theme3', 'str:50'), 
              (u'age', 'int:5'), 
              (u'area', 'float:10.1')]
tvy_name = 'totvol'
snk_epsg = 3005 # ESPG:3005 corresponds to NAD83/BC Albers
sns.set_style('dark')
hdt_path = '%s/hdt' % dat_path
#coast_bn = ['tsa01', 'tsa02', 'tsa03'] # FIX ME: totally bogus (check GIS data for real coast TSA codes) 
#_p_slashburn = read_pslashburn('%s/p_slashburn.csv' % dat_path)
#p_slashburn = lambda bn: _p_slashburn[bn]
# oe = operability expression
# ogi = old growth index
oe_harvest = '_age >= 40 and _age <= 999'
action_params = {'harvest':{'oe':oe_harvest,
                            'mask':('?', '1', '?', '?'),
                            'is_harvest':True,
                            'targetage':0}}
util = 0.85

proc_mode = 'seq_tsa'
deg_mode = 'classic'
obj_mode = 'min_harea'
run_sda = 1 # set to True for production runs
fm_horizon = horizon
cap_age = 900

twopass_cacut_factor = 1.50
run_sda = True
#basenames = ['tsa31'] #, 'tsa31']

raster_d = 250

kwargs = []
if 1:
    hdt = {bn:pickle.load(open('%s/hdt_%s.pkl' % (hdt_path, bn), 'rb')) for bn in basenames}
    basecodes = compile_basecodes(hdt, basenames, theme_cols)
    kwargs.append({'basenames':basenames,
                   'model_name':'foo',
                   'model_path':dat_path,
                   'base_year':base_year,
                   'yld_path':yld_path,
                   'tif_path':tif_path,
                   'horizon':horizon,
                   'period_length':period_length,
                   'max_age':max_age,
                   'basecodes':basecodes,
                   'action_params':action_params,
                   'hdt':hdt,
                   'add_null_action':True,
                   'tvy_name':tvy_name,
                   'compile_actions':True,
                   'yields_x_unit':yields_x_unit,
                   'yields_period_length':yields_period_length,
                   'verbose':1})
            

    
else:    
    for i, bn in enumerate(basenames):
        hdt = {bn:pickle.load(open('%s/hdt_%s.pkl' % (hdt_path, bn), 'rb'))}
        kwargs.append({'basenames':[bn],
                       'model_name':'foo',
                       'model_path':dat_path,
                       'base_year':base_year,
                       'yld_path':yld_path,
                       'tif_path':tif_path,
                       'horizon':horizon,
                       'period_length':period_length,
                       'max_age':max_age,
                       'basecodes':compile_basecodes(hdt, [bn], theme_cols),
                       'action_params':action_params,
                       'hdt':hdt,
                       'add_null_action':True,
                       'tvy_name':tvy_name,
                       'compile_actions':True,
                       'yields_x_unit':yields_x_unit,
                       'yields_period_length':yields_period_length,
                       'verbose':1})


####
#import ipyparallel as ipp
#import os
#import time
#c = ipp.Client()
#c[:].use_cloudpickle()
#v = c.load_balanced_view()
#v.map(os.chdir, [os.getcwd()]*len(c.ids))
#amr = v.map_async(lambda kwargs: do_the_thing(**kwargs), kwargs)
####

def bootstrap_forestmodel_kwargs():
    return bootstrap_forestmodel(**kwargs[0])


def simulate_harvest(fm, basenames, year):
    bootstrap_areas(fm, basenames, tif_path, hdt, year, new_dts=False)
    #fm.compile_actions()
    #fm.reset_actions()
    fm.initialize_areas()
    #fm.grow()
    schedule_harvest(fm, basenames, target_path=target_path)
    sda(fm, basenames, 1, tif_path, hdt, sda_mode=sda_mode)


#def bootstrap_forestmodel_from_kwargslist(i=0):
#    return bootstrap_forestmodel(**kwargs[int(i)])

#def schedule_harvest_kwargslist(fm, i=0):
#    return schedule_harvest(fm, kwargs[0]['basenames'], target_path=target_path)
    
#def run_sda_kwargslist(fm, i=0):
#    return sda(fm, basenames, 1, kwargs[i]['tif_path'], kwargs[i]['hdt'])

