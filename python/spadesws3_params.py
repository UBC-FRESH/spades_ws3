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
#from ws3.common import clean_vector_data, reproject_vector_data, rasterize_stands, hash_dt, warp_raster
from ws3.common import clean_vector_data, reproject_vector_data, hash_dt, warp_raster
########################################################################################################

#print(os.getcwd())
#print(os.sys.path)

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

from spadesws3 import clean_shapefiles, rasterize_inventory, read_basenames, compile_basecodes, bootstrap_forestmodel, bootstrap_areas, schedule_harvest_optimize, schedule_harvest_areacontrol, sda

# configure paths and global variables
scenario_name = 'base'

sda_mode = 'randblk' # 'randpxl'
obj_mode = 'min_harea' # 'max_harea'
horizon = 2
period_length = 10
yields_period_length = 10
yields_x_unit = 'years'
time_step = 1
max_age = 1000
try:
    dat_path = os.path.abspath(r.getPaths()['inputPath'])
except:
    dat_path = '../../../input'
target_path = join(dat_path, 'targets.csv')
#yld_path = dat_path # '%s/yld.csv' % dat_path
tolerance = 10.
clean_inv = False
rasterize_inv = False
shp_path = lambda bn: '%s/gis/shp/%s.shp' % (dat_path, bn)
tif_path = lambda bn: '%s/tif/%s' % (dat_path, bn)
yld_path = '%s/yld' % dat_path
shp_name = 'stands'
age_col = 'age'
theme_cols = ['fma', 'thlb', 'ygrp', 'bcov']
compress = 'lzw'
dtype = rasterio.uint8
base_year = 2015
prop_names = [u'fma', u'thlb', u'ygrp', u'bcov', u'stAge', u'area_ha_pl']
prop_types = [(u'theme0', 'str:5'),
              (u'theme1', 'str:1'),
              (u'theme2', 'str:20'), 
              (u'theme3', 'str:5'), 
              (u'age', 'int:5'), 
              (u'area', 'float:10.1')]
tvy_name = 'totvol'
snk_epsg = 3005 # ESPG:3005 corresponds to NAD83/BC Albers
sns.set_style('dark')
hdt_path = '%s/hdt' % dat_path
oe_harvest = '_age >= 40 and _age <= 999'
oe_fire = '_age >= 100 and _age <= 999'
action_params = {'harvest':{'oe':oe_harvest,
                            'mask':('?', '1', '?', '?'),
                            'is_harvest':True,
                            'targetage':0},
                 'fire':{'oe':oe_fire,
                            'mask':('?', '?', '?', '?'),
                            'is_harvest':False,
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

raster_d = 90

try:
    hdt = {bn:pickle.load(open('%s/hdt_%s.pkl' % (hdt_path, bn), 'rb')) for bn in basenames}
except:
    hdt = {}
    
def kwargs():
    basecodes = compile_basecodes(hdt, basenames, theme_cols)
    kwargs = {'basenames':basenames,
              'model_name':'foo',
              'model_path':dat_path,
              'base_year':int(base_year),
              'yld_path':yld_path,
              'tif_path':tif_path,
              'horizon':int(horizon),
              'period_length':int(period_length),
              'max_age':int(max_age),
              'basecodes':basecodes,
              'action_params':action_params,
              'hdt':hdt,
              'add_null_action':True,
              'tvy_name':tvy_name,
              'compile_actions':True,
              'yields_x_unit':yields_x_unit,
              'yields_period_length':int(yields_period_length),
              'verbose':1}
    return kwargs    


def bootstrap_forestmodel_kwargs():
    return bootstrap_forestmodel(**kwargs())


def simulate_harvest(fm, basenames, year,
                     mode='optimize',
                     scenario_name='base',
                     cbird=None,
                     target_masks=None, 
                     target_areas=None,
                     target_scalefactors=None,
                     mask_area_thresh=0.,
                     verbose=False):
    bootstrap_areas(fm, basenames, tif_path, yld_path, hdt, year, new_dts=False)
    fm.reset()
    print('cbird1', cbird)
    if mode == 'optimize':
        schedule_harvest_optimize(fm, basenames, scenario_name=scenario_name, target_scalefactors=target_scalefactors, target_path=target_path, util=util, cbird=cbird)
    elif mode == 'areacontrol':
        schedule_harvest_areacontrol(fm, 
                                     target_masks=target_masks, 
                                     target_areas=target_areas, 
                                     target_scalefactors=target_scalefactors,
                                     mask_area_thresh=mask_area_thresh,
                                     verbose=verbose)
    else: # bad mode value
        raise ValueError('Bad mode value')
    acode_map = {acode:'projected_%s' % acode for acode in fm.actions.keys() if acode not in ['null']}
    sda(fm, basenames, 1, tif_path, hdt, sda_mode=sda_mode, acode_map=acode_map)
