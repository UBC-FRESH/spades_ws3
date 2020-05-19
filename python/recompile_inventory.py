#!/usr/bin/env python
# coding: utf-8

# This notebook demonstrates the geodata processing pipeline functions built into `ws3` and `spadesws3`.
# 
# Basically, the geodata processing pipeline works in two stages. 
# 
# The first stage imports vegetation resource inventory (VRI) vector polygon data and cleans the data (i.e., fix problems with and potentially simplify polygon geometry, discard unwanted attribute data columns, rename remaining attribute data columns, patch bad data values), reproject vectory geometry to target CRS, and exports clean copies of reprojected vector data layers to disk.  
# 
# The second stage imports the cleaned vector datasets and compiles multi-band GeoTIFF files that are compatible with input data requirements of the `spades_ws3` SpaDES module (i.e., the harvesting module).

# Set a few key variables before sourcing the `spadesws3_params.py` module code. We do not really need to source most of the code in that module, but this way we are certain that the Python environment we are using is almost identical to the environment set up inside the `reticulate` Python bubble in SpaDES when we load the `spades_ws3` module.

# In[1]:


dat_path = '../../../input'
from spadesws3 import read_basenames
basenames = read_basenames(dat_path+'/basenames.txt')
get_ipython().run_line_magic('run', '-i spadesws3_params')


# Before running this notebook, you should clone the dev branch of `ws3` and use setuptools "development mode" to deploy symlinks to this cloned repo to `opt/tljh/user/lib/python3.6/site-packages/` by running `sudo -H python setup.py develop` from inside the cloned `ws3` directory, _from a Jupyter terminal_ (TLJH uses its own Python environment, so any deps have to be installed to there not to the system Python environment). 

# In[2]:


ws3.__path__


# In[3]:


#basenames = ['tsa08', 'tsa16', 'tsa24', 'tsa40', 'tsa41'] # RIA landbase
#basenames = ['tsa10'] # Kalum TSA is the smallest in the set, so good for testing.
snk_epsg = 3005 # BC Albers
tolerance = 10.
prop_names = [u'THLB', u'AU', u'LdSpp', u'Age2015', u'Shape_Area']
prop_types = [(u'theme0', 'str:10'),
              (u'theme1', 'str:1'),
              (u'theme2', 'str:5'), 
              (u'theme3', 'str:50'), 
              (u'age', 'int:5'), 
              (u'area', 'float:10.1')]
update_area_prop = 'area'
do_clean_shapefiles = True
#pixel_width = 100.
pixel_width = 250.


# Run stage 1 (i.e., clean the original VRI vector datasets).
# 
# `clean_shapefiles` is implemented in local module `spadesws3.py`, which basically just calls `ws3.common.clean_vector_data` for reach TSA in `basenames`. 

# In[4]:


if do_clean_shapefiles:
    clean_shapefiles(basenames, gdb_path, shp_path, snk_epsg, prop_names, prop_types, tolerance, update_area_prop)


# Run stage 2 (i.e., compile initial inventory multi-band GeoTIFF files for use with `spades_ws3` SpaDES module).
# 
# `rasterize_inventory` is implemented in local module `spadesws3.py`, which basically just calls `ws3.common.rasterize_stands` for reach TSA in `basenames`. 

# In[5]:


_ = rasterize_inventory(basenames, shp_path, tif_path, hdt_path, theme_cols, age_col, period_length, base_year,
                    cap_age=None, d=pixel_width, verbose=True)


# In[ ]:




