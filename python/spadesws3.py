import rasterio
import numpy as np
import pandas as pd
#import pickle
#import shutil

import ws3


def bootstrap_areas(fm, age_raster, hdt_raster, hdt_decode, pxa=1.):
    
  
    print('bootstrap_areas', basenames)
    if not year:
        for bn in basenames:
            print('copying', '%s/inventory_init.tif' % rst_path(bn), 
                  '%s/inventory_%i.tif' % (rst_path(bn), fm.base_year))
            shutil.copyfile('%s/inventory_init.tif' % rst_path(bn), 
                            '%s/inventory_%i.tif' % (rst_path(bn), fm.base_year))
        year = fm.base_year
    for dt in fm.dtypes.values(): # yuck
        dt.reset_areas(0)
        dt.reset_areas()
    for bn in basenames:
        _sumarea = 0.
        with rasterio.open('%s/inventory_%i.tif' % (rst_path(bn), year), 'r') as src:
            pxa = pow(src.transform.a, 2) * 0.0001 # pixel area (hectares)
            bh, ba = src.read(1), src.read(2)
            for h, dt in hdt[bn].items():
                ra = ba[np.where(bh == h)] # match themes hash value
                if new_dts:
                    fm.dtypes[dt] = ws3.forest.DevelopmentType(dt, fm)
                for age in np.unique(ra):
                    area = len(ra[np.where(ra == age)]) * pxa
                    _sumarea += area
                    fm.dtypes[dt].area(0, age, area)
        print('bootstrap_areas', bn, year, pxa, _sumarea)


def bootstrap_forestmodel(model_name, 
                          model_path, 
                          base_year, 
                          age_raster,
                          hdt_raster,
                          hdt_map,
                          hdt_func=ws3.common.hash_dt,
                          tif_path,
                          blk_raster=None,
                          horizon=ws3.common.HORIZON_DEFAULT, 
                          period_length=ws3.common.PERIOD_LENTH_DEFAULT, 
                          max_age=ws3.common.MAX_AGE_DEFAULT,
                          pxa=1.): 
    """
    Instantiate a new `ws3.forest.ForestModel` object, and initialize it
    from Woodstock-formatted input files. 

    :param numpy.ndarray age_raster: 2D `numpy.ndarray` of `numpy.int32` values 
        reprenting stand age values in forest initial inventory. 
    :param numpy.ndarray hdt_raster: 2D `numpy.ndarray` of `numpy.int32` values 
        reprenting hashed development type key tuple values in forest initial inventory. 
    :param dict hdt_map: Development type hash value decoding map, in the form
        returned by `ws3.common.hash_dt`.
    :param function hdt_func: Function used to hash development type keys.
    :param string tif_path: Path to directory where we can write temporary 
    :param numpy.ndarray blk_raster: 2D `numpy.ndarray` of `numpy.int32` values 
        reprenting harvest block ID values. If not specified, will be automatically 
        generated from hdt and age raster data.
    :param str model_name: The name of model.
    :param str model_path: The path to input data of model.
    :param int base_year: The base year of teh model.
    :param int horizon: The simulation horizon of the model.
    :param int horizon: The length of the simulation period.
    :param int max_age: The maximum age considered in the model.
    
    Initial inventory data is provided in the form of two 2D arrays of raster pixel values
    (representing hashed development type keys and stand ages) and a `dict`
    used to expand hash codes into development type keys (i.e., tuples of string values).
    
    raster pixel 
    data representing hashed `ws3.forest.DevelopmentType` key tuple values
    """
    fm = ws3.forest.ForestModel(model_name=model_name,
                                model_path=model_path,
                                base_year=base_year,
                                horizon=horizon,
                                period_length=period_length,
                                max_age=max_age)
    fm.import_landscape_section()
    bootstrap_areas(fm, age_raster, hdt_raster, hdt_decode, pxa=1.)
    fm.import_yields_section()
    fm.import_actions_section()
    fm.import_transition_section()
    #fm.reset() # do we need this?
    src_path = '%s/inventory_%i.tif' % (tif_path, fm.base_year)
    acode_map = {acode:'projected_%s' % acode for acode in fm.actions.keys() if acode not in ['null']}
    fr = ws3.spatial.ForestRaster(hdt_map=hdt_map,
                                  hdt_func=hdt_func,
                                  src_path=src_path,
                                  snk_path=snk_path,
                                  acode_map=acode_map,
                                  forestmodel=fm,
                                  base_year=fm.base_year,
                                  horizon=fm.horizon,
                                  period_length=fm.period_length,
                                  time_step=1)
    fm.forestraster = fr
    return fm
  
  
