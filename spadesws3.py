import rasterio
import numpy as np
import pandas as pd
import ws3.common, ws3.forest, ws3.spatial
import sys

this = sys.modules[__name__]

action_decode = {}

def init_action_decode(d):
    for acode, raster_value in d.items():
        this.action_decode[int(raster_value)] = acode
    print(this.action_decode)

#def map_raster_actions(ws3_acode, raster_value):
#    this.action_decode[int(raster_value)] = ws3_acode


def transition_raster_actions(raster_data, hdtk_layer_index=0, age_layer_index=1, action_layer_index=3):
    xk, xg, xc = raster_data[hdtk_layer_index], raster_data[age_layer_index], raster_data[action_layer_index]
    for vc in np.unique(xc):
        if vc not in this.action_decode: 
            continue
        else:
            acode = this.action_decode[vc]
        ic = np.where(xc == vc)
        #print()
        #print('found raster action', acode, ic)
        #print(xg)
        for vk in np.unique(xk[ic]):
            source_dtk = this.hdtk_decode[vk]
            ick = np.where(xk[ic] == vk)
            for vg in np.unique(xg[ic][ick]):
                source_age = vg
                ickg = np.where((xc == vc) & (xk == vk) & (xg == vg))
                if source_dtk not in this.forestmodel.dtypes:
                    this.forestmodel.dtypes[source_dtk] = this.forestmodel.create_dtype_fromkey(source_dtk)
                error_code, _, target_dt = this.forestmodel.apply_action(source_dtk, acode, 1, source_age, 0.)
                if error_code > 0: continue
                #assert error_code == 0 # else badness
                assert len(target_dt) == 1 and target_dt[0][1] == 1. # else badness
                target_dtk, _, target_age = target_dt[0]
                target_hdtk = ws3.common.hash_dt(target_dtk)
                #print(ickg)
                #print(vk, target_hdtk, source_age, target_age)
                #print(source_dtk, this.hdtk_decode[target_hdtk])
                #print('before value update', xk[ickg], xg[ickg])
                xk[ickg] = target_hdtk
                xg[ickg] = target_age
                #print('after value update', xk[ickg], xg[ickg])


def import_inventory(shp_path, tif_path, tif_filename, theme_cols, age_col='age', age_divisor=1., 
                     cap_age=None, d=100., verbose=False):
    this.hdtk_decode = ws3.common.rasterize_stands(shp_path=shp_path, 
                                                   tif_path='%s/%s' % (tif_path, tif_filename),
                                                   theme_cols=theme_cols,
                                                   age_col=age_col,
                                                   age_divisor=age_divisor,
                                                   cap_age=cap_age,
                                                   d=d,
                                                   extra_bands=1, # for SpaDES modules to leave action requests
                                                   verbose=verbose)


def update_areas(age_raster, hdtk_raster, pixel_area=1.):
    for hashcode, dtk in this.hdtk_decode.items():
        ages = age_raster[np.where(hdtk_raster == hashcode)]
        this.forestmodel.dtypes[dtk] = this.forestmodel.create_dtype_fromkey(dtk)
        #if len(ages) and dtk not in this.forestmodel.dtypes: # create new development types as needed
        #    #this.forestmodel.dtypes[dtk] = ws3.forest.DevelopmentType(dtk, this.forestmodel)
        #    this.forestmodel.dtypes[dtk] = this.forestmodel.create_dtype_fromkey(dtk)
        for age in np.unique(ages):
            area = len(ages[np.where(ages == age)]) * pixel_area
            this.forestmodel.dtypes[dtk].area(0, age, area)


def raster_memfile(raster_data, raster_profile):
    memfile = rasterio.MemoryFile()
    with memfile.open(**raster_profile) as snk:
        snk.write(raster_data)
    return memfile


def bootstrap_forestmodel(model_name, 
                          model_path, 
                          base_year, 
                          tif_path,
                          tif_filename,
                          hdt_func=ws3.common.hash_dt,
                          horizon=ws3.common.HORIZON_DEFAULT, 
                          period_length=ws3.common.PERIOD_LENGTH_DEFAULT, 
                          max_age=ws3.common.MAX_AGE_DEFAULT):
    this.forestmodel = ws3.forest.ForestModel(model_name=model_name,
                                              model_path=model_path,
                                              base_year=base_year,
                                              horizon=horizon,
                                              period_length=period_length,
                                              max_age=max_age)
    this.forestmodel.import_landscape_section()
    this.forestmodel.import_yields_section()
    this.forestmodel.import_actions_section()
    this.forestmodel.import_transitions_section()
    this.forestmodel.reset_actions()
    print('reading inventory raster data', tif_filename)
    with rasterio.open('%s/%s' % (tif_path, tif_filename), 'r') as src:
        raster_data = src.read()
        raster_profile = src.profile
        ## debug
        this.__raster_data = raster_data
        this.__raster_profile = raster_profile
        #assert False
    transition_raster_actions(raster_data)
    t = raster_profile['transform']
    pixel_area = t[0] * t[4] * -0.0001 # hectares, assuming pixel dimensions in meters
    update_areas(raster_data[1], raster_data[0], pixel_area)
    this.forestmodel.reset()
    acode_map = {acode:'projected_%s' % acode for acode in this.forestmodel.actions.keys() if acode not in ['null']}
    this.forestmodel.forestraster = ws3.spatial.ForestRaster(hdt_map=this.hdtk_decode,
                                                             hdt_func=hdt_func,
                                                             src_path=raster_memfile(raster_data, raster_profile), #'%s/%s' % (tif_path, tif_filename), # to do: swap out with memfile
                                                             snk_path=tif_path,
                                                             acode_map=acode_map,
                                                             forestmodel=this.forestmodel,
                                                             base_year=base_year, # to do: tweak ForestRaster constructor to skip this (read from fm)
                                                             horizon=1,
                                                             period_length=1,
                                                             time_step=1,
                                                             piggyback_acodes={},
                                                             spades_mode=True)
  
def simulate_actions(scheduler_mode, targets=None):
    match scheduler_mode:
        case 'heuristic':
            schedule_actions_heuristic(targets)
        case 'optimize':
            schedule_actions_optimize()
        case _:
            print('invalid scheduler mode')
            assert False # badness


def schedule_actions_heuristic(targets):
    for target in targets:
        #this.forestmodel.areaselector.operate(1, target['acode'], target['area'], target['mask'])
        this.forestmodel.areaselector.operate(1, verbose=1, **target)
    this.forestmodel.forestraster.allocate_schedule()


def schedule_actions_optimize():
    pass


##################################################
# below this line: cut-and-paste stuff or junk  
def schedule_harvest_heuristic(fm, period=1, acode='harvest', util=0.85, 
                               target_masks=None, target_areas=None, target_scalefactors=None,
                               mask_area_thresh=0.,
                               verbose=0):
    fm.reset_actions()
    if not target_areas:
        if not target_masks: # default to AU-wise THLB 
            au_vals = []
            au_agg = []
            for au in fm.theme_basecodes(2):
                mask = '? 1 %s ?' % au
                masked_area = fm.inventory(0, mask=mask)
                if masked_area > mask_area_thresh:
                    au_vals.append(au)
                else:
                    au_agg.append(au)
                    if verbose > 0:
                        print('adding to au_agg', mask, masked_area)
            if au_agg:
                fm._themes[2]['areacontrol_au_agg'] = au_agg 
                au_vals.append('areacontrol_au_agg')
            target_masks = ['? 1 %s ?' % au for au in au_vals]
        #print(target_masks)
        #assert False
        target_areas = []
        for i, mask in enumerate(target_masks): # compute area-weighted mean CMAI age for each masked DT set
            masked_area = fm.inventory(0, mask=mask, verbose=verbose)
            if not masked_area: continue
            r = sum((fm.dtypes[dtk].ycomp('totvol').mai().ytp().lookup(0) * fm.dtypes[dtk].area(0)) for dtk in fm.unmask(mask))
            r /= masked_area
            #awr = []
            #dtype_keys = fm.unmask(mask)
            #for dtk in dtype_keys:
            #    dt = fm.dtypes[dtk]
            #    awr.append(dt.ycomp('totvol').mai().ytp().lookup(0) * dt.area(0))
            #r = sum(awr)  / masked_area
            asf = 1. if not target_scalefactors else target_scalefactors[i]  
            ta = (1/r) * masked_area * asf
            target_areas.append(ta)
    for mask, target_area in zip(target_masks, target_areas):
        if verbose > 0:
            print('calling areaselector', period, acode, target_area, mask)
        fm.areaselector.operate(period, acode, target_area, mask=mask, verbose=verbose)
    sch = fm.compile_schedule()
    return sch