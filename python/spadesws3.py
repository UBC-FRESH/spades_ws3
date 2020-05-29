import rasterio
import numpy as np
import pandas as pd
import shutil
try:
    import cPickle as pickle
except:
    import pickle
import ws3

def read_basenames(path):
    return [line.lower().strip().split(' ')[0] 
        for line in open(path, 'r') if not line.startswith('#')]
            
def cmp_c_z(fm, path, expr):
    """
    Compile objective function coefficient (given ForestModel instance, 
    leaf-to-root-node path, and expression to evaluate).
    """
    result = 0.
    for t, n in enumerate(path, start=1):
        d = n.data()
        if fm.is_harvest(d['acode']):
            result += fm.compile_product(t, expr, d['acode'], [d['dtk']], d['age'], coeff=False)
    return result


def cmp_c_cflw(fm, path, expr, mask=None): # product, all harvest actions
    """
    Compile flow constraint coefficient (given ForestModel instance, 
    leaf-to-root-node path, expression to evaluate, and optional mask).
    """
    result = {}
    for t, n in enumerate(path, start=1):
        d = n.data()
        if mask and not fm.match_mask(mask, d['dtk']): continue
        if fm.is_harvest(d['acode']):
            result[t] = fm.compile_product(t, expr, d['acode'], [d['dtk']], d['age'], coeff=False)
    return result


def cmp_c_caa(fm, path, expr, acodes, mask=None): # product, named actions
    """
    Compile flow constraint coefficient (given ForestModel instance, 
    leaf-to-root-node path, expression to evaluate, list of action codes, and optional mask).
    """
    result = {}
    for t, n in enumerate(path, start=1):
        d = n.data()
        if mask and not fm.match_mask(mask, d['dtk']): continue
        if d['acode'] in acodes:
            result[t] = fm.compile_product(t, expr, d['acode'], [d['dtk']], d['age'], coeff=False)
    return result

def _gen_scen_base(fm, basenames, name='base', util=0.85, param_funcs=None, harvest_acode='harvest',  
                   tvy_name='totvol', toffset=0, obj_mode='max_hvol', target_path='./input/targets.csv',
                   max_tp=2074, cacut=None, mask=None):
    from functools import partial
    acodes = ['null', harvest_acode]  
    vexpr = '%s * %0.2f' % (tvy_name, util)
    if obj_mode == 'max_hvol':
        sense = ws3.opt.SENSE_MAXIMIZE 
        zexpr = vexpr
    elif obj_mode == 'min_harea':
        sense = ws3.opt.SENSE_MINIMIZE 
        zexpr = '1.'
    else:
        raise ValueError('Invalid obj_mode: %s' % obj_mode)
    if not param_funcs:
        df_targets = pd.read_csv(target_path).set_index(['tsa', 'year'])
        param_funcs = {}
        param_funcs['cvcut'] = lambda bn, t: float(df_targets.loc[bn, t]['vcut']) if t <= max_tp else float(df_targets.loc[bn, max_tp]['vcut'])
        param_funcs['cabrn'] = lambda bn, t: float(df_targets.loc[bn, t]['abrn']) if t <= max_tp else float(df_targets.loc[bn, max_tp]['abrn'])
        param_funcs['cflw_acut_e'] = lambda bn, t: df_targets.loc[bn, t]['cflw_acut_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cflw_acut_e']
        param_funcs['cgen_vcut_e'] = lambda bn, t: df_targets.loc[bn, t]['cgen_vcut_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cgen_vcut_e']
        param_funcs['cgen_acut_e'] = lambda bn, t: df_targets.loc[bn, t]['cgen_vcut_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cgen_vcut_e']
        param_funcs['cgen_abrn_e'] = lambda bn, t: df_targets.loc[bn, t]['cgen_abrn_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cgen_abrn_e']
    coeff_funcs = {'z':partial(cmp_c_z, expr=zexpr)}
    coeff_funcs.update({'cacut_%s' % bn:partial(cmp_c_caa, expr='1.', acodes=[harvest_acode], mask=(bn, '?', '?', '?')) 
                        for bn in basenames})
    coeff_funcs.update({'cvcut_%s' % bn:partial(cmp_c_caa, expr=vexpr, acodes=[harvest_acode], mask=(bn, '?', '?', '?')) 
                        for bn in basenames})
    T = fm.periods# [fm.base_year+(t-1)*fm.period_length for t in fm.periods]
    cflw_e, cgen_data = {}, {}
    #foo = {bn:{t:(bn, t+toffset) for t in T} for bn in basenames}
    #print(T)
    #assert False
    #cflw_ebn = {bn:({t:param_funcs['cflw_acut_e'](bn, fm.base_year+(t-1)*fm.period_length+toffset) for t in T}, 1) for bn in basenames}
    #cflw_e.update({'cacut_%s'%bn:cflw_ebn[bn] for bn in basenames})
    for bn in basenames:
        #print(df_targets.loc[bn])
        cgen_data.update({'cvcut_%s' % bn:{'lb':{t:param_funcs['cvcut'](bn, fm.base_year+(t-1)*fm.period_length+toffset) *
                                                 (1. - param_funcs['cgen_vcut_e'](bn, fm.base_year+(t-1)*fm.period_length+toffset))
                                               for t in T}, 
                                         'ub':{t:param_funcs['cvcut'](bn, fm.base_year+(t-1)*fm.period_length+toffset) for t in T}}})
        if cacut:
            cgen_data.update({'cacut_%s' % bn:{'lb':{t:param_funcs['cacut'](bn, fm.base_year+(t-1)*fm.period_length)*
                                                   (1. - param_funcs['cgen_acut_e'](bn, fm.base_year+(t-1)*fm.period_length)) for t in T}, 
                                             'ub':{t:param_funcs['cacut'](bn, fm.base_year+(t-1)*fm.period_length) for t in T}}})
    #print(cflw_e)
    fm._tmp = {}
    fm._tmp['param_funcs'] = param_funcs
    fm._tmp['cgen_data'] = cgen_data
    return fm.add_problem(name, coeff_funcs, cflw_e, cgen_data=cgen_data, acodes=acodes, sense=sense, mask=mask)


def gen_scen(fm, basenames, name, util, param_funcs, toffset=0, obj_mode='max_hvol', cacut=None, mask=None, target_path='./input/targets.csv'):
    dsp = {'base':_gen_scen_base}
    return dsp[name](fm, basenames, name, util, param_funcs=param_funcs, toffset=toffset, obj_mode=obj_mode, cacut=cacut, mask=mask, target_path=target_path)


def unhash_ij(problem):
    r = {}
    for i, tree in problem.trees.items():
        for path in tree.paths():
            j = tuple(n.data('acode') for n in path)
            r['x_%i' % hash((i, j))] = i, j
    return r


def bootstrap_themes(fm, theme_cols=['theme0', 'theme1', 'theme2', 'theme3'], 
                     basecodes=[[], [], [], []], aggs=[{}, {}, {}, {}], verbose=False):
    for ti, t in enumerate(theme_cols):
        fm.add_theme(t, basecodes=basecodes[ti], aggs=aggs[ti])
    fm.nthemes = len(theme_cols)

    
def bootstrap_areas(fm, basenames, rst_path, hdt, year=None, new_dts=True):
    #fm.dtypes = {}
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
        #print('bootstrap_areas', bn, year)
        #print('%s/inventory_%i.tif' % (rst_path(bn), year))
        with rasterio.open('%s/inventory_%i.tif' % (rst_path(bn), year), 'r') as src:
            pxa = pow(src.transform.a, 2) * 0.0001 # pixel area (hectares)
            bh, ba = src.read(1), src.read(2)
            for h, dt in hdt[bn].items():
                #from IPython import embed; embed()
                ra = ba[np.where(bh == h)] # match themes hash value
                if new_dts:
                    fm.dtypes[dt] = ws3.forest.DevelopmentType(dt, fm)
                #else:
                #    fm.dtypes[dt].reset_areas() 
                for age in np.unique(ra):
                    area = len(ra[np.where(ra == age)]) * pxa
                    _sumarea += area
                    #print(bn, dt, age, area)
                    fm.dtypes[dt].area(0, age, area)
        print('bootstrap_areas', bn, year, pxa, _sumarea)

                    
def bootstrap_yields(fm, yld_path, theme_cols=['AU', 'LDSPP'], spcode='SPCode', 
                     startp_col='Wdks', x_max=360, y_cols=None, 
                     period_length=10, x_unit='periods', tvy_name='totvol'):
    y_cols = ['X%i' % i for i in range(0, x_max, period_length)]
    df = pd.read_csv(yld_path, usecols=theme_cols+[spcode, startp_col]+y_cols)
    df [theme_cols[1]] = df[theme_cols[1]].str.lower().str.replace(r'[- ]+', '_')
    for i in (0, 1): df[theme_cols[i]] = df[theme_cols[i]].astype(str)
    df = df.set_index(theme_cols)
    x_factor = 1 if x_unit == 'periods' else period_length
    period_length = period_length if x_unit == 'periods' else 1
    for t1, t2 in df.index.values: # assuming exactly one yield curve per unique combination of AU and LDSPP
        mask = ('?', '?', t1, t2)
        dt_keys = fm.unmask(mask)
        if not dt_keys: continue
        r = df.loc[t1, t2]
        yname = str.lower(r[spcode])
        points = [((x+r[startp_col])*x_factor, r[y]) for x, y in enumerate(y_cols)]
        c = fm.register_curve(ws3.core.Curve(yname, points=points, type='a', is_volume=True, xmax=fm.max_age, period_length=period_length))
        fm.yields.append((mask, 'a', [(yname, c)]))
        fm.ynames.add(yname)
        for dtk in dt_keys: fm.dtypes[dtk].add_ycomp('a', yname, c)
    # add total volume curve ###
    expr = '_SUM(%s)' % ', '.join(fm.ynames)
    fm.yields.append((('?', '?', '?', '?'), 'c', [(tvy_name, expr)]))
    fm.ynames.add(tvy_name)
    for dtk in fm.dtypes.keys(): fm.dtypes[dtk].add_ycomp('c', tvy_name, expr)


def bootstrap_actions(fm, action_params):
    for acode in action_params:
        ap = action_params[acode]
        mask, oe, is_harvest, targetage = ap['mask'], ap['oe'], ap['is_harvest'], ap['targetage']
        target = [(mask, 1.0, None, None, None, None, None)]
        fm.actions[acode] = ws3.forest.Action(acode, targetage=targetage, is_harvest=is_harvest)
        fm.oper_expr[acode] = {mask:oe}
        fm.transitions[acode] = {mask:{'':target}}
        for dtk in fm.unmask(mask):
            dt = fm.dtypes[dtk]
            dt.oper_expr[acode] = [oe]
            for age in range(1, fm.max_age):
                if not dt.is_operable(acode, 1, age): continue
                fm.dtypes[dtk].transitions[acode, age] = target

                
def bootstrap_forestmodel(basenames, model_name, model_path, base_year, yld_path, tif_path, horizon, 
                          period_length, max_age, basecodes, action_params, hdt,
                          add_null_action=True, tvy_name='totvol', compile_actions=True,
                          yields_x_unit='periods', yields_period_length=None, verbose=0):
    if not yields_period_length: yields_period_length = period_length
    from ws3.forest import ForestModel
    fm = ForestModel(model_name=model_name, 
                     model_path=model_path,
                     base_year=base_year,
                     horizon=horizon,     
                     period_length=period_length,
                     max_age=max_age)
    bootstrap_themes(fm, basecodes=basecodes)    
    bootstrap_areas(fm, basenames, tif_path, hdt)
    bootstrap_yields(fm, yld_path, tvy_name=tvy_name, period_length=yields_period_length, x_unit=yields_x_unit)
    bootstrap_actions(fm, action_params)
    if add_null_action: fm.add_null_action()
    fm.compile_actions()
    fm.reset_actions()
    fm.initialize_areas()
    fm.grow()
    return fm


def clean_shapefiles(basenames, gdb_path, shp_path, snk_epsg, prop_names, prop_types, tolerance, update_area_prop=''):
    import pathlib
    import fiona
    from ws3.common import clean_vector_data
    #from os import listdir, remove
    #from os.path import isfile, join
    for bn in basenames:
        print('cleaning GDB', gdb_path(bn))
        if not pathlib.Path(shp_path(bn)).exists(): 
            pathlib.Path(shp_path(bn)).mkdir()
        snk1_path, snk2_path = clean_vector_data(gdb_path(bn), shp_path(bn), 'stands', prop_names, 
                                                 tolerance=tolerance, max_records=None, 
                                                 theme0=bn, prop_types=prop_types, dst_epsg=snk_epsg, 
                                                 update_area_prop=update_area_prop)
        with fiona.open(gdb_path(bn)) as src0, fiona.open(snk1_path) as src1, fiona.open(snk2_path) as src2:
            print('Polygons in original dataset', len(src0))
            print('Polygons in clean dataset', len(src1))
            print('Uncleanable polygons', len(src2))

            
def rasterize_inventory(basenames, shp_path, tif_path, hdt_path, theme_cols, age_col, period_length, base_year,
                        cap_age=None, d=100., verbose=True):
    hdt = {}
    for bn in basenames:
        kwargs = {'shp_path':'%s/stands.shp' % shp_path(bn), 
                  'tif_path':'%s/inventory_init.tif' % tif_path(bn), 
                  'theme_cols':theme_cols, 
                  'age_col':age_col, 
                  'age_divisor':period_length,
                  'cap_age':cap_age,
                  'verbose':verbose,
                  'd':d}
        hdt[bn] = ws3.common.rasterize_stands(**kwargs)
        pickle.dump(hdt[bn], open('%s/hdt_%s.pkl' % (hdt_path, bn), 'wb'))
    return hdt


def compile_basecodes(hdt, basenames, theme_cols):
    import numpy as np
    bc1 = {bn:[list(np.unique(x)) for x in zip(*hdt[bn].values())] for bn in basenames}
    bc2 = [set() for _ in range(len(theme_cols))]
    for bn in basenames:
        for i in range(len(theme_cols)):
            bc2[i].update(bc1[bn][i])
    basecodes = [list(bc2[i]) for i in range(len(theme_cols))]
    return basecodes


def schedule_harvest_optimize(fm, basenames, scenario_name='base', util=0.85, param_funcs=None, 
                              target_path='./input/targets.csv', obj_mode='min_harea', mask=None):
    import gurobipy as grb
    p = gen_scen(fm, basenames, scenario_name, util, param_funcs=param_funcs, toffset=0, obj_mode=obj_mode, mask=mask, target_path=target_path)
    m = p.solve()
    if m.status != grb.GRB.OPTIMAL:
        print('Model not optimal.')
        return None
    sch = fm.compile_schedule(p)
    fm.reset_actions()
    fm.initialize_areas()
    fm.apply_schedule(sch, 
                      force_integral_area=True, 
                      override_operability=True,
                      fuzzy_age=True,
                      recourse_enabled=True,
                      verbose=False,
                      compile_c_ycomps=True)
    return sch


def schedule_harvest_areacontrol(fm, period=1, acode='harvest', util=0.85, 
                                 target_masks=None, target_areas=None, target_scalefactors=None, 
                                 verbose=False):
    print(target_masks, target_areas, target_scalefactors) # debug
    fm.reset_actions()
    if not target_areas:
        target_masks = ['? 1 ? ?'] if not target_masks else target_masks
        target_areas = []
        for i, mask in enumerate(target_masks): # compute area-weighted mean CMAI age for each masked DT set
            awr = []
            dtype_keys = fm.unmask(mask)
            for dtk in dtype_keys:
                dt = fm.dtypes[dtk]
                area = dt.area(0)
                cmai_age = dt.ycomp('totvol').mai().ytp().lookup(0)
                awr.append(area * cmai_age)
            r = sum(awr)  / fm.inventory(0, mask=mask)
            asf = 1. if not target_scalefactors else target_scalefactors[i]  
            ta = (1/r) * fm.inventory(0, mask=mask) * asf
            target_areas.append(ta)
    for mask, target_area in zip(target_masks, target_areas):
        fm.areaselector.operate(period, acode, target_area, mask=mask, verbose=verbose)
    sch = fm.compile_schedule()
    return sch


def sda(fm, basenames, time_step, tif_path, hdt, acode_map=None, nthresh=10, sda_mode='randblk'):
    from pathlib import Path
    from ws3.spatial import ForestRaster
    from ws3.common import hash_dt
    import os
    if acode_map is None:
        acode_map = {'harvest':'projected_harvest'}
    def cmp_fr_kwargs(bn):
        tmp_path = os.path.split(tif_path(bn))[0]
        _tif_path = '%s/%s' % (tmp_path, bn)
        if not Path(_tif_path).exists():
            Path(_tif_path).mkdir()
        fr_kwargs = {'hdt_map':hdt[bn], 
                     'hdt_func':hash_dt, 
                     'src_path':'%s/inventory_%i.tif' % (tif_path(bn), fm.base_year),
                     'snk_path':_tif_path,
                     'acode_map':acode_map,
                     'forestmodel':fm,
                     'horizon':fm.horizon,
                     'period_length':fm.period_length,
                     'time_step':time_step,
                     'base_year':fm.base_year,
                     'piggyback_acodes':{}}
        return fr_kwargs
    for bn in basenames:
        print('SDA for TSA', bn)
        mask = (bn, '?', '?', '?')
        fr = ForestRaster(**cmp_fr_kwargs(bn))
        fr.allocate_schedule(mask=mask, verbose=1, sda_mode=sda_mode, nthresh=nthresh)
        fr.cleanup()


def pickle_forestmodel(fm, scenario_name, basename):
    pickle.dump(fm, open('dat/out/%s_%s_fm.pkl' % (scenario_name, basename), 'wb'))

    
def pickle_schedule(sch, scenario_name, basename):
    pickle.dump(sch, open('dat/out/%s_%s_sch.pkl' % (scenario_name, basename), 'wb'))


def unpickle_forestmodel(scenario_name, basename):
    return pickle.load(open('dat/out/%s/%s_%s_fm.pkl' % (scenario_name, scenario_name, basename), 'rb'))


def unpickle_schedule(scenario_name, basename):
    return pickle.load(open('dat/out/%s/%s_%s_sch.pkl' % (scenario_name, scenario_name, basename), 'rb'))
