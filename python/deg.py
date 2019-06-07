import os
########################################################################################################
# gep_dev = True
# if gep_dev:
#     ws3_path = '../lib/ws3'
#     try:
#         _path = os.path.abspath(os.path.join(os.path.dirname(__file__), ws3_path))
#     except:
#         _path = os.path.abspath(os.path.join(os.path.dirname('__file__'), ws3_path))
#     if not _path in os.sys.path:
#         os.sys.path.insert(0, _path)
import ws3
from ws3.forest import ForestModel, Action
from ws3.spatial import ForestRaster
from ws3.common import clean_vector_data, reproject_vector_data, rasterize_stands, hash_dt, warp_raster
########################################################################################################
import rasterio
import numpy as np
import pandas as pd
try:
   import cPickle as pickle
except:
   import pickle


def read_basenames(path):
    return [line.lower().strip().split(' ')[0] 
        for line in open(path, 'r') if not line.startswith('#')]

def get_expr_c(expr, var_indices):
    for i in range(expr.size()):
        dvar = expr.getVar(i)
        yield expr.getCoeff(i), var_indices[dvar], dvar
        
def get_matrix_c(m, unhash_ij):
    dvars = m.getVars()
    constrs = m.getConstrs()
    var_indices = {v: i for i, v in enumerate(dvars)}
    for row_idx, constr in enumerate(constrs):
        #print constr, constr.slack
        for coeff, col_idx, dvar in get_expr_c(m.getRow(constr), var_indices):
            i, j = unhash_ij[dvar.VarName]
            yield row_idx, col_idx, constr.constrName, i, j, coeff, dvar.X

def cmp_c_z(fm, path, expr):
    result = 0.
    for t, n in enumerate(path, start=1):
        d = n.data()
        if fm.is_harvest(d['acode']):
            result += fm.compile_product(t, expr, d['acode'], [d['dtk']], d['age'], coeff=False)
    return result

def cmp_c_cflw(fm, path, expr, mask=None): # product, all harvest actions
    result = {}
    for t, n in enumerate(path, start=1):
        d = n.data()
        #if mask and d['dtk'] not in fm.unmask(mask): continue
        if mask and not fm.match_mask(mask, d['dtk']): continue
        if fm.is_harvest(d['acode']):
            result[t] = fm.compile_product(t, expr, d['acode'], [d['dtk']], d['age'], coeff=False)
    return result

def cmp_c_caa(fm, path, expr, acodes, mask=None): # product, named actions
    result = {}
    for t, n in enumerate(path, start=1):
        d = n.data()
        #if mask and d['dtk'] not in fm.unmask(mask): continue
        if mask and not fm.match_mask(mask, d['dtk']): continue
        if d['acode'] in acodes:
            result[t] = fm.compile_product(t, expr, d['acode'], [d['dtk']], d['age'], coeff=False)
    return result


def _gen_scen_base_maxhv(fm, basenames, name, util, harvest_acode='harvest', fire_acode='fire', 
                         tvy_name='totvol', toffset=0, target_path='dat/targets.csv', obj_mode='max_hvol'):
    from functools import partial
    acodes = ['null', harvest_acode]
    if fire_acode: acodes.append(fire_acode)
    df_targets = pd.read_csv(target_path).set_index(['tsa', 'period'])
    #cvcut = lambda bn, t: float(df_targets.loc[bn, t]['vcut']) if t < 5 else float(df_targets.loc[bn, 4]['vcut'])
    cabrn = lambda bn, t: float(df_targets.loc[bn, t]['abrn']) if t < 5 else float(df_targets.loc[bn, 4]['abrn'])
    cflw_acut_e = lambda bn, t: df_targets.loc[bn, t]['cflw_acut_e'] if t < 5 else df_targets.loc[bn, 4]['cflw_acut_e']
    cflw_vcut_e = lambda bn, t: df_targets.loc[bn, t]['cgen_vcut_e'] if t < 5 else df_targets.loc[bn, 4]['cgen_vcut_e']
    #cgen_vcut_e = lambda bn, t: df_targets.loc[bn, t]['cgen_vcut_e'] if t < 5 else df_targets.loc[bn, 4]['cgen_vcut_e']
    cgen_abrn_e = lambda bn, t: df_targets.loc[bn, t]['cgen_abrn_e'] if t < 5 else df_targets.loc[bn, 4]['cgen_abrn_e']
    vexpr = '%s * %0.2f' % (tvy_name, util)
    if obj_mode == 'max_hvol':
        sense = ws3.opt.SENSE_MAXIMIZE 
        zexpr = vexpr
    elif obj_mode == 'min_harea':
        sense = ws3.opt.SENSE_MINIMIZE 
        zexpr = '1.'
    else:
        raise ValueError('Invalid obj_mode: %s' % obj_mode)
    coeff_funcs = {'z':partial(cmp_c_z, expr=zexpr)}
    coeff_funcs.update({'cacut_%s'%bn:partial(cmp_c_caa, expr='1.', acodes=[harvest_acode], mask=(bn, '?', '?', '?')) 
                        for bn in basenames})
    coeff_funcs.update({'cvcut_%s'%bn:partial(cmp_c_caa, expr=vexpr, acodes=[harvest_acode], mask=(bn, '?', '?', '?')) 
                        for bn in basenames})
    if fire_acode:
        coeff_funcs.update({'cabrn_%s'%bn:partial(cmp_c_caa, expr='1.', acodes=[fire_acode], mask=(bn, '?', '?', '?')) 
                            for bn in basenames})
    T = fm.periods
    cflw_e, cgen_data = {}, {}
    cflw_ebn_vcut = {bn:({t:cflw_vcut_e(bn, t+toffset) for t in T}, 1) for bn in basenames}
    cflw_ebn_acut = {bn:({t:cflw_acut_e(bn, t+toffset) for t in T}, 1) for bn in basenames}
    cflw_e.update({'cvcut_%s'%bn:cflw_ebn_vcut[bn] for bn in basenames})
    cflw_e.update({'cacut_%s'%bn:cflw_ebn_acut[bn] for bn in basenames})
    for bn in basenames:
        #cgen_data.update({'cvcut_%s'%bn:{'lb':{t:cvcut(bn, t+toffset)*(1.-cgen_vcut_e(bn, t+toffset)) for t in T}, 
        #                                 'ub':{t:cvcut(bn, t+toffset) for t in T}}})
        if fire_acode:
            cgen_data.update({'cabrn_%s'%bn:{'lb':{t:cabrn(bn, t+toffset)*(1.-cgen_abrn_e(bn, t+toffset)) for t in T}, 
                                             'ub':{t:cabrn(bn, t+toffset) for t in T}}})
    return fm.add_problem(name, coeff_funcs, cflw_e, cgen_data=cgen_data, acodes=acodes, sense=sense)



def _gen_scen_base(fm, basenames, name, util, harvest_acode='harvest', fire_acode='fire', 
                   tvy_name='totvol', toffset=0, target_path='dat/targets.csv', obj_mode='max_hvol',
                   max_tp=6, cacut=None):
    from functools import partial
    acodes = ['null', harvest_acode]
    if fire_acode: acodes.append(fire_acode)
    df_targets = pd.read_csv(target_path).set_index(['tsa', 'period'])
    cvcut = lambda bn, t: float(df_targets.loc[bn, t]['vcut']) if t <= max_tp else float(df_targets.loc[bn, max_tp]['vcut'])
    cabrn = lambda bn, t: float(df_targets.loc[bn, t]['abrn']) if t <= max_tp else float(df_targets.loc[bn, max_tp]['abrn'])
    cflw_acut_e = lambda bn, t: df_targets.loc[bn, t]['cflw_acut_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cflw_acut_e']
    cgen_vcut_e = lambda bn, t: df_targets.loc[bn, t]['cgen_vcut_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cgen_vcut_e']
    cgen_acut_e = lambda bn, t: df_targets.loc[bn, t]['cgen_vcut_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cgen_vcut_e']
    cgen_abrn_e = lambda bn, t: df_targets.loc[bn, t]['cgen_abrn_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cgen_abrn_e']
    vexpr = '%s * %0.2f' % (tvy_name, util)
    if obj_mode == 'max_hvol':
        sense = ws3.opt.SENSE_MAXIMIZE 
        zexpr = vexpr
    elif obj_mode == 'min_harea':
        sense = ws3.opt.SENSE_MINIMIZE 
        zexpr = '1.'
    else:
        raise ValueError('Invalid obj_mode: %s' % obj_mode)
    coeff_funcs = {'z':partial(cmp_c_z, expr=zexpr)}
    coeff_funcs.update({'cacut_%s'%bn:partial(cmp_c_caa, expr='1.', acodes=[harvest_acode], mask=(bn, '?', '?', '?')) 
                        for bn in basenames})
    coeff_funcs.update({'cvcut_%s'%bn:partial(cmp_c_caa, expr=vexpr, acodes=[harvest_acode], mask=(bn, '?', '?', '?')) 
                        for bn in basenames})
    if fire_acode:
        coeff_funcs.update({'cabrn_%s'%bn:partial(cmp_c_caa, expr='1.', acodes=[fire_acode], mask=(bn, '?', '?', '?')) 
                            for bn in basenames})
    T = fm.periods
    cflw_e, cgen_data = {}, {}
    cflw_ebn = {bn:({t:cflw_acut_e(bn, t+toffset) for t in T}, 1) for bn in basenames}
    cflw_e.update({'cacut_%s'%bn:cflw_ebn[bn] for bn in basenames})
    for bn in basenames:
        cgen_data.update({'cvcut_%s'%bn:{'lb':{t:cvcut(bn, t+toffset)*(1.-cgen_vcut_e(bn, t+toffset)) for t in T}, 
                                         'ub':{t:cvcut(bn, t+toffset) for t in T}}})
        if cacut:
            cgen_data.update({'cacut_%s'%bn:{'lb':{t:cacut(bn, t)*(1. - cgen_acut_e(bn, t)) for t in T}, 
                                             'ub':{t:cacut(bn, t) for t in T}}})
        if fire_acode:
            cgen_data.update({'cabrn_%s'%bn:{'lb':{t:cabrn(bn, t+toffset)*(1.-cgen_abrn_e(bn, t+toffset)) for t in T}, 
                                             'ub':{t:cabrn(bn, t+toffset) for t in T}}})
        #if fire_acode:
        #    cgen_data.update({'cabrn_%s'%bn:{'lb':{t:cabrn(bn, t+toffset)*(1.-cgen_abrn_e(bn, t+toffset)) for t in T}, 
        #                                     'ub':{t:cabrn(bn, t+toffset)*(1.+cgen_abrn_e(bn, t+toffset)) for t in T}}})
    return fm.add_problem(name, coeff_funcs, cflw_e, cgen_data=cgen_data, acodes=acodes, sense=sense)


def gen_scen(fm, basenames, name, util, fire_acode, toffset=0, obj_mode='max_hvol', cacut=None):
    dsp = {'base':_gen_scen_base,
           'more':_gen_scen_base,
           'rogh':_gen_scen_base,
           'base_maxhv':_gen_scen_base_maxhv}
    #print(name)
    return dsp[name](fm, basenames, name, util, fire_acode=fire_acode, toffset=toffset,
                     obj_mode=obj_mode, cacut=cacut)

def unhash_ij(problem):
    r = {}
    for i, tree in problem.trees.items():
        for path in tree.paths():
            j = tuple(n.data('acode') for n in path)
            r['x_%i' % hash((i, j))] = i, j
    return r

def read_pslashburn(path):
    result = {}
    for line in open(path, 'r'):
        if line.strip().startswith('#'): continue
        k, v = line.lower().strip().split(',')
        result[k] = float(v)
    return result

def bootstrap_themes(fm, theme_cols=['theme0', 'theme1', 'theme2', 'theme3'], 
                     basecodes=[[], [], [], []], aggs=[{}, {}, {}, {}], verbose=False):
    for ti, t in enumerate(theme_cols):
        fm.add_theme(t, basecodes=basecodes[ti], aggs=aggs[ti])
    fm.nthemes = len(theme_cols)
        
def bootstrap_areas(fm, basenames, rst_path, hdt):
    for bn in basenames:
        with rasterio.open(rst_path(bn), 'r') as src:
            #print bn, src.meta
            pxa = pow(src.transform.a, 2) * 0.0001 # pixel area (hectares)
            bh, ba = src.read(1), src.read(2)
            for h, dt in hdt[bn].items():
                ra = ba[np.where(bh == h)] # match themes hash value
                fm.dtypes[dt] = ws3.forest.DevelopmentType(dt, fm)
                for age in np.unique(ra):
                    area = len(ra[np.where(ra == age)]) * pxa
                    fm.dtypes[dt].area(0, age, area)
            
def bootstrap_yields(fm, yld_path, theme_cols=['AU', 'LDSPP'], spcode='SPCode', 
                     startp_col='Wdks', x_max=360, y_cols=None, 
                     period_length=10, tvy_name='totvol'):
    y_cols = ['X%i' % i for i in range(0, x_max, period_length)]
    df = pd.read_csv(yld_path, usecols=theme_cols+[spcode, startp_col]+y_cols)
    df [theme_cols[1]] = df[theme_cols[1]].str.lower().str.replace(r'[- ]+', '_')
    for i in (0, 1): df[theme_cols[i]] = df[theme_cols[i]].astype(str)
    df = df.set_index(theme_cols)
    for t1, t2 in df.index.values: # assuming exactly one yield curve per unique combination of AU and LDSPP
        #print t1, t2
        mask = ('?', '?', t1, t2)
        dt_keys = fm.unmask(mask)
        if not dt_keys: continue
        #print t1, t2, len(dt_keys)
        r = df.loc[t1, t2]
        yname = str.lower(r[spcode])
        points = [(x+r[startp_col], r[y]) for x, y in enumerate(y_cols)]
        c = fm.register_curve(ws3.core.Curve(yname, points=points, type='a', is_volume=True, 
                                             xmax=fm.max_age, period_length=period_length))
        fm.yields.append((mask, 'a', [(yname, c)]))
        fm.ynames.add(yname)
        for dtk in dt_keys: fm.dtypes[dtk].add_ycomp('a', yname, c)
    # add total volume curve ###
    expr = '_SUM(%s)' % ', '.join(fm.ynames)
    #c_tv = wm.register_curve(ws3.core.Curve(tvy_name))
    fm.yields.append((('?', '?', '?', '?'), 'c', [(tvy_name, expr)]))
    fm.ynames.add(tvy_name)
    for dtk in fm.dtypes.keys(): fm.dtypes[dtk].add_ycomp('c', tvy_name, expr)

def bootstrap_ogi(fm, tvy_name='totvol', ra1_type='cmai', ra2_type='cyld', rc1=[1., 0.], rc2=[1., 0.], max_y=1.,
                  mask=None, yname='ogi', period_length=10):
    """
    Adds a yield component to each development type expressing "old-growthedness".
    f(x) = 0 on the interval [0, ra1*rc1[0]+rc1[1]].
    f(x) is linearly interpolated on the interval [ra2*rc2[0]+rc2[1], ra1*rc1[0]+rc1[1]]
    f(x) = 1 on the interval [ra2*rc2[0]+rc2[1], inf].
    ra1 defaults to age at which total volume MAI curve culminates.
    ra2 defaults to age at which total volume curve culminates.
    """
    mask = mask if mask else ('?', '?', '?', '?')
    fm.ynames.add(yname)
    for dtk in fm.unmask(mask):
        dt = fm.dtypes[dtk]
        yldca = dt.ycomp(tvy_name).ytp().lookup(0)
        maica = dt.ycomp(tvy_name).mai().ytp().lookup(0)
        ra1 = maica if ra1_type=='cmai' else yldca
        ra2 = maica if ra2_type=='cmai' else yldca
        points = [(0, 0.), 
                  (int(ra1*rc1[0]+rc1[1]), 0.), 
                  (int(ra2*rc2[0]+rc2[1]), max_y),
                  (fm.max_age, max_y)]
        #print(dtk, points)
        c = fm.register_curve(ws3.core.Curve(yname, points=points, type='a', is_volume=False, 
                                             xmax=fm.max_age, period_length=period_length))
        #print(dtk, c.points())
        #assert False
        _mask = (mask[0], '?', dtk[2], dtk[3])
        fm.yields.append((_mask, 'a', [(yname, c)]))
        dt.add_ycomp('a', yname, c)
            
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
            #print('_', acode, dtk)
            for age in range(1, fm.max_age):
                #print('__')
                if not dt.is_operable(acode, 1, age): continue
                fm.dtypes[dtk].transitions[acode, age] = target
                
def bootstrap_forestmodel(basenames, model_name, model_path, yld_path, tif_path, horizon, 
                          period_length, max_age, basecodes, action_params, hdt,
                          add_null_action=True, tvy_name='totvol', compile_actions=True, compile_ogi=False, verbose=0,
                          **kwargs):
    if verbose > 0: print('instantiating ForestModel object')
    fm = ForestModel(model_name=model_name, 
                     model_path=model_path, 
                     horizon=horizon,     
                     period_length=period_length,
                     max_age=max_age)
    if verbose > 0: print('bootstrap_themes')
    bootstrap_themes(fm, basecodes=basecodes)    
    if verbose > 0: print('bootstrap_areas')
    print(tif_path)
    bootstrap_areas(fm, basenames, tif_path, hdt)
    if verbose > 0: print('bootstrap_yields')
    bootstrap_yields(fm, yld_path, tvy_name=tvy_name)
    if compile_ogi:
        if verbose > 0: print('bootstrap_ogi')
        bootstrap_ogi(fm, ra2_type='cmai', rc2=[1.5, 0.])
    if verbose > 0: print('bootstrap_actions')
    bootstrap_actions(fm, action_params)
    if verbose > 0: print('adding null actions')
    if add_null_action: fm.add_null_action()
    if verbose > 0: print('compiling actions')
    fm.compile_actions()
    if verbose > 0: print('finalizing ForestModel config')
    fm.reset_actions()
    fm.initialize_areas()
    fm.grow()
    return fm


def clean_shapefiles():
    for bn in basenames:
        print('cleaning GDB', gdb_path(bn), shp_path(bn))
        if not Path(shp_path(bn)).exists(): 
            Path(shp_path(bn)).mkdir()
        snk1_path, snk2_path = clean_vector_data(gdb_path(bn), shp_path(bn), '_stands', prop_names, 
                                                 tolerance=tolerance, max_records=None, 
                                                 theme0=bn, prop_types=prop_types, dst_epsg=None)
        with fiona.open(gdb_path(bn)) as src0, fiona.open(snk1_path) as src1, fiona.open(snk2_path) as src2:
            print('Polygons in original dataset', len(src0))
            print('Polygons in clean dataset', len(src1))
            print('Uncleanable polygons', len(src2))
        reproject_vector_data(shp_path(bn)+'/_stands.shp', shp_path(bn)+'/stands.shp', snk_epsg)
        _path = '%s/shp/%s.shp' % (dat_path, bn)
        for f in [f for f in listdir(_path) if isfile(join(_path, f)) and f.startswith('_stands')]: 
            os.remove('%s/%s' % (_path, f))
            
def rasterize_inventory(basenames, shp_path, tif_path, hdt_path, theme_cols, age_col, period_length, cap_age=None, d=100., verbose=True):
    hdt = {}
    for bn in basenames:
        kwargs = {'shp_path':'%s/stands.shp' % shp_path(bn), 
                  'tif_path':tif_path(bn), 
                  'theme_cols':theme_cols, 
                  'age_col':age_col, 
                  'age_divisor':period_length,
                  'cap_age':cap_age,
                  'verbose':verbose,
                  'd':d}
        hdt[bn] = rasterize_stands(**kwargs)
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

def debug_single_pixel(fm, basenames, tif_path, hdt, horizon, period_length, base_year, p_slashburn):
    import copy
    from pathlib import Path
    from ws3.common import hash_dt

    fm.reset_actions()
    fm.initialize_areas()

    done = False
    dtype_keys = copy.copy(list(fm.dtypes.keys()))
    while not done: # Find operable developement type.
        dtk = dtype_keys.pop()
        dt = fm.dtypes[dtk]
        is_operable = dt.is_operable('harvest', 1)
        if is_operable:
            lo, hi = is_operable
            for age, area in dt._areas[1].items():
                if age >= lo and area >= 10.:
                    done = True
                    continue
        else:
            continue
    print('found match', dtk, age, fm.dtypes[dtk].area(1, age))

    fm.apply_action(dtk, 'harvest', 1, age, area=10.)

    debug_this = 0
    da = 0.
    fudge = 1.0
    acode_map = {'harvest':'projected_harvest', 
                 'fire':'projected_fire', 
                 'slashburn': 'projected_slashburn'}
    for bn in basenames[:]:
        print('processing', bn)
        mask = (bn, '?', '?', '?')
        tmp_path = os.path.split(tif_path(bn))[0]
        _tif_path = '%s/%s' % (tmp_path, bn)
        if not Path(_tif_path).exists():
            Path(_tif_path).mkdir()
        kwargs = {'hdt_map':hdt[bn], 
                  'hdt_func':hash_dt, 
                  'src_path':tif_path(bn),
                  'snk_path':_tif_path,
                  'acode_map':acode_map,
                  'forestmodel':fm,
                  'horizon':horizon,
                  'period_length':period_length,
                  'base_year':base_year,
                  'piggyback_acodes':{'harvest':[('slashburn', p_slashburn(bn))]}}
        if debug_this:
            fr = ForestRaster(**kwargs)
            fr.allocate_schedule(mask=mask, da=da, fudge=fudge, verbose=True, sda_mode='randpxl')
            fr.cleanup()
        else:
            with ForestRaster(**kwargs) as fr:
                fr.allocate_schedule(mask=mask, da=da, fudge=fudge, verbose=True, sda_mode='randpxl')
    #with rasterio.open('dat/tif/tsa10/projected_harvest_2015.tif') as src:
    #    x = src.read()
    #np.unique(x), x.sum()

    

def roll_horizon(fm):
    fm.grow(1, False)
    clobber_areas(fm, 2)
    fm.initialize_areas()
    fm.grow()
    fm.reset_actions()
    
    
def clobber_areas(fm, from_period):
    import copy
    for dt in list(fm.dtypes.values()):
        dt._areas[0] = copy.copy(dt._areas[from_period])            
    
    
def copy_areas(fm, period=0):
    import copy
    return {dtk:copy.deepcopy(fm.dtypes[dtk]._areas[period]) for dtk in list(fm.dtypes.keys())}


def reset_areas(fm, areas):
    for dt in list(fm.dtypes.values()):
        dt._areas[0] = areas[dt.key]     

        
def compile_schedule_rhrs(fm, sch):
    result = []
    for irh in range(len(sch)):
        for dtype_key, age, area, acode, period, etype in sch[irh]:
            if period > 1: break
            result.append((dtype_key, age, area, acode, irh+1, etype))
    return result


def burn_aspatial(fm, basenames, period, burn_area, acode='fire'):
    for bn in basenames:
        mask = (bn, '?', '?', '?')
        p_burn = burn_area[bn] / fm.operable_area(acode, period, mask=mask)
        #print('burn_aspatial', bn, burn_area[bn], fm.operable_area(acode, period, mask=mask), p_burn)
        for dtk in fm.unmask(mask):
            dt = fm.dtypes[dtk]
            ages = list(dt._areas[period].keys())
            for age in ages:
                if dt.is_operable(acode, period, age):
                    #print('apply_action', dtk, acode, period, age, dt.area(period, age) * p_burn)
                    fm.apply_action(dtk, acode, period, age, dt.area(period, age) * p_burn)

                    
    
def sda(fm, basenames, base_year, period_length, horizon, tif_path, hdt, acode_map=None,
        p_slashburn_path='dat/p_slashburn.csv', nthresh=10, sda_mode='randblk',
        **kwargs):
    from pathlib import Path
    _p_slashburn = read_pslashburn(p_slashburn_path)
    p_slashburn = lambda bn: _p_slashburn[bn]
    if acode_map is None:
        acode_map = {'harvest':'projected_harvest', 
                     'fire':'projected_fire', 
                     'slashburn': 'projected_slashburn'}
    def cmp_fr_kwargs(bn):
        tmp_path = os.path.split(tif_path(bn))[0]
        _tif_path = '%s/%s' % (tmp_path, bn)
        if not Path(_tif_path).exists():
            Path(_tif_path).mkdir()
        fr_kwargs = {'hdt_map':hdt[bn], 
                     'hdt_func':hash_dt, 
                     'src_path':tif_path(bn),
                     'snk_path':_tif_path,
                     'acode_map':acode_map,
                     'forestmodel':fm,
                     'horizon':horizon,
                     'period_length':period_length,
                     'base_year':base_year,
                     'piggyback_acodes':{'harvest':[('slashburn', p_slashburn(bn))]}}
        return fr_kwargs
    for bn in basenames:
        print('SDA for TSA', bn)
        mask = (bn, '?', '?', '?')
        fr = ForestRaster(**cmp_fr_kwargs(bn))
        fr.allocate_schedule(mask=mask, verbose=1, sda_mode=sda_mode, nthresh=nthresh)
        fr.cleanup()
        #with ForestRaster(**cmp_fr_kwargs(bn)) as fr:
        #    fr.allocate_schedule(mask=mask, verbose=True, sda_mode='randpxl')

def run_deg_classic(fm, basenames, scenario_name, util, tif_path, hdt, period_length=10, horizon=6, 
                    base_year=2015, target_path='dat/targets.csv', sda_mode='randblk', obj_mode='max_hvol', fire_acode='fire'):
    import gurobipy as grb
    p = gen_scen(fm, basenames, scenario_name, util, fire_acode=fire_acode, toffset=0, obj_mode=obj_mode)
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
 
                    
def run_deg_rhrs(fm, basenames, scenario_name, util, tif_path, hdt, period_length=10, horizon=6, fm_horizon=2,
                 base_year=2015, target_path='dat/targets.csv', sda_mode='randblk', obj_mode='max_hvol',
                 pcut=1.0, pbrn=1.0, fm_lookahead=2, max_tp=6, twopass_cacut_factor=None, harvest_acode='harvest'):
    import gurobipy as grb
    import random
    fm.set_horizon(fm_horizon)
    fm.initialize_areas()
    fm.grow()             
    fm.reset_actions()
    areas0 = copy_areas(fm)
    df_targets = pd.read_csv(target_path).set_index(['tsa', 'period'])
    cabrn = lambda bn, t: float(df_targets.loc[bn, t]['abrn']) if t <= max_tp else float(df_targets.loc[bn, max_tp]['abrn'])
    sch = []
    for irh in range(horizon):
        print('Running iteration', irh)
        if fm_horizon > (horizon - irh):
            fm_horizon = horizon + fm_lookahead - irh - 1
            print(irh, 'truncating fm_horizon to', fm_horizon)
            fm.set_horizon(fm_horizon)
        burn_area = {bn:cabrn(bn, irh+1)*pbrn for bn in basenames}
        burn_aspatial(fm, basenames, 1, burn_area)
        sch_brn = fm.compile_schedule()
        _areas = copy_areas(fm)
        clobber_areas(fm, 1)
        p = gen_scen(fm, basenames, scenario_name, util, fire_acode=None, toffset=irh, obj_mode=obj_mode)
        m = p.solve()
        if m.status != grb.GRB.OPTIMAL:
            print('Model not optimal.')
            return None
        sch_cut = fm.compile_schedule(p)
        if twopass_cacut_factor:
            # Run the model a second time, this time activating harvest area constraints set to a user-defined
            # proportion of periodic harvest area in the base run. The point of this is to force simulation of
            # a harvest volume:area ratio that does not tend to occur when obj_mode is set to one of the two
            # available modes (i.e., 'max_hvol' or 'min_harea').
            fm.reset_actions()
            fm.apply_schedule(sch_cut, 
                              fuzzy_age=True,
                              recourse_enabled=True,
                              verbose=False,
                              scale_area=pcut)
            _cacut = {bn:[fm.compile_product(t, '1.', harvest_acode, fm.unmask('%s ? ? ?' % bn))*twopass_cacut_factor
                          for t in fm.periods]
                      for bn in basenames}
            print('activating cacut constraints', fm.horizon, fm.periods, _cacut)
            cacut = lambda bn, t: _cacut[bn][t-1]
            p = gen_scen(fm, basenames, scenario_name, util, fire_acode=None, toffset=irh, obj_mode=obj_mode, cacut=cacut)
            m = p.solve()
            if m.status != grb.GRB.OPTIMAL:
                print('Model not optimal.')
                return None
            sch_cut = fm.compile_schedule(p)
        reset_areas(fm, _areas)
        fm.reset_actions()
        fm.apply_schedule(sch_brn, 
                          max_period=1,
                          force_integral_area=False, 
                          override_operability=False,
                          fuzzy_age=False,
                          recourse_enabled=False,
                          verbose=False,
                          compile_c_ycomps=True)
        fm.apply_schedule(sch_cut, 
                          max_period=1,
                          force_integral_area=True, 
                          override_operability=True,
                          fuzzy_age=True,
                          recourse_enabled=True,
                          verbose=False,
                          compile_c_ycomps=True,
                          scale_area=pcut)
        sch.append(fm.compile_schedule())
        roll_horizon(fm)
    fm.set_horizon(horizon)
    reset_areas(fm, areas0)
    fm.reset_actions()
    fm.initialize_areas()
    sch_ = compile_schedule_rhrs(fm, sch)
    fm.apply_schedule(sch_, 
                      force_integral_area=True, 
                      override_operability=True,
                      fuzzy_age=True,
                      recourse_enabled=True,
                      verbose=False,
                      compile_c_ycomps=True)
    return sch_

def run_deg(fm, basenames, scenario_name, util, tif_path, hdt, horizon, deg_mode='classic', fm_horizon=2, sda_mode='randblk',
            obj_mode='min_harea', pcut=1.0, pbrn=1.0, twopass_cacut_factor=None, fire_acode='fire',
            **kwargs):
    if deg_mode == 'rhrs':
        sch = run_deg_rhrs(fm, basenames, scenario_name, util, tif_path, hdt, horizon=horizon, fm_horizon=fm_horizon,
                           sda_mode=sda_mode, obj_mode=obj_mode, pcut=pcut, pbrn=pbrn, twopass_cacut_factor=twopass_cacut_factor)
    elif deg_mode == 'classic':
        sch = run_deg_classic(fm, basenames, scenario_name, util, tif_path, hdt, horizon, sda_mode=sda_mode, obj_mode=obj_mode, fire_acode=fire_acode)
    return sch
    
def zip_output(basenames, scenario_name):
    import shutil
    for bn in basenames:
        print(shutil.make_archive('dat/out/%s_%s_distout' % (scenario_name, bn), 'zip', 'dat/tif', bn, True))
        
def pickle_forestmodel(fm, scenario_name, basename):
    pickle.dump(fm, open('dat/out/%s_%s_fm.pkl' % (scenario_name, basename), 'wb'))

    
def pickle_schedule(sch, scenario_name, basename):
    pickle.dump(sch, open('dat/out/%s_%s_sch.pkl' % (scenario_name, basename), 'wb'))

def compile_schedule_from_geotiffs(sn, bn, hdt):
    print('compiling area report data for: scenario %s, %s' % (scenario_name, bn))
    with zipfile.ZipFile('dat/out/%s/%s_%s_distout.zip' % (sn, sn, bn), 'r') as zf:
        zf.extractall('dat/out/%s' % sn)
    for dn in ['fire', 'harvest']:
        for dy in range(60):
            yr = 2015 + dy
            with rasterio.open('dat/out/%s/%s/projected_%s_%i.tif' % (sn, bn, dn, yr)) as src:
                n = src.read()[0]

def do_the_thing(basenames, model_name, model_path, yld_path, tif_path, shp_path, hdt_path, theme_cols, age_col, horizon, period_length,
                 max_age, basecodes, action_params, hdt, add_null_action, tvy_name, compile_actions,
                 compile_ogi, verbose, scenario_name, util, deg_mode, fm_horizon, sda_mode, obj_mode,
                 base_year, fm=None, run_sda=True, cmp_fm_only=False, rasterize_inv=True, cap_age=90, 
                 twopass_cacut_factor=None):
    from deg import bootstrap_forestmodel, pickle_forestmodel, pickle_schedule, rasterize_inventory
    if rasterize_inv: rasterize_inventory(basenames, shp_path, tif_path, hdt_path, theme_cols, age_col, period_length, cap_age=cap_age)
    try:
        fm = bootstrap_forestmodel(basenames=basenames, 
                                   model_name=model_name, 
                                   model_path=model_path, 
                                   yld_path=yld_path, 
                                   tif_path=tif_path, 
                                   horizon=horizon,
                                   period_length=period_length, 
                                   max_age=max_age, 
                                   basecodes=basecodes, 
                                   action_params=action_params, 
                                   hdt=hdt, 
                                   add_null_action=add_null_action,
                                   tvy_name=tvy_name, 
                                   compile_actions=compile_actions, 
                                   compile_ogi=compile_ogi, 
                                   verbose=verbose)
        pickle_forestmodel(fm, scenario_name, basenames[0])
        if cmp_fm_only:
            return 0
        #print('foo bar')
        sch = run_deg(fm=fm, 
                      basenames=basenames, 
                      scenario_name=scenario_name, 
                      util=util, 
                      tif_path=tif_path, 
                      hdt=hdt, 
                      horizon=horizon, 
                      deg_mode=deg_mode, 
                      fm_horizon=fm_horizon, 
                      sda_mode=sda_mode, 
                      obj_mode=obj_mode,
                      twopass_cacut_factor=twopass_cacut_factor)
        pickle_schedule(sch, scenario_name, basenames[0])
        if run_sda:
            sda(fm=fm, 
                basenames=basenames, 
                base_year=base_year, 
                period_length=period_length, 
                horizon=horizon, 
                tif_path=tif_path, 
                hdt=hdt, 
                sda_mode=sda_mode)
        zip_output(basenames, scenario_name)
        return 0
    except Exception as e:
        return e

def _do_the_thing(basenames, model_name, model_path, yld_path, tif_path, shp_path, hdt_path, horizon, period_length,
                  max_age, basecodes, action_params, hdt, add_null_action, tvy_name, compile_actions,
                  compile_ogi, verbose, scenario_name, util, deg_mode, fm_horizon, sda_mode, obj_mode,
                  base_year, fm=None, run_sda=True, cmp_fm_only=False, rasterize_inv=True, cap_age=90, 
                  twopass_cacut_factor=None):
    from deg import bootstrap_forestmodel, pickle_forestmodel, pickle_schedule, rasterize_inventory
    if rasterize_inv: rasterize_inventory(basenames, shp_path, tif_path, hdt_path, theme_cols, age_col, period_length, cap_age=cap_age)
    fm = bootstrap_forestmodel(basenames=basenames, 
                               model_name=model_name, 
                               model_path=model_path, 
                               yld_path=yld_path, 
                               tif_path=tif_path, 
                               horizon=horizon,
                               period_length=period_length, 
                               max_age=max_age, 
                               basecodes=basecodes, 
                               action_params=action_params, 
                               hdt=hdt, 
                               add_null_action=add_null_action,
                               tvy_name=tvy_name, 
                               compile_actions=compile_actions, 
                               compile_ogi=compile_ogi, 
                               verbose=verbose)
    pickle_forestmodel(fm, scenario_name, basenames[0])
    if cmp_fm_only:
        return 0
    #print('foo bar')
    sch = run_deg(fm=fm, 
                  basenames=basenames, 
                  scenario_name=scenario_name, 
                  util=util, 
                  tif_path=tif_path, 
                  hdt=hdt, 
                  horizon=horizon, 
                  deg_mode=deg_mode, 
                  fm_horizon=fm_horizon, 
                  sda_mode=sda_mode, 
                  obj_mode=obj_mode,
                  twopass_cacut_factor=twopass_cacut_factor)
    pickle_schedule(sch, scenario_name, basenames[0])
    if run_sda:
        sda(fm=fm, 
            basenames=basenames, 
            base_year=base_year, 
            period_length=period_length, 
            horizon=horizon, 
            tif_path=tif_path, 
            hdt=hdt, 
            sda_mode=sda_mode)
    zip_output(basenames, scenario_name)
    
def deg_diagnostic1(fm, scenario_name, basenames, sch, util=0.85, figsize=(10, 4)):
    import matplotlib.pyplot as plt
    from pandas import DataFrame as DF
    exprs = {'harvest':'totvol * %0.2f' % util, 'fire':'1.'}
    _acodes = 'harvest', 'fire'
    for bn in basenames:
        fig, ax = plt.subplots(figsize=figsize)
        mask = (bn, '?', '?', '?')
        dtype_keys = fm.unmask(mask)
        _area = {_acode:{t:0.  for t in fm.periods} for _acode in _acodes}
        for dtk, age, area, acode, period, etype in sch:
            if dtk in dtype_keys and acode in _acodes:
                _area[acode][period] += area
        d = {'t':[], 'harvest':[], 'fire':[], 
             'acut_ws3':[], 
             'acut_sch':list(_area['harvest'].values())}
        for t in fm.periods:
            d['t'].append(t)
            d['acut_ws3'].append(fm.compile_product(t, '1.', 'harvest', dtype_keys=dtype_keys))
            for acode in _acodes:
                expr = exprs[acode]
                d[acode].append(fm.compile_product(t, expr, acode, dtype_keys=dtype_keys))
        #print(d)
        df = DF(data=d).set_index('t')
        df['unit_hvol'] = df['harvest'] / df['acut_ws3']
        df = df.astype(int)
        df[['acut_ws3', 'acut_sch']].plot.bar(ax=ax)
        print(df)
        print(df.unit_hvol.mean())
        df.to_csv('dat/out/%s/%s_%s_diag1.csv' % (scenario_name, scenario_name, bn))
        fig.savefig('dat/out/%s/%s_%s_diag1.png' % (scenario_name, scenario_name, bn), format='png')

        
def deg_diagnostic2(fm, scenario_name, basenames, sch, tif_path, figsize=(10, 10), D=('harvest', 'fire'), period_length=10):
    import matplotlib.pyplot as plt
    from pandas import DataFrame as DF
    fig, ax = plt.subplots(3, 1, figsize=figsize)
    for bn in basenames:
        _sch = [x for x in sch if x[0][0] == bn]
        df = pd.DataFrame(_sch, columns=['themevals', 'age', 'area', 'action', 'period', 'foo'])
        #df = df.query('period <= 6')
        df['period'] = df['period']*period_length - period_length + 1
        df['area'] /= float(period_length)
        sch_area = df.groupby(by=['action', 'period'])['area'].sum()
        sda_area = {d:[0] for d in D}
        for i, d in enumerate(D):
            for j in range(fm.horizon*period_length):
                path = '%s/projected_%s_%i.tif' % (tif_path(bn)[:-4], d, 2015+j)
                with rasterio.open(path) as src:
                    sda_area[d].append(src.read().sum())
        for i, d in enumerate(D):
            if d in sch_area.index:
                sch_area.loc[d].plot(style='o', ax=ax[i], label='Targeted (%s, %s)' % (bn, d))
            ax[i].plot(sda_area[d], '*', label='Scheduled (%s, %s)' % (bn, d))
            ax[i].legend()
            ax[i].set_title('Validation of targeted and scheduled area')
            plt.xlabel('Period (years)')
            plt.ylabel('Area (hectares)')
        for d in D:
            if d in sch_area.index:
                T = [x for x in sch_area.loc[d] for i in range(fm.period_length)]
            S = sda_area[d][1:]
            E = [1.-(s/t) for t, s in zip(T, S)]
            plt.plot(E, 'o', label=d)
        ax[2].set_title('Relative target shortfall')
        ax[2].set_xlabel('Period (years)')
        ax[2].set_ylabel('Proportion of target not scheduled')
        ax[2].legend()
        fig.savefig('dat/out/%s/%s_%s_diag2.png' % (scenario_name, scenario_name, bn), format='png')
    
    
def deg_diagnostic3(fm, scenario_name, basenames, distname, figsize=(15, 6)):
    import matplotlib.pyplot as plt
    from pandas import DataFrame as DF
    n = 4
    fig, ax = plt.subplots(1, 4, figsize=figsize, sharex=True, sharey=True)
    for ni, bn in enumerate(basenames):
        for j in [0, 1]:
            ax[j].set_yticklabels([])
            ax[j].set_xticklabels([])
            ax[j].set_title('%s %s (%i)' % (bn, distname, 2015+j))
            path = '%s/%s/projected_%s_%i.tif' % (os.path.split(tif_path(bn))[0], bn, distname, 2015+j)
            rasterio.plot.show(rasterio.open(path), ax=ax[j])
        for j in [2, 3]:
            ax[j].set_yticklabels([])
            ax[j].set_xticklabels([])
            year = 2015+horizon-3+(j-1)
            ax[j].set_title('%s %s (%i)' % (bn, distname, year))
            path = '%s/%s/projected_%s_%i.tif' % (os.path.split(tif_path(bn))[0], bn, distname, year)
            rasterio.plot.show(rasterio.open(path), ax=ax[j])
    plt.tight_layout()
    
def unpickle_forestmodel(scenario_name, basename):
    return pickle.load(open('dat/out/%s/%s_%s_fm.pkl' % (scenario_name, scenario_name, basename), 'rb'))

def unpickle_schedule(scenario_name, basename):
    return pickle.load(open('dat/out/%s/%s_%s_sch.pkl' % (scenario_name, scenario_name, basename), 'rb'))
