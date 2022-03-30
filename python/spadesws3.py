import math
import rasterio
import numpy as np
import pandas as pd
import shutil
try:
    import cPickle as pickle
except:
    import pickle
import ws3
import re


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
            result += fm.compile_product(t, expr, acode=d['acode'], dtype_keys=[d['dtk']], age=d['age'], coeff=False)
    return result


def cmp_c_zbird(fm, path):
    """
    Compile bird AF objective function coefficient (given ForestModel instance, 
    leaf-to-root-node path).
    """
    result = 0.
    for t, n in enumerate(path, start=1):
        d = n.data()
        age = d['age']
        if d['acode'] != 'null':
            _, tprop, tyield, tage, _, _, _ = fm.dt(d['dtk']).transitions[d['acode'], d['age']][0]
            assert tprop == 1.
            age = fm.resolve_targetage(d['dtk'], tyield, d['age'], tage, d['acode'])
        result += fm.dt(d['dtk']).ycomp('birdaf')[age] * d['area']
        #result += fm.inventory(t, yname='birdaf', age=n.data('age'), dtype_keys=[n.data()['dtk']])
    return result


def _cmp_c_cbird(fm, path, mask=None):
    """
    Compile bird AF coefficients (given ForestModel instance, leaf-to-root-node path, mask).
    """
    result = {}
    for t, n in enumerate(path, start=1):
        d = n.data()
        if mask and not fm.match_mask(mask, d['dtk']): continue
        age = d['age']
        if d['acode'] != 'null':
            _, tprop, tyield, tage, _, _, _ = fm.dt(d['dtk']).transitions[d['acode'], d['age']][0]
            assert tprop == 1.
            age = fm.resolve_targetage(d['dtk'], tyield, d['age'], tage, d['acode'])
        result[t] = fm.inventory(t, yname='birdaf', age=age, dtype_keys=[d['dtk']])
    return result

def cmp_c_cbird(fm, path, mask=None):
    """
    Compile bird AF coefficients (given ForestModel instance, leaf-to-root-node path, mask).
    """
    result = {}
    for t, n in enumerate(path, start=1):
        d = n.data()
        if mask and not fm.match_mask(mask, d['dtk']): continue
        age = d['age']
        if d['acode'] != 'null':
            _, tprop, tyield, tage, _, _, _ = fm.dt(d['dtk']).transitions[d['acode'], d['age']][0]
            assert tprop == 1.
            age = fm.resolve_targetage(d['dtk'], tyield, d['age'], tage, d['acode'])
        result[t] = fm.dt(d['dtk']).ycomp('birdaf')[age] * d['area']
    return result

    
# no good
#def cmp_c_zbird(fm, path):
#    """
#    Compile bird AF objective function coefficient (given ForestModel instance, 
#    leaf-to-root-node path).
#    """
#    result = fm.inventory(len(path), yname='birdaf', dtype_keys=[path[-1].data()['dtk']])
#    return result


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
            result[t] = fm.compile_product(t, expr, acode=d['acode'], dtype_keys=[d['dtk']], age=d['age'], coeff=False)
    return result


def cmp_c_caa(fm, path, expr, acodes, mask=None): # product, named actions
    """
    Compile actioned area constraint coefficient (given ForestModel instance, 
    leaf-to-root-node path, expression to evaluate, list of action codes, and optional mask).
    """
    result = {}
    for t, n in enumerate(path, start=1):
        d = n.data()
        if mask and not fm.match_mask(mask, d['dtk']): continue
        if d['acode'] in acodes:
            result[t] = fm.compile_product(t, expr, acode=d['acode'], dtype_keys=[d['dtk']], age=d['age'], coeff=False)
    return result
    
    
def cmp_c_cia(fm, path, yname, mask=None): # product, named actions
    """
    Compile inventory area constraint coefficient (given ForestModel instance, 
    leaf-to-root-node path, yield component name, and optional mask).
    """
    result = {}
    for t, n in enumerate(path, start=1):
        d = n.data()
        if mask and not fm.match_mask(mask, d['dtk']): continue
        result[t] = fm.inventory(t, yname, dtype_keys=[d['dtk']])
    return result
    

def _gen_scen_base(fm, basenames, name='base', util=0.85, param_funcs=None, target_scalefactors=None, harvest_acode='harvest', fire_acode='fire', 
                   tvy_name='totvol', toffset=0, obj_mode='min_harea', target_path='./input/targets.csv',
                   max_tp=2020, cacut=None, mask=None):
    fm.foo3 = target_scalefactors
    from functools import partial
    acodes = ['null', harvest_acode, fire_acode]  
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
        target_scalefactors = {bn:1. for bn in basenames} if not target_scalefactors else target_scalefactors
        df_targets = pd.read_csv(target_path).set_index(['tsa', 'year'])
        param_funcs = {}
        param_funcs['cvcut'] = lambda bn, t: float(df_targets.loc[bn, t]['vcut']) * target_scalefactors[bn] if t <= max_tp else float(df_targets.loc[bn, max_tp]['vcut']) * target_scalefactors[bn]
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
    if fire_acode:
        for i in ['0', '1']:
            coeff_funcs.update({'cabrn-thlb%s_%s' % (i, bn):partial(cmp_c_caa, expr='1.', acodes=[fire_acode], mask=(bn, i, '?', '?')) 
                                for bn in basenames})
    T = fm.periods# [fm.base_year+(t-1)*fm.period_length for t in fm.periods]
    cflw_e, cgen_data = {}, {}
    #foo = {bn:{t:(bn, t+toffset) for t in T} for bn in basenames}
    #print(T)
    #assert False
    cflw_ebn = {bn:({t:param_funcs['cflw_acut_e'](bn, fm.base_year+(t-1)*fm.period_length+toffset) for t in T}, fm.periods[-1]) for bn in basenames}
    cflw_e.update({'cacut_%s'%bn:cflw_ebn[bn] for bn in basenames})
    for bn in basenames:
        #print(df_targets.loc[bn])
        cgen_data.update({'cvcut_%s' % bn:{'lb':{t:param_funcs['cvcut'](bn, fm.base_year+(t-1)*fm.period_length+toffset) * fm.period_length *
                                                 (1. - param_funcs['cgen_vcut_e'](bn, fm.base_year+(t-1)*fm.period_length+toffset))
                                               for t in T}, 
                                         'ub':{t:param_funcs['cvcut'](bn, fm.base_year+(t-1)*fm.period_length+toffset) * fm.period_length for t in T}}})
        if cacut:
            cgen_data.update({'cacut_%s' % bn:{'lb':{t:param_funcs['cacut'](bn, fm.base_year+(t-1)*fm.period_length) * fm.period_length *
                                                     (1. - param_funcs['cgen_acut_e'](bn, fm.base_year+(t-1)*fm.period_length)) for t in T}, 
                                               'ub':{t:param_funcs['cacut'](bn, fm.base_year+(t-1)*fm.period_length) * fm.period_length for t in T}}})
        if fire_acode:
            for i in ['0', '1']:
                p = fm.inventory(0, mask='? %s ? ?' % i) / fm.inventory(0, mask='? ? ? ?')
                cgen_data.update({'cabrn-thlb%s_%s' % (i, bn):{'lb':{t:param_funcs['cabrn'](bn, fm.base_year+(t-1)*fm.period_length) * p * fm.period_length *
                                                                     (1. - param_funcs['cgen_abrn_e'](bn, fm.base_year+(t-1)*fm.period_length)) for t in T}, 
                                                               'ub':{t:param_funcs['cabrn'](bn, fm.base_year+(t-1)*fm.period_length) * p * fm.period_length for t in T}}})
                #fm.cgen_data = cgen_data
                #assert False
    #print(cflw_e)
    fm._tmp = {}
    fm._tmp['param_funcs'] = param_funcs
    fm._tmp['cgen_data'] = cgen_data
    return fm.add_problem(name, coeff_funcs, cflw_e, cgen_data=cgen_data, acodes=acodes, sense=sense, mask=mask)


def _gen_scen_maxharv_fire(fm, basenames, name='maxharv_fire', util=0.85, param_funcs=None, target_scalefactors=None,
                           harvest_acode='harvest', fire_acode='fire', 
                           tvy_name='totvol', toffset=0, obj_mode='max_hvol', target_path='./input/targets.csv',
                           max_tp=2020, cbird=None, mask=None, cacut=None):
    from functools import partial
    cvcut = cbird
    fm._tmp = {}
    T = fm.periods
    acodes = ['null', harvest_acode, fire_acode]  
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
        param_funcs = {}
        #assert target_scalefactors
        if target_scalefactors:
            print(target_scalefactors)
            # run maxharv scenario
            schedule_harvest_optimize(fm, basenames, 'maxharv_fire')
            ub = {(bn, t):fm.compile_product(t, vexpr, mask=(bn, '?', '?', '?')) for t in T for bn in basenames}
            # interpolate intermediate RHS vectors
            cvcut = {(bn, t):(ub[bn, t] * target_scalefactors[bn]) for t in T for bn in basenames}
            print(ub)
            print(cvcut)
            #print('foo1')
            param_funcs['cvcut'] = lambda bn, t: cvcut[(bn, t)] 
            #fm._tmp['cvcut'] = cvcut
            #fm._tmp['param_funs'] = param_funcs
            print('xxxxxxxxxxxxxxxxxxxxxxxxxx')
        #target_scalefactors = {bn:1. for bn in basenames} if not target_scalefactors else target_scalefactors
        df_targets = pd.read_csv(target_path).set_index(['tsa', 'year'])
        param_funcs['cabrn'] = lambda bn, t: float(df_targets.loc[bn, t]['abrn']) if t <= max_tp else float(df_targets.loc[bn, max_tp]['abrn'])
        param_funcs['cflw_acut_e'] = lambda bn, t: df_targets.loc[bn, t]['cflw_acut_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cflw_acut_e']
        param_funcs['cflw_vcut_e'] = lambda bn, t: df_targets.loc[bn, t]['cflw_vcut_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cflw_vcut_e']
        param_funcs['cgen_abrn_e'] = lambda bn, t: df_targets.loc[bn, t]['cgen_abrn_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cgen_abrn_e']
        param_funcs['cgen_vcut_e'] = lambda bn, t: df_targets.loc[bn, t]['cgen_vcut_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cgen_vcut_e']
        
        if cvcut:
            if isinstance(cvcut, float): # bootstrap constraint RHS values
                c = cvcut
                schedule_harvest_optimize(fm, basenames, 'maxharv_fire')
                ub = {(bn, t):fm.compile_product(t, expr=vexpr, mask=(bn, '?', '?', '?')) for t in T for bn in basenames}
                lb = {(bn, t):0. for t in T for bn in basenames}
                cvcut = {(bn, t):((ub[bn, t] - lb[bn, t]) * c) + lb[bn, t]  for t in T for bn in basenames}
                fm._tmp['cvcut'] = cvcut, c, ub, lb
            param_funcs['cvcut'] = lambda bn, t: cvcut[(bn, t)]  if (bn, t) in cvcut else 0. 
    coeff_funcs = {'z':partial(cmp_c_z, expr=zexpr)}
    coeff_funcs.update({'cacut_%s' % bn:partial(cmp_c_caa, expr='1.', acodes=[harvest_acode], mask=(bn, '?', '?', '?')) for bn in basenames})
    coeff_funcs.update({'cvcut_%s' % bn:partial(cmp_c_caa, expr=vexpr, acodes=[harvest_acode], mask=(bn, '?', '?', '?')) for bn in basenames})
    if fire_acode:
        for i in ['0', '1']:
            coeff_funcs.update({'cabrn-thlb%s_%s' % (i, bn):partial(cmp_c_caa, expr='1.', acodes=[fire_acode], mask=(bn, i, '?', '?')) 
                                for bn in basenames})
    cflw_e, cgen_data = {}, {}
    fm._tmp['foo'] = param_funcs, toffset, T, basenames, df_targets
    cflw_ebn = {bn:({t:param_funcs['cflw_acut_e'](bn, fm.base_year+(t-1)*fm.period_length+toffset) for t in T}, fm.periods[-1]) for bn in basenames}
    cflw_e.update({'cacut_%s'%bn:cflw_ebn[bn] for bn in basenames})
    cflw_ebn = {bn:({t:param_funcs['cflw_vcut_e'](bn, fm.base_year+(t-1)*fm.period_length+toffset) for t in T}, fm.periods[-1]) for bn in basenames}
    cflw_e.update({'cvcut_%s'%bn:cflw_ebn[bn] for bn in basenames})
    for bn in basenames:
        if target_scalefactors:
            fm._tmp['param_funs'] = param_funcs
            cgen_data.update({'cvcut_%s' % bn:{'lb':{t:0. for t in T},
                                               'ub':{t:param_funcs['cvcut'](bn, t) for t in T}}})
            print('cgen_data cvcut:', cgen_data['cvcut_%s' % bn])
        #if cvcut:
        #    cgen_data.update({'cvcut_%s' % bn:{'lb':{t:param_funcs['cvcut'](bn, t) for t in T}, 
        #                                       'ub':{t:np.inf for t in T}}})
        #else:
        #    cgen_data.update({'cvcut_%s' % bn:{'lb':{t:0. for t in T}, 
        #                                       'ub':{t:np.inf for t in T}}})

        if fire_acode:
            for i in ['0', '1']:
                p = fm.inventory(0, mask='? %s ? ?' % i) / fm.inventory(0, mask='? ? ? ?')
                cgen_data.update({'cabrn-thlb%s_%s' % (i, bn):{'lb':{t:param_funcs['cabrn'](bn, fm.base_year+(t-1)*fm.period_length) * p * fm.period_length *
                                                                     (1. - param_funcs['cgen_abrn_e'](bn, fm.base_year+(t-1)*fm.period_length)) for t in T}, 
                                                               'ub':{t:param_funcs['cabrn'](bn, fm.base_year+(t-1)*fm.period_length) * p * fm.period_length for t in T}}})
    fm._tmp['param_funcs'] = param_funcs
    fm._tmp['cgen_data'] = cgen_data
    fm._tmp['cflw_e'] = cflw_e
    fm._tmp['foo'] = coeff_funcs, sense, zexpr, acodes, mask
    return fm.add_problem(name, coeff_funcs, cflw_e, cgen_data=cgen_data, acodes=acodes, sense=sense, mask=mask)


def _gen_scen_maxharv(fm, basenames, name='maxharv', util=0.85, param_funcs=None, target_scalefactors=None, harvest_acode='harvest', 
                      tvy_name='totvol', toffset=0, obj_mode='max_hvol', target_path='./input/targets.csv',
                      max_tp=2020, cbird=None, mask=None, cacut=None):
    from functools import partial
    fm._tmp = {}
    T = fm.periods
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
        param_funcs = {}
        param_funcs['cflw_vcut_e'] = lambda bn, t: 0.1
    coeff_funcs = {'z':partial(cmp_c_z, expr=zexpr)}
    coeff_funcs.update({'cvcut_%s' % bn:partial(cmp_c_caa, expr=vexpr, acodes=[harvest_acode], mask=(bn, '?', '?', '?')) for bn in basenames})
    cflw_e, cgen_data = {}, {}
    fm._tmp['foo'] = param_funcs, toffset, T, basenames
    cflw_ebn = {bn:({t:param_funcs['cflw_vcut_e'](bn, fm.base_year+(t-1)*fm.period_length+toffset) for t in T}, fm.periods[-1]) for bn in basenames}
    cflw_e.update({'cvcut_%s'%bn:cflw_ebn[bn] for bn in basenames})
    fm._tmp['param_funcs'] = param_funcs
    fm._tmp['cgen_data'] = cgen_data
    fm._tmp['cflw_e'] = cflw_e
    fm._tmp['foo'] = coeff_funcs, sense, zexpr, acodes, mask
    return fm.add_problem(name, coeff_funcs, cflw_e, cgen_data=cgen_data, acodes=acodes, sense=sense, mask=mask)


def _gen_scen_maxbird(fm, basenames, name='maxharv', util=0.85, param_funcs=None, target_scalefactors=None, harvest_acode='harvest', 
                      tvy_name='totvol', toffset=0, obj_mode='max_hvol', target_path='./input/targets.csv',
                      max_tp=2020, cacut=None, mask=None, **kwargs):
    from functools import partial
    acodes = ['null', harvest_acode]
    vexpr = '%s * %0.2f' % (tvy_name, util)
    T = fm.periods
    sense = ws3.opt.SENSE_MAXIMIZE 
    if not param_funcs:
        target_scalefactors = {bn:1. for bn in basenames} if not target_scalefactors else target_scalefactors
        df_targets = pd.read_csv(target_path).set_index(['tsa', 'year'])
        param_funcs = {}
        param_funcs['cflw_acut_e'] = lambda bn, t: df_targets.loc[bn, t]['cflw_acut_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cflw_acut_e']
        param_funcs['cflw_vcut_e'] = lambda bn, t: df_targets.loc[bn, t]['cflw_vcut_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cflw_vcut_e']
        #param_funcs['cgen_acut_e'] = lambda bn, t: df_targets.loc[bn, t]['cgen_vcut_e'] if t <= max_tp else df_targets.loc[bn, max_tp]['cgen_vcut_e']
    coeff_funcs = {'z':cmp_c_zbird}
    coeff_funcs.update({'cacut_%s' % bn:partial(cmp_c_caa, expr='1.', acodes=[harvest_acode], mask=(bn, '?', '?', '?')) for bn in basenames})
    coeff_funcs.update({'cvcut_%s' % bn:partial(cmp_c_caa, expr=vexpr, acodes=[harvest_acode], mask=(bn, '?', '?', '?')) for bn in basenames})
    cflw_e, cgen_data = {}, {}
    cflw_ebn = {bn:({t:param_funcs['cflw_acut_e'](bn, fm.base_year+(t-1)*fm.period_length+toffset) for t in T}, fm.periods[-1]) for bn in basenames}
    cflw_e.update({'cacut_%s'%bn:cflw_ebn[bn] for bn in basenames})
    cflw_ebn = {bn:({t:param_funcs['cflw_vcut_e'](bn, fm.base_year+(t-1)*fm.period_length+toffset) for t in T}, fm.periods[-1]) for bn in basenames}
    cflw_e.update({'cvcut_%s'%bn:cflw_ebn[bn] for bn in basenames})
    return fm.add_problem(name, coeff_funcs, cflw_e, cgen_data=cgen_data, acodes=acodes, sense=sense, mask=mask)


def gen_scen(fm, basenames, name, util, param_funcs, target_scalefactors=None, toffset=0, obj_mode='max_hvol',
             cacut=None, mask=None, target_path='./input/targets.csv', cbird=None):
    #print('cbird3', cbird)
    dsp = {'base':_gen_scen_base,
           'maxharv':_gen_scen_maxharv,
           'maxbird':_gen_scen_maxbird,
           'maxharv_fire':_gen_scen_maxharv_fire}
    return dsp[name](fm, basenames, name, util, param_funcs=param_funcs, target_scalefactors=target_scalefactors,
                         toffset=toffset, obj_mode=obj_mode, cacut=cacut, mask=mask, target_path=target_path, cbird=cbird)


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


def custom_round(x, base=1):
    return int(base * round(float(x)/base))

    
def bootstrap_areas(fm, basenames, rst_path, yld_path, hdt, year=None, new_dts=True):
    #print('xxxxx', type(yld_path))
    #au_table = pd.read_csv('%s/au_table.csv' % yld_path).set_index('au_id')
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
            n = 0
            for h, dt in hdt[bn].items():
                ##################################################
                # debug
                #if dt not in [('tsa41', '1', '4102007', '100')]: continue
                #if dt not in [('tsa41', '0', '4103008', '1211')]: continue
                #if dt not in [('tsa41', '0', '4101002', '104')]: continue
                #if dt not in [('tsa41', '1', '4101006', '105')]: continue
                #if dt not in [('tsa41', '1', '4103000', '1201')]: continue
                ##################################################
                #dt = (dt[0], dt[1], dt[2], 's%04d' % int(au_table.loc[int(dt[2])].canfi_species))
                ra = ba[np.where(bh == h)] # match themes hash value
                if new_dts:
                    fm.dtypes[dt] = ws3.forest.DevelopmentType(dt, fm)
                for age in np.unique(ra):
                    #_age = custom_round(age, fm.period_length)
                    _age = age
                    ##################################################
                    # debug
                    #if _age not in [200]: continue
                    ##################################################
                    area = len(ra[np.where(ra == age)]) * pxa
                    _sumarea += area
                    fm.dtypes[dt].area(0, _age, area)
        print('bootstrap_areas', bn, year, pxa, _sumarea)

def jw_string_match(query, choices):
    import jellyfish
    best_match = None
    highest_jw = 0
    for choice in choices:
        current_score = jellyfish.jaro_winkler(query, choice)
        if(current_score > highest_jw):
            highest_jw = current_score
            best_match = choice
    return best_match


def bootstrap_yields(fm, basenames, yld_path, spcode='bcov', 
                    x_max=350, period_length=10., tvy_name='totvol', x_unit='years', province_codes=['bc', 'ab', 'sk', 'mb']):
    def f(x, a, b, c):
        ########################################################
        # This would be the expression to use according to KL readme file:
        #   return a * x ** (b * math.exp(-c * x)) 
        # ... but that yields nonsense output (volumes always near 0).
        # So I replaced with this expression (similar to yield curve form I used in RIA project):
        return a * pow(x,b) * math.exp(-c * x) 
    
    vf = np.vectorize(f)
    #au_table = pd.read_csv('%s/au_table.csv' % yld_path).set_index('au_id')
    gydfs = {pc:pd.read_csv('%s/%s_gy.csv' % (yld_path, pc)).set_index('ygrp') for pc in province_codes}
    key_maps = {pc:{i[:-3]:i for i in gydfs[pc].index} for pc in province_codes}
    for bn in basenames:
        gydf = gydfs[bn[:2]]
        key_map = key_maps[bn[:2]]
        for ygrp in fm.theme_basecodes(2):
            curve_id = ygrp
            mask = ('?', '?', str(curve_id), '?')
            dt_keys = fm.unmask(mask)
            if not dt_keys: continue
            best_match = key_map[jw_string_match(ygrp, key_map.keys())]
            gy = gydf.loc[best_match].iloc[0]
            X = [x for x in range(0, x_max+period_length, period_length)]
            #print(bn, ygrp, period_length)
            #print(X)
            #print(gy.c_a, gy.c_b, gy.c_c)
            #print(gy.d_a, gy.d_b, gy.d_c)
            cY = vf(X, gy.c_a, gy.c_b, gy.c_c)
            dY = vf(X, gy.d_a, gy.d_b, gy.d_c)
            #print(cY)
            #print(dY)
            #assert False
            # add c volume curve
            yname = 'cvol'
            points = zip(X, cY)
            curve = fm.register_curve(ws3.core.Curve(yname, points=points, type='a', 
                                      is_volume=True, xmax=fm.max_age, period_length=period_length))
            fm.yields.append((mask, 'a', [(yname, curve)]))
            fm.ynames.add(yname)
            for dtk in dt_keys: 
                fm.dtypes[dtk].add_ycomp('a', yname, curve)
            # add d volume curve
            yname = 'dvol'
            points = zip(X, dY)
            curve = fm.register_curve(ws3.core.Curve(yname, points=points, type='a', 
                                      is_volume=True, xmax=fm.max_age, period_length=period_length))
            fm.yields.append((mask, 'a', [(yname, curve)]))
            fm.ynames.add(yname)
            for dtk in dt_keys: fm.dtypes[dtk].add_ycomp('a', yname, curve)
    # add total volume curve ###
    mask = ('?', '?', '?', '?')
    expr = '_SUM(%s)' % ', '.join(fm.ynames)
    fm.yields.append((mask, 'c', [(tvy_name, expr)]))
    fm.ynames.add(tvy_name)
    for dtk in fm.dtypes.keys(): fm.dtypes[dtk].add_ycomp('c', tvy_name, expr)  


            
def bootstrap_yields_(fm, yld_path, spcode='canfi_species', 
                        x_max=350, period_length=10., tvy_name='totvol', x_unit='years'):
    #print('yyy', yld_path)
    au_table = pd.read_csv('%s/au_table.csv' % yld_path).set_index('au_id')
    curve_table = pd.read_csv('%s/curve_table.csv' % yld_path)
    curve_points_table = pd.read_csv('%s/curve_points_table.csv' % yld_path).set_index('curve_id')
    #print(au_table.shape)

    # add constants (for bird AF coefficients)    
    c1 = fm.constants['birdaf_swvol_coeff'] =     +2.181e-3
    c2 = fm.constants['birdaf_hwvol_coeff'] =     -1.176e-2
    c3 = fm.constants['birdaf_age_coeff'] =       -8.235e-4
    c4 = fm.constants['birdaf_intercept_coeff'] = +7.594e-1

    swvol_ynames = ['s0105', 's0204', 's0101', 's0304', 's0100', 's0104']
    hwvol_ynames = ['s1201', 's1211']

    for au_id, au_row in au_table.iterrows():
        curve_id = au_row.unmanaged_curve_id # if not is_managed else au_row.managed_curve_id
        mask = ('?', '?', str(curve_id), '?')
        dt_keys = fm.unmask(mask)
        if not dt_keys: continue
        
        # add volume curve
        yname = 's%04d' % int(au_row.canfi_species)
        points = [(r.x, r.y) for _, r in curve_points_table.loc[curve_id].iterrows() if not r.x % period_length and r.x <= x_max]
        curve = fm.register_curve(ws3.core.Curve(yname, points=points, type='a', is_volume=True, xmax=fm.max_age, period_length=period_length))
        fm.yields.append((mask, 'a', [(yname, curve)]))
        fm.ynames.add(yname)
        for dtk in dt_keys: 
            fm.dtypes[dtk].add_ycomp('a', yname, curve)
        
        # add birdaf curve
        curve = curve * (c1 if yname in swvol_ynames else c2) + fm.common_curves['ages'] * c3 + c4
        points = [(x, math.exp(y)) for x, y in curve.points()]
        #####################################################################################################################
        # remove negative y values (else negative values will "cancel out" positive values when rolled up to landscape level)
        # TO DO: confirm that this is "the right thing to do" (i.e., consistent with statistical interpretation of bird AF)
        #points = [(x, max(0., y)) for x, y in curve.points()]
        #####################################################################################################################
        curve.add_points(points=points, compile_y=True) 
        yname = 'birdaf'
        fm.yields.append((mask, 'a', [(yname, curve)]))
        fm.ynames.add(yname)
        for dtk in dt_keys: 
            fm.dtypes[dtk].add_ycomp('a', yname, curve)
                
    mask = ('?', '?', '?', '?')
                
    # add total volume curve ###
    expr = '_SUM(%s)' % ', '.join(fm.ynames)
    fm.yields.append((mask, 'c', [(tvy_name, expr)]))
    fm.ynames.add(tvy_name)
    for dtk in fm.dtypes.keys(): fm.dtypes[dtk].add_ycomp('c', tvy_name, expr)  
    
    ## add softwood volume curve
    #yname = 'swvol'
    #expr = '_SUM(%s)' % ', '.join(swvol_ynames)
    #fm.yields.append((mask, 'c', [(yname, expr)]))
    #fm.ynames.add(yname)
    #for dtk in fm.dtypes.keys(): fm.dtypes[dtk].add_ycomp('c', yname, expr)
    
    # add hardwood volume curve
    #yname = 'hwvol'
    #expr = '_SUM(%s)' % ', '.join(hwvol_ynames)
    #fm.yields.append((mask, 'c', [(yname, expr)]))
    #fm.ynames.add(yname)
    #for dtk in fm.dtypes.keys(): fm.dtypes[dtk].add_ycomp('c', yname, expr)  
    
    
if 0: # old definition (for reference... delete when done)                    
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
    #print('xxx', yld_path)
    bootstrap_areas(fm, basenames, tif_path, yld_path, hdt)
    bootstrap_yields(fm, basenames, yld_path, tvy_name=tvy_name, period_length=yields_period_length, x_unit=yields_x_unit)
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

def rasterize_stands(shp_path, tif_path, theme_cols, age_col, blk_col='', age_divisor=1., d=100.,
                     dtype=rasterio.int32, compress='lzw', round_coords=True,
                     value_func=lambda x: re.sub(r'(-| )+', '_', str(x).lower()), cap_age=None, round_age=1.,
                     verbose=False):
    """
    Rasterize vector stand data.
    """
    import fiona
    from rasterio.features import rasterize
    if verbose: print('rasterizing', shp_path)
    if dtype == rasterio.int32: 
        nbytes = 4
    else:
        raise TypeError('Data type not implemented: %s' % dtype)
    hdt = {}
    shapes = [[], [], []]
    crs = None
    with fiona.open(shp_path, 'r') as src:
        crs = src.crs
        b = src.bounds #(x_min, y_min, x_max, y_max)
        w, h = b[2] - b[0], b[3] - b[1]
        m, n = int((h - (h%d) + d) / d), int((w - (w%d) + d) /  d)
        W = b[0] - (b[0]%d) if round_coords else b[0]
        N = b[1] - (b[1]%d) +d*m if round_coords else b[1] + d*m
        transform = rasterio.transform.from_origin(W, N, d, d)
        for i, f in enumerate(src):
            fp = f['properties']
            fp['thlb'] = 1
            dt = tuple(value_func(fp[t]) for t in theme_cols)
            h = ws3.common.hash_dt(dt, dtype, nbytes)
            hdt[h] = dt
            try:
                age = np.int32(math.ceil(fp[age_col]/float(age_divisor)))
                age = np.int32(custom_round(age, round_age))
            except:
                #######################################
                # DEBUG
                # print(i, fp)                
                #######################################
                if fp[age_col] == None: 
                    age = np.int32(1)
                else:
                    raise ValueError('Bad age value in record %i: %s' % (i, str(fp[age_col])))
            if cap_age and age > cap_age: age = cap_age
            #try:
            #    assert age > 0
            #except:
            #    if fp[age_col] == 0:
            #        age = np.int32(1)
            #    else:
            #        print('bad age', age, fp[age_col], age_divisor)
            #        raise
            blk = i if not blk_col else fp[blk_col]
            shapes[0].append((f['geometry'], h))   # themes
            shapes[1].append((f['geometry'], age)) # age
            shapes[2].append((f['geometry'], blk)) # block identifier
    #rst_path = shp_path[:-4]+'.tif' if not rst_path else rst_path
    nodata_value = -2147483648
    kwargs = {'out_shape':(m, n), 'transform':transform, 'dtype':dtype, 'fill':nodata_value}
    r = np.stack([rasterize(s, **kwargs) for s in shapes])
    kwargs = {'driver':'GTiff', 
              'width':n, 
              'height':m, 
              'count':3, 
              'crs':crs,
              'transform':transform,
              'dtype':dtype,
              'nodata':nodata_value,
              'compress':compress}
    #print(shp_path)
    #print(src.crs)
    #print(kwargs)
    with rasterio.open(tif_path, 'w', **kwargs) as snk:
        snk.write(r[0], indexes=1)
        snk.write(r[1], indexes=2)
        snk.write(r[2], indexes=3)
    return hdt
            
def rasterize_inventory(basenames, shp_path, tif_path, hdt_path, theme_cols, age_col, base_year,
                        age_divisor=1., round_age=1., cap_age=None, d=90., verbose=True):
    hdt = {}
    for bn in basenames:
        kwargs = {'shp_path':'%s' % shp_path(bn), 
                  'tif_path':'%s/inventory_init.tif' % tif_path(bn), 
                  'theme_cols':theme_cols, 
                  'age_col':age_col, 
                  'age_divisor':age_divisor,
                  'round_age':round_age,
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


def schedule_harvest_optimize(fm, basenames, scenario_name='maxharv', util=0.85, target_scalefactors=None, param_funcs=None, 
                              target_path='./input/targets.csv', obj_mode='max_hvol', mask=None, cbird=None):
    import gurobipy as grb
    #print('cbird2', cbird)
    p = gen_scen(fm, basenames, scenario_name, util, target_scalefactors=target_scalefactors, param_funcs=param_funcs, toffset=0, 
                 obj_mode=obj_mode, mask=mask, target_path=target_path, cbird=cbird)
    m = p.solve()
    if m.status != grb.GRB.OPTIMAL:
        print('Model not optimal.')
        return None

    sch = fm.compile_schedule(p)
    fm.reset()
    fm.apply_schedule(sch, 
                      force_integral_area=False, 
                      override_operability=False,
                      fuzzy_age=False,
                      recourse_enabled=False,
                      verbose=False,
                      compile_c_ycomps=True)
    return sch


def schedule_harvest_areacontrol(fm, period=1, acode='harvest', util=0.85, 
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


def sda(fm, basenames, time_step, tif_path, hdt, acode_map=None, nthresh=0, sda_mode='randpxl'):
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
        fr.allocate_schedule(mask=mask, verbose=0, sda_mode=sda_mode, nthresh=nthresh, max_period=None, shuffle_dy=True)
        fr.cleanup()


def pickle_forestmodel(fm, scenario_name, basename):
    pickle.dump(fm, open('dat/out/%s_%s_fm.pkl' % (scenario_name, basename), 'wb'))

    
def pickle_schedule(sch, scenario_name, basename):
    pickle.dump(sch, open('dat/out/%s_%s_sch.pkl' % (scenario_name, basename), 'wb'))


def unpickle_forestmodel(scenario_name, basename):
    return pickle.load(open('dat/out/%s/%s_%s_fm.pkl' % (scenario_name, scenario_name, basename), 'rb'))


def unpickle_schedule(scenario_name, basename):
    return pickle.load(open('dat/out/%s/%s_%s_sch.pkl' % (scenario_name, scenario_name, basename), 'rb'))
