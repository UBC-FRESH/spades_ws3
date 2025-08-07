import rasterio
import numpy as np
import pandas as pd
import shutil
try:
    import cPickle as pickle
except:
    import pickle
import ws3

from concurrent.futures import ProcessPoolExecutor, as_completed


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

def _gen_scen_base(fm, basenames, name,
                   util=0.85, obj_mode='max_hv', tvy_name='totvol', harvest_acode='harvest',
                   cgen_hv=None, cgen_ha=None, cgen_hv_e=None, cgen_ha_e=None, cgen_e_default=0.01,
                   cflw_hv=True, cflw_ha=True, cflw_hv_e=None, cflw_ha_e=None, cflw_e_default=0.05,
                   mask=None, workers=1):
    from functools import partial
    acodes = ['null', harvest_acode]
    vexpr = '%s * %0.2f' % (tvy_name, util)
    T = fm.periods
    cflw_e, cgen_data = {}, {}
    if obj_mode == 'max_hv':
        sense = ws3.opt.SENSE_MAXIMIZE 
        zexpr = vexpr
    elif obj_mode == 'min_ha':
        sense = ws3.opt.SENSE_MINIMIZE 
        zexpr = '1.'
    else:
        raise ValueError('Invalid obj_mode: %s' % obj_mode)
    if not cgen_hv_e:
        cgen_hv_e = {(bn, t):cgen_e_default for bn in basenames for t in T}
    if not cgen_ha_e:
        cgen_ha_e = {(bn, t):cgen_e_default for bn in basenames for t in T}
    if not cflw_hv_e:
        cflw_hv_e = {(bn, t):cflw_e_default for bn in basenames for t in T}
    if not cflw_ha_e:
        cflw_ha_e = {(bn, t):cflw_e_default for bn in basenames for t in T}
    coeff_funcs = {'z':partial(cmp_c_z, expr=zexpr)}
    coeff_funcs.update({'cgen_hv_%s' % bn:partial(cmp_c_caa, expr=vexpr, acodes=[harvest_acode], mask=(bn, '?', '?', '?')) 
                        for bn in basenames})
    coeff_funcs.update({'cgen_ha_%s' % bn:partial(cmp_c_caa, expr='1.', acodes=[harvest_acode], mask=(bn, '?', '?', '?')) 
                        for bn in basenames})
    coeff_funcs.update({'cflw_hv_%s' % bn:partial(cmp_c_caa, expr=vexpr, acodes=[harvest_acode], mask=(bn, '?', '?', '?')) 
                        for bn in basenames})
    coeff_funcs.update({'cflw_ha_%s' % bn:partial(cmp_c_caa, expr='1.', acodes=[harvest_acode], mask=(bn, '?', '?', '?')) 
                        for bn in basenames})
    for bn in basenames:
        if cgen_hv:
            cgen_data.update({'cgen_hv_%s' % bn:{'lb':{t:cgen_hv[bn, t] * (1. - cgen_hv_e[bn, t])  for t in T},
                                                 'ub':{t:cgen_hv[bn, t] for t in T}}})
        else: # flow constraints only guaranteed to be feasible if no arbitrary cgen constraints on harvest volume
            cflw_e.update({'cflw_hv_%s' % bn:({t:cflw_hv_e[bn, t] for t in T}, 1)})

        if cgen_ha:
            cgen_data.update({'cgen_ha_%s' % bn:{'lb':{t:cgen_ha[bn, t] * (1. - cgen_ha_e[bn, t])  for t in T},
                                                 'ub':{t:cgen_ha[bn, t] for t in T}}})

        cflw_e.update({'cflw_ha_%s' % bn:({t:cflw_ha_e[bn, t] for t in T}, 1)})
    p = fm.add_problem(name, 
                       coeff_funcs, 
                       cflw_e=cflw_e, 
                       cgen_data=cgen_data, 
                       acodes=acodes, 
                       sense=sense,
                       mask=mask,
                       workers=workers)
    return p

def _build_unit_problem_task(args):
    """
    Worker function for ProcessPoolExecutor to build a Problem for a single management unit.

    Parameters
    ----------
    args : tuple
        (fm, basenames, name, util, obj_mode, tvy_name,
         cgen_hv, cgen_ha, cgen_hv_e, cgen_ha_e, cgen_e_default,
         cflw_hv, cflw_ha, cflw_hv_e, cflw_ha_e, cflw_e_default,
         unit_theme_idx, unit, workers)

    Returns
    -------
    tuple
        (unit_name, ws3.opt.Problem)
    """
    (fm, basenames, name, util, obj_mode, tvy_name,
     cgen_hv, cgen_ha, cgen_hv_e, cgen_ha_e, cgen_e_default,
     cflw_hv, cflw_ha, cflw_hv_e, cflw_ha_e, cflw_e_default,
     unit_theme_idx, unit, workers) = args

    dsp = {'base': _gen_scen_base}

    # Compute mask for this unit
    local_mask = tuple(
        unit if idx == unit_theme_idx else "?"
        for idx in range(len(fm._themes))
    )

    # Build the Problem for this unit
    problem = dsp['base'](
        fm=fm,
        basenames=[unit],  # Only this unit's base
        name=f"{name}_{unit}",
        util=util,
        obj_mode=obj_mode,
        tvy_name=tvy_name,
        cgen_hv=cgen_hv,
        cgen_ha=cgen_ha,
        cgen_hv_e=cgen_hv_e,
        cgen_ha_e=cgen_ha_e,
        cgen_e_default=cgen_e_default,
        cflw_hv=cflw_hv,
        cflw_ha=cflw_ha,
        cflw_hv_e=cflw_hv_e,
        cflw_ha_e=cflw_ha_e,
        cflw_e_default=cflw_e_default,
        mask=local_mask,
        workers=workers
    )

    return (unit, problem)

def gen_scen(
    fm,
    basenames,
    name,
    util=0.85,
    obj_mode='max_hv',
    tvy_name='totvol',
    cgen_hv=None,
    cgen_ha=None,
    cgen_hv_e=None,
    cgen_ha_e=None,
    cgen_e_default=0.01,
    cflw_hv=True,
    cflw_ha=True,
    cflw_hv_e=None,
    cflw_ha_e=None,
    cflw_e_default=0.05,
    mask=None,
    mgmt_unit_theme=None,
    workers=1
):
    """
    Dispatcher for generating ws3 optimization problems.

    Parameters
    ----------
    fm : ForestModel
        Forest model instance to generate problems from.
    basenames : list
        List of base development type names to include (typically per-unit).
    name : str
        Scenario name.
    util : float
        Utilization factor for volume outputs.
    obj_mode : str
        'max_hv' or 'min_ha' objective.
    tvy_name : str
        Output name used for harvest volume.
    cgen_hv : dict or None
        General constraints (harvest volume) as {(basename, t): value}.
    cgen_ha : dict or None
        General constraints (harvest area).
    cgen_hv_e, cgen_ha_e : dict or None
        Even-flow epsilon definitions for volume/area.
    cgen_e_default : float
        Default epsilon for general constraints.
    cflw_hv, cflw_ha : bool
        Whether to include flow constraints on harvest volume or area.
    cflw_hv_e, cflw_ha_e : dict or None
        Epsilon dictionaries for flow constraints.
    cflw_e_default : float
        Default epsilon for flow constraints.
    mask : tuple or None
        Optional mask for filtering development types.
    mgmt_unit_theme : int or None
        Theme index for per-management-unit decomposition.
    workers : int
        Number of CPU cores to use for problem generation.

    Returns
    -------
    ws3.opt.Problem or dict
        Single Problem (default) or dict of {unit_name: Problem} if mgmt_unit_theme is specified.
    """
    dsp = {'base': _gen_scen_base}

    # --- Single-unit behavior ---
    if mgmt_unit_theme is None:
        return dsp[name](
            fm=fm,
            basenames=basenames,
            name=name,
            util=util,
            obj_mode=obj_mode,
            tvy_name=tvy_name,
            cgen_hv=cgen_hv,
            cgen_ha=cgen_ha,
            cgen_hv_e=cgen_hv_e,
            cgen_ha_e=cgen_ha_e,
            cgen_e_default=cgen_e_default,
            cflw_hv=cflw_hv,
            cflw_ha=cflw_ha,
            cflw_hv_e=cflw_hv_e,
            cflw_ha_e=cflw_ha_e,
            cflw_e_default=cflw_e_default,
            mask=mask,
            workers=workers
        )

    # --- Multi-unit decomposition ---
    unit_theme_idx = int(mgmt_unit_theme)
    unit_codes = fm._theme_basecodes[unit_theme_idx]
    print(f"gen_scen: multi-unit mode on theme {unit_theme_idx}, units={unit_codes}")

    max_outer_workers = min(len(unit_codes), workers if workers > 1 else 1)
    workers_per_unit = max(1, workers // max_outer_workers)
    problems = {}

    args_list = []
    for unit in unit_codes:
        # --- Filter cgen_hv for this unit ---
        if cgen_hv is not None:
            unit_cgen_hv = {
                (bn, t): val
                for (bn, t), val in cgen_hv.items()
                if bn == unit
            }
        else:
            unit_cgen_hv = None

        # Build args tuple
        args_list.append((
            fm, basenames, name, util, obj_mode, tvy_name,
            unit_cgen_hv, cgen_ha, cgen_hv_e, cgen_ha_e, cgen_e_default,
            cflw_hv, cflw_ha, cflw_hv_e, cflw_ha_e, cflw_e_default,
            unit_theme_idx, unit, workers_per_unit
        ))

    # Launch tasks in parallel if requested
    if max_outer_workers > 1:
        with ProcessPoolExecutor(max_workers=max_outer_workers) as executor:
            futures = {executor.submit(_build_unit_problem_task, args): args[17] for args in args_list}
            for fut in as_completed(futures):
                unit, problem = fut.result()
                problems[unit] = problem
    else:
        for args in args_list:
            unit, problem = _build_unit_problem_task(args)
            problems[unit] = problem

    return problems
    
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
    
def bootstrap_areas(fm, basenames, rst_path, hdt, year=None, new_dts=True):
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

def bootstrap_yields(fm, yld_path, spcode='canfi_species', 
                     x_max=350, period_length=10., tvy_name='totvol', x_unit='years'):
    au_table = pd.read_csv('%s/au_table.csv' % yld_path).set_index('au_id')
    curve_table = pd.read_csv('%s/curve_table.csv' % yld_path)
    curve_points_table = pd.read_csv('%s/curve_points_table.csv' % yld_path).set_index('curve_id')
    print(au_table.shape)
    for au_id, au_row in au_table.iterrows():
        yname = 's%04d' % int(au_row.canfi_species)
        for is_managed in [0]:
            curve_id = au_row.unmanaged_curve_id if not is_managed else au_row.managed_curve_id
            mask = ('?', '?', str(curve_id), '?')
            dt_keys = fm.unmask(mask)
            if not dt_keys: continue
            points = [(r.x, r.y) for _, r in curve_points_table.loc[curve_id].iterrows() if not r.x % period_length and r.x <= x_max]
            c = fm.register_curve(ws3.core.Curve(yname, points=points, type='a', is_volume=True, xmax=fm.max_age, period_length=period_length))
            fm.yields.append((mask, 'a', [(yname, c)]))
            fm.ynames.add(yname)
            for dtk in dt_keys: 
                fm.dtypes[dtk].add_ycomp('a', yname, c)
    # add total volume curve ###
    expr = '_SUM(%s)' % ', '.join(fm.ynames)
    fm.yields.append((('?', '?', '?', '?'), 'c', [(tvy_name, expr)]))
    fm.ynames.add(tvy_name)
    for dtk in fm.dtypes.keys(): fm.dtypes[dtk].add_ycomp('c', tvy_name, expr)
                    
def bootstrap_yields_(fm, yld_path, theme_cols=['AU', 'LDSPP'], spcode='SPCode', 
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

def _solve_stage(fm, problems, stage_label, workers=1, warm_starts=None):
    """
    Solve one optimization stage for single or multi-unit problems in parallel if requested.

    Parameters
    ----------
    fm : ForestModel
        The forest model instance (used for schedule compilation if serial).
    problems : ws3.opt.Problem or dict
        Single Problem object or dict of {unit_name: Problem}.
    stage_label : str
        Label printed in logs to identify the stage.
    workers : int, optional
        Number of cores to use for parallel solving across units. Default = 1.
    warm_starts : dict[str, list[float]] or list[float] or None
        Optional warm start solutions:
            - For multi-unit: {unit_name: col_value_list}
            - For single-unit: list of col values
            - None to disable warm starts

    Returns
    -------
    list
        Flattened list of schedule tuples sorted by period.
    """
    schedules = []

    # --- Multi-unit mode ---
    if isinstance(problems, dict):
        units = list(problems.keys())
        n_units = len(units)
        print(f"{stage_label}: solving {n_units} unit problems")

        # Parallel if >1 unit and >1 worker
        outer_workers = min(workers, n_units)
        if outer_workers > 1:
            # Prepare (unit, problem, warm_start_vector) tuples
            args_list = []
            for unit in units:
                ws_vec = warm_starts[unit] if warm_starts and unit in warm_starts else None
                args_list.append((fm, unit, problems[unit], ws_vec))

            with ProcessPoolExecutor(max_workers=outer_workers) as executor:
                futures = {
                    executor.submit(_solve_unit_task_with_warm_start, args): args[1] for args in args_list
                }
                for fut in as_completed(futures):
                    unit, schedule = fut.result()
                    assert schedule
                    schedules.extend(schedule)
        else:
            # Serial solve
            for unit, prob in problems.items():
                ws_vec = warm_starts[unit] if warm_starts and unit in warm_starts else None
                print(f"{stage_label} solving for unit {unit}")
                prob.solve(threads=0, warm_start=ws_vec)
                schedule = fm.compile_schedule(prob)
                assert schedule
                schedules.extend(schedule)

    else:
        # --- Single problem mode ---
        print(f"{stage_label}: solving single problem")
        ws_vec = warm_starts if warm_starts is not None and warm_starts.any() else None
        problems.solve(threads=0, warm_start=ws_vec)
        schedules = fm.compile_schedule(problems)
        assert schedules

    # Sort merged schedules by period (index 4 in tuple)
    schedules.sort(key=lambda x: x[4])
    return schedules

def _solve_unit_task_with_warm_start(args):
    """
    Top-level helper for ProcessPoolExecutor to solve a single Problem with optional warm start.
    Returns (unit_name, schedule_list).
    """
    fm, unit, problem, warm_start = args
    print(f"Solving problem for unit {unit} (warm start: {warm_start is not None})")
    problem.solve(threads=0, warm_start=warm_start)
    schedule = fm.compile_schedule(problem)
    return (unit, schedule)

def schedule_harvest_optimize(
    fm,
    basenames,
    scenario_name='base',
    tvy_name='totvol',
    util=0.85,
    p_max_hv={},
    mask=None,
    mgmt_unit_theme=None,
    workers=1
):
    """
    Run a two-stage harvest scheduling optimization:

    Stage 1
        Maximize even-flow harvest volumes across the planning horizon.

    Stage 2
        Minimize harvest area subject to the harvest volume constraints
        obtained from Stage 1.

    Parameters
    ----------
    fm : ForestModel
        The forest model instance to optimize.
    basenames : list
        List of base development type names to include.
    scenario_name : str
        Name of the scenario (used for problem names and reporting).
    tvy_name : str
        Name of the output used to measure harvest volume (objective).
    util : float
        Utilization factor (0.0 - 1.0) applied to volume outputs.
    p_max_hv : dict
        Optional dict mapping basename → scaling factor for harvest volume limits.
    mask : tuple or None
        Optional mask for filtering development types.
    mgmt_unit_theme : int or None
        Optional theme index for per-management-unit decomposition.
    workers : int
        Number of CPU cores to use for parallel problem generation.

    Returns
    -------
    list
        Compiled Stage 2 schedule as a list of tuples sorted by period.
    """
    ########################################
    # Stage 1: Maximize harvest volume (even-flow)
    ########################################
    print('schedule_harvest_optimize: stage 1, generating problem')
    p1 = gen_scen(
        fm=fm,
        basenames=basenames,
        name=scenario_name,
        util=util,
        mgmt_unit_theme=mgmt_unit_theme,
        workers=workers
    )

    schedules_stage1 = _solve_stage(
        fm, p1, "schedule_harvest_optimize: stage 1", workers
    )

    # --- Build warm start vectors for Stage 2 ---
    if isinstance(p1, dict):
        warm_starts = {
            unit: np.array([var.val or 0.0 for var in problem._vars.values()], dtype=np.float64)
            for unit, problem in p1.items()
        }
    else:
        warm_starts = np.array([var.val or 0.0 for var in p1._vars.values()], dtype=np.float64)

    # --- Reset and apply Stage 1 schedule to forest model ---
    fm.reset()
    fm.apply_schedule(
        schedules_stage1,
        force_integral_area=True,
        override_operability=True,
        fuzzy_age=True,
        recourse_enabled=True,
        verbose=False,
        compile_c_ycomps=True
    )

    ########################################
    # Build harvest volume constraints for Stage 2
    ########################################
    vexpr = f"{tvy_name} * {util:0.2f}"
    hv_coeffs = {bn: p_max_hv.get(bn, 1.0) for bn in basenames}
    cgen_hv = {
        (bn, t): fm.compile_product(
            t, vexpr,
            dtype_keys=fm.unmask((bn, '?', '?', '?'))
        ) * hv_coeffs[bn]
        for bn in basenames
        for t in fm.periods
    }

    ########################################
    # Stage 2: Minimize harvest area
    ########################################
    print('schedule_harvest_optimize: stage 2, generating problem')
    p2 = gen_scen(
        fm=fm,
        basenames=basenames,
        name=scenario_name,
        util=util,
        obj_mode='min_ha',
        cgen_hv=cgen_hv,
        mgmt_unit_theme=mgmt_unit_theme,
        workers=workers
    )

    schedules_stage2 = _solve_stage(
        fm, p2, "schedule_harvest_optimize: stage 2", workers,
        warm_starts=warm_starts
    )

    # --- Reset and apply Stage 2 schedule ---
    fm.reset()
    fm.apply_schedule(
        schedules_stage2,
        force_integral_area=True,
        override_operability=True,
        fuzzy_age=True,
        recourse_enabled=True,
        verbose=False,
        compile_c_ycomps=True
    )

    return schedules_stage2

def schedule_harvest_areacontrol(fm, period=None, acode='harvest', util=0.85, 
                                 target_scalefactors=None,
                                 mask_area_thresh=0.,
                                 verbose=0):
    fm.reset() #fm.reset_actions()
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
    target_areas = []
    for i, mask in enumerate(target_masks): # compute area-weighted mean CMAI age for each masked DT set
        masked_area = fm.inventory(0, mask=mask, verbose=verbose)
        if not masked_area: continue
        r = sum((fm.dtypes[dtk].ycomp('totvol').mai().ytp().lookup(0) * fm.dtypes[dtk].area(0)) for dtk in fm.unmask(mask))
        r /= masked_area
        _target_scalefactor = 1.
        if target_scalefactors and isinstance(target_scalefactors, dict):
            for _mask in [(bn, '1', '?', '?') for bn in target_scalefactors]:
                try:
                    if fm.match_mask(_mask, fm.unmask(mask)[0]): _target_scalefactor = target_scalefactors[_mask[0]]
                except:
                    pass
        asf = _target_scalefactor  
        ta = (1/r) * masked_area * asf * fm.period_length
        target_areas.append(ta)
    periods = fm.periods if not period else [period]
    for period in periods:
        for mask, target_area in zip(target_masks, target_areas):
            if verbose > 0:
                print('calling areaselector', period, acode, target_area, mask)
            fm.areaselector.operate(period, acode, target_area, mask=mask, verbose=verbose)
    sch = fm.compile_schedule()
    return sch

def sda(fm, basenames, time_step, tif_path, hdt, acode_map=None, nthresh=10, 
        sda_mode='randblk', horizon=1, verbose=False):
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
                     'horizon':horizon,
                     'period_length':fm.period_length,
                     'time_step':time_step,
                     'base_year':fm.base_year,
                     'piggyback_acodes':{}}
        return fr_kwargs
    for bn in basenames:
        print('SDA for TSA', bn)
        mask = (bn, '?', '?', '?')
        fr = ForestRaster(**cmp_fr_kwargs(bn))
        fr.allocate_schedule(mask=mask, verbose=verbose, sda_mode=sda_mode, nthresh=nthresh)
        fr.cleanup()

def pickle_forestmodel(fm, scenario_name, basename):
    pickle.dump(fm, open('dat/out/%s_%s_fm.pkl' % (scenario_name, basename), 'wb'))

    
def pickle_schedule(sch, scenario_name, basename):
    pickle.dump(sch, open('dat/out/%s_%s_sch.pkl' % (scenario_name, basename), 'wb'))

def unpickle_forestmodel(scenario_name, basename):
    return pickle.load(open('dat/out/%s/%s_%s_fm.pkl' % (scenario_name, scenario_name, basename), 'rb'))

def unpickle_schedule(scenario_name, basename):
    return pickle.load(open('dat/out/%s/%s_%s_sch.pkl' % (scenario_name, scenario_name, basename), 'rb'))
    return pickle.load(open('dat/out/%s/%s_%s_sch.pkl' % (scenario_name, scenario_name, basename), 'rb'))
