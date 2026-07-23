# %%
import numpy as np
import os, sys, tempfile

sys.path.append('../../physion/src') # add src code directory for physion
sys.path.append('../../utils') 

import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')

import plots_fcts as pt_fcts
import params
import tools

from scipy import stats
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import sem
from itertools import product
import datetime
from sklearn.metrics import auc

# %%

# TO LOOP OVER NWB FILES WITH VISUAL STIMULUS --- DRIFITING GRATING ---  multisession
pname = "drifting-grating"
varied_parameter = 'contrast'
contrast_values = [0.2, 0.6, 1.]

folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna',
                        'PN_cond-NDNF-CB1_WT-vs-KD', 'final_NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol=pname)

#%% SET PARAMETERS

state_metric = 'speed' # 'speed' or 'pupil' or 'speed & pupil'

stat_test_props = params.stat_test_props.copy()
stat_test_props['interval_post'] = params.interval_post[pname]
stat_test_props['sign'] = 'positive'

dFoF_options = params.dFoF_options
dFoF_options['with_computed_neuropil_fact'] = True

ep_props = dict(quantities=['dFoF', 'Deconvolved', 'running_speed'],
                #prestim_duration=1.5,
                dt_sampling=params.dt_sampling,
                verbose=False)

savepath_fig = os.path.join(r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\PN\figures", pname)
savepath_excel = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\PN\excels"
savepath_data = os.path.join(r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\PN\data", pname)

os.makedirs(savepath_fig, exist_ok=True)
os.makedirs(savepath_data, exist_ok=True)
os.makedirs(savepath_excel, exist_ok=True)

# %% COMPUTATION

# Initialization
DATASET["virus"], DATASET['alpha'], DATASET["nb_neurons"] = [], [], []
for c in contrast_values:
    DATASET["nb_responsive_neurons_%s" % c] = []
    DATASET['mean_dFoF_%s' % c] = []

included_mice = None
dFoF = None
deconvolved = None
run_means = None
viruses = ['sgRosa', 'sgCnr1']

# LOOP OVER SESSIONS

for i, filename in enumerate(DATASET['files']):
#for i, filename in enumerate([DATASET['files'][0]]):
    
    data = physion.analysis.read_NWB.Data(filename, verbose=False)

    # determine virus        
    if 'sgRosa' in data.nwbfile.virus:
        virus = 'sgRosa'
    elif 'sgCnr1' in data.nwbfile.virus:
        virus = 'sgCnr1'
    else :
        raise ValueError("Virus not identified in session %s" % filename)
    DATASET["virus"].append(virus)
    
    print(i+1, '--', filename, '--', data.nROIs)

    if 'dFoF' in ep_props['quantities']:
        data.build_dFoF(**dFoF_options, verbose=True)
    DATASET["alpha"].append(data.neuropil_correction_factor)
    
    DATASET["nb_neurons"].append(data.nROIs)

    if data.nROIs>0:

        if 'Deconvolved' in ep_props['quantities']:
            data.build_Deconvolved()
        if 'pupil_diameter' in ep_props['quantities']:
            data.build_pupil_diameter()
        if 'facemotion' in ep_props['quantities']:
            data.build_facemotion()
        if 'running_speed' in ep_props['quantities']: 
            data.build_running_speed()

        # Build episode data
        ep = physion.analysis.episodes.build.EpisodeData(data, 
                                                         **ep_props,
                                                         protocol_name=pname)

        if varied_parameter in ep.varied_parameters:

            ######### 1) Define behavioral states #########
            states_names, states_filters = tools.define_trials_arousal_state(ep, cond=state_metric)

            ######### 2) identify visually-responsive cells #########

            evokedStats = ep.pre_post_statistics(\
                                                stat_test_props,
                                                response_args=params.response_args,
                                                response_significance_threshold=params.response_significance_threshold,
                                                loop_over_cells=True,
                                                repetition_keys=['repeat'],
                                                verbose=False
                                                )

            # find responsive ROIs for this contrast (from summary stats)
            for k, vparam in enumerate(contrast_values):
                DATASET["nb_responsive_neurons_%s" % vparam].append(np.sum(evokedStats['significant'][:, k], axis=0) if evokedStats['significant'][:, k].size != 0 else 0)

                window = (ep.t>stat_test_props['interval_post'][0]) & (ep.t<stat_test_props['interval_post'][1])
                if evokedStats['significant'][:, k].size != 0:
                    DATASET['mean_dFoF_%s' % vparam].append(np.mean(ep.dFoF[:, evokedStats['significant'][:, k], :][:, :, window]))
                else :
                    DATASET['mean_dFoF_%s' % vparam].append(np.nan)

            ######### 3) Get episodes' traces #########

            # INITILIAZE DICTIONARIES TO STORE RESPONSES AND BEHAVIORAL QUANTITIES
            if included_mice is None:

                vparam_values = contrast_values
                included_mice  = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}
                dFoF = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}
                deconvolved = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}
                run_means = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}

            for k, (vparam, vparam_exact) in enumerate(zip(contrast_values, ep.varied_parameters[varied_parameter])):
                                        
                vparam_cond = (getattr(ep, varied_parameter)==vparam_exact)

                for state, state_filter in zip(states_names, states_filters):

                    key = f"{virus}-{state}-{vparam}"

                    if (evokedStats['significant'][:, k].size != 0) and \
                        (np.sum(evokedStats['significant'][:, k], axis=0)>=params.NMIN_ROIS) and \
                            (np.sum(vparam_cond & state_filter) >= params.NMIN_EPISODES):
                        
                        print("cond: %s-%s -> included %i ROIs and %i episodes" % 
                            (virus, state, np.sum(evokedStats['significant'][:, k]), np.sum(vparam_cond & state_filter)))
                        print("    -> %i ROIs out of %i ROIs are responsive" % 
                            (np.sum(evokedStats['significant'][:, k], axis=0), evokedStats['significant'][:, k].shape[0]))
                        
                        dFoF[key].append(
                            ep.dFoF[vparam_cond & state_filter][:, evokedStats['significant'][:, k], :])
                        
                        deconvolved[key].append(
                            ep.Deconvolved[vparam_cond & state_filter][:, evokedStats['significant'][:, k], :])
                        
                        included_mice[key].append(DATASET['subjects'][i])

                        run_means[key].append(ep.running_speed[vparam_cond & state_filter, :])
                    
                    else:
                        print("cond: %s-%s -> [XX] response not included (%i ROIs, %i eps)" % 
                            (virus,state, np.sum(evokedStats['significant'][:, k]), np.sum(vparam_cond & state_filter)))

            ###############################

        else :
            for c in contrast_values:
                DATASET["nb_responsive_neurons_%s" % c].append(np.nan)
                DATASET['mean_dFoF_%s' % c].append(np.nan)
    
    else :
        print("session %s has no ROIs, excluded from analysis" % filename)
        for c in contrast_values:
            DATASET["nb_responsive_neurons_%s" % c].append(0)
            DATASET['mean_dFoF_%s' % c].append(np.nan)
    
    print('')
                
dFoF = tools.remove_empty_sessions(dFoF)
deconvolved = tools.remove_empty_sessions(deconvolved)
run_means = tools.remove_empty_sessions(run_means)

# SAVE IN NPY FILES
np.save(os.path.join(savepath_data, pname + '_dFoF_means.npy'), dFoF)
np.save(os.path.join(savepath_data, pname + '_deconvolved_means.npy'), deconvolved)
np.save(os.path.join(savepath_data, pname + '_run_means.npy'), run_means)
np.save(os.path.join(savepath_data, pname + '_included_mice.npy'), included_mice)

# BUILD DATAFRAME
import pandas as pd 

DATASET['session'] = np.array([os.path.basename(file)[:-4] for file in DATASET['files']])

df = pd.DataFrame(DATASET).drop(columns=['protocol_ids', 'protocols', 'files']).set_index('session')
df.dropna(axis=0, thresh=len(contrast_values)*4, inplace=True)
for c in contrast_values:
    df['perc_resp_%s' % c] = df['nb_responsive_neurons_%s' % c] / df['nb_neurons'] * 100

excel_filename = pname + '_summary_data_' + 'PYR' + '.xlsx'
df.to_excel(os.path.join(savepath_excel, excel_filename))

df

#%% 1.a Averaged dF/F0 over responsive ROIs and over episodes across virus and behavioral states (std over sessions)

fig, AX = pt_fcts.plot_average_response(ep.t, dFoF, 
                                        viruses, states_names, vparam_values, varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS)

firgurename = pname + '_average_res_dfof_'+ 'PYR' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")

#%% 1.b Averaged dF/F0 with baseline substracted over responsive ROIs and over episodes across virus and behavioral states (std over sessions)
baselineCond = (ep.t>stat_test_props['interval_pre'][0]) & (ep.t<stat_test_props['interval_pre'][1])

fig, AX = pt_fcts.plot_average_response(ep.t, dFoF, 
                                        viruses, states_names, vparam_values, varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS, 
                                        baselineSubtraction=True, baselineCond=baselineCond)

firgurename = pname + '_average_res_dfof_bsl_sub_'+ 'PYR' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")

#%% 1.c Averaged deconvolved trace over responsive ROIs and over episodes across virus and behavioral states (std over sessions)

fig, AX = pt_fcts.plot_average_response(ep.t, deconvolved, 
                                        viruses, states_names, vparam_values, varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS)
for i in range(len(contrast_values)):
    AX[i][0].set_ylabel('deconvolved')

firgurename = pname + '_average_res_deconvolved_' + 'PYR' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")

#%% 2. Pie chart of responsive neurons
percentages = {}

for v in viruses:
    percentages[v] = []
    for c in contrast_values:
        percentages[v].append((df[df['virus'] == v]['perc_resp_%s' %c ]).to_list())

    percentages[v] = np.array(percentages[v]).T

fig, AX = pt_fcts.pie_chart_responsive_neurons(percentages, viruses, vparam_values)

firgurename = pname + '_pie_charts_' + 'PYR' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")