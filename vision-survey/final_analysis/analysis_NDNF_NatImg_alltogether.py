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
pname = 'Natural-Images-4-repeats'
varied_parameter = 'Image-ID'

folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna',
                        'NDNF-cond-CB1', 'NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol=pname)

#%% SET PARAMETERS

pos_stat_test_props = params.stat_test_props.copy()
pos_stat_test_props['interval_post'] = params.interval_post[pname]
pos_stat_test_props['sign'] = 'positive'

neg_stat_test_props = params.stat_test_props.copy()
neg_stat_test_props['interval_post'] = params.interval_post[pname]
neg_stat_test_props['sign'] = 'negative'

dFoF_options = params.dFoF_options
dFoF_options['with_computed_neuropil_fact'] = False

ep_props = dict(quantities=['dFoF', 'Deconvolved'] ,
                #prestim_duration=1.5,
                dt_sampling=params.dt_sampling,
                verbose=False)

savepath_fig = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\NDNF\figures\natural_images"
savepath_data = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\NDNF\excels"

# %% COMPUTATION

# Initialization
DATASET["virus"], DATASET['alpha'], DATASET["nb_neurons"] = [], [], []
DATASET["nb_responsive_neurons_pos"], DATASET["nb_responsive_neurons_neg"] = [], []
DATASET['mean_dFoF_pos'], DATASET['mean_dFoF_neg'] = [], []

included_mice_pos = None
included_mice_neg = None
dFoF_pos = None
dFoF_neg = None
deconvolved_pos = None
deconvolved_neg = None
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

        ######### 1) Define behavioral states #########
        states_names, states_filters = ['all'], [np.ones(ep.repeat.shape[0], dtype=bool)]

        ######### 2) identify visually-responsive cells #########

        pos_evokedStats = ep.pre_post_statistics(\
                                            pos_stat_test_props,
                                            response_args=params.response_args,
                                            response_significance_threshold=params.response_significance_threshold,
                                            loop_over_cells=True,
                                            repetition_keys=[varied_parameter, 'repeat'],
                                            verbose=False
                                            )
        
        neg_evokedStats = ep.pre_post_statistics(\
                                            neg_stat_test_props,
                                            response_args=params.response_args,
                                            response_significance_threshold=params.response_significance_threshold,
                                            loop_over_cells=True,
                                            repetition_keys=[varied_parameter, 'repeat'],
                                            verbose=False
                                            )

        # find responsive ROIs for this contrast (from summary stats)
        DATASET["nb_responsive_neurons_pos"].append(np.sum(pos_evokedStats['significant'], axis=0) if pos_evokedStats['significant'].size != 0 else 0)
        DATASET["nb_responsive_neurons_neg"].append(np.sum(neg_evokedStats['significant'], axis=0) if neg_evokedStats['significant'].size != 0 else 0)

        window = (ep.t>pos_stat_test_props['interval_post'][0]) & (ep.t<pos_stat_test_props['interval_post'][1])
        if pos_evokedStats['significant'].size != 0: 
            DATASET['mean_dFoF_pos'].append(np.mean(ep.dFoF[:, pos_evokedStats['significant'], :][:, :, window]))
        else :
            DATASET['mean_dFoF_pos'].append(np.nan)
        
        if neg_evokedStats['significant'].size != 0: 
            DATASET['mean_dFoF_neg'].append(np.mean(ep.dFoF[:, neg_evokedStats['significant'], :][:, :, window]))
        else :
            DATASET['mean_dFoF_neg'].append(np.nan)

        ######### 3) Get episodes' traces #########

        # INITILIAZE DICTIONARIES TO STORE RESPONSES AND BEHAVIORAL QUANTITIES
        if included_mice_pos is None:

            included_mice_pos  = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            included_mice_neg = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            dFoF_pos = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            dFoF_neg = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            deconvolved_pos = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            deconvolved_neg = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}

        for state, state_filter in zip(states_names, states_filters):

            if (pos_evokedStats['significant'].size != 0) and \
                (np.sum(pos_evokedStats['significant'], axis=0)>=params.NMIN_ROIS) and \
                    (np.sum(state_filter) >= params.NMIN_EPISODES):
                
                dFoF_pos[f"{virus}-{state}"].append(
                    ep.dFoF[state_filter][:, pos_evokedStats['significant'], :])
                
                
                deconvolved_pos[f"{virus}-{state}"].append(
                    ep.Deconvolved[state_filter][:, pos_evokedStats['significant'], :])
                
                included_mice_pos[f"{virus}-{state}"].append(DATASET['subjects'][i])

            else:
                print("cond: %s-%s -> [XX] response not included (%i ROIs, %i eps)" % 
                      (virus,state, np.sum(pos_evokedStats['significant']), np.sum(state_filter)))
            
            if (neg_evokedStats['significant'].size != 0) and \
                (np.sum(neg_evokedStats['significant'], axis=0)>=params.NMIN_ROIS) and \
                    (np.sum(state_filter) >= params.NMIN_EPISODES):
                
                dFoF_neg[f"{virus}-{state}"].append(
                    ep.dFoF[state_filter][:, neg_evokedStats['significant'], :])
                
                deconvolved_neg[f"{virus}-{state}"].append(
                    ep.Deconvolved[state_filter][:, neg_evokedStats['significant'], :])
                
                included_mice_neg[f"{virus}-{state}"].append(DATASET['subjects'][i])

            else:
                print("cond: %s-%s -> [XX] response not included (%i ROIs, %i eps)" % 
                      (virus,state, np.sum(neg_evokedStats['significant']), np.sum(state_filter)))

        ###############################

    else :
        print("session %s has no ROIs, excluded from analysis" % filename)
        DATASET["nb_responsive_neurons"].append(0)
        DATASET['mean_dFoF'].append(np.nan)
    
    print('')
                
dFoF_pos = tools.remove_empty_sessions(dFoF_pos)
deconvolved_pos = tools.remove_empty_sessions(deconvolved_pos)
dFoF_neg = tools.remove_empty_sessions(dFoF_neg)
deconvolved_neg = tools.remove_empty_sessions(deconvolved_neg)
#%% BUILD DATAFRAME
import pandas as pd 

DATASET['session'] = np.array([os.path.basename(file)[:-4] for file in DATASET['files']])

df = pd.DataFrame(DATASET).drop(columns=['protocol_ids', 'protocols', 'files']).set_index('session')
df['perc_resp_pos'] = df['nb_responsive_neurons_pos'] / df['nb_neurons'] * 100
df['perc_resp_neg'] = df['nb_responsive_neurons_neg'] / df['nb_neurons'] * 100

excel_filename = 'natimg_summary_data_' + 'NDNF' + '.xlsx'
df.to_excel(os.path.join(savepath_data, excel_filename))

df

#%% 1.a Averaged dF/F0 over positively responsive ROIs and over episodes across virus and behavioral states (std over sessions)

fig, AX = pt_fcts.plot_average_response(ep.t, dFoF_pos, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice_pos, params.NMIN_SESSIONS)
fig.tight_layout()
fig.set_size_inches(3, 4)

firgurename = 'natimg_average_pos_res_dfof_'+ 'NDNF' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg')

#%% 1.b Averaged dF/F0 with baseline substracted over positively responsive ROIs and over episodes across virus and behavioral states (std over sessions)
baselineCond = (ep.t>pos_stat_test_props['interval_pre'][0]) & (ep.t<pos_stat_test_props['interval_pre'][1])
fig, AX = pt_fcts.plot_average_response(ep.t, dFoF_pos, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice_pos, params.NMIN_SESSIONS, 
                                        baselineSubtraction=True, baselineCond=baselineCond)

fig.tight_layout()
fig.set_size_inches(3, 4)

firgurename = 'natimg_average_pos_res_dfof_bsl_sub_'+ 'NDNF' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg')

#%% 1.c Averaged deconvolved trace over positively responsive ROIs and over episodes across virus and behavioral states (std over sessions)

fig, AX = pt_fcts.plot_average_response(ep.t, deconvolved_pos, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice_pos, params.NMIN_SESSIONS)
AX[0].set_ylabel('deconvolved')
fig.tight_layout()
fig.set_size_inches(3, 4)

firgurename = 'natimg_average_pos_res_deconvolved_' + 'NDNF' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg')

#%% 2.a Averaged dF/F0 over negatively responsive ROIs and over episodes across virus and behavioral states (std over sessions)

fig, AX = pt_fcts.plot_average_response(ep.t, dFoF_neg, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice_neg, params.NMIN_SESSIONS)
fig.tight_layout()
fig.set_size_inches(3, 4)

firgurename = 'natimg_average_neg_res_dfof_' + 'NDNF' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg')

#%% 2.b Averaged dF/F0 with baseline substracted over negatively responsive ROIs and over episodes across virus and behavioral states (std over sessions)
baselineCond = (ep.t>pos_stat_test_props['interval_pre'][0]) & (ep.t<pos_stat_test_props['interval_pre'][1])

fig, AX = pt_fcts.plot_average_response(ep.t, dFoF_neg, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice_neg, params.NMIN_SESSIONS, 
                                        baselineSubtraction=True, baselineCond=baselineCond)

fig.tight_layout()
fig.set_size_inches(3, 4)

firgurename = 'natimg_average_neg_res_dfof_bsl_sub_' + 'NDNF' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg')

#%% 2.c Averaged deconvolved trace over negatively responsive ROIs and over episodes across virus and behavioral states (std over sessions)

fig, AX = pt_fcts.plot_average_response(ep.t, deconvolved_neg, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice_neg, params.NMIN_SESSIONS)
AX[0].set_ylabel('deconvolved')
fig.tight_layout()
fig.set_size_inches(3, 4)

firgurename = 'natimg_average_neg_res_deconvolved_' + 'NDNF' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg')

#%% 3. Pie chart of responsive neurons
percentages_pos = {}
percentages_neg = {}

for v in viruses:
    percentages_pos[v] = (df[df['virus'] == v]['perc_resp_pos']).to_list()
    percentages_neg[v] = (df[df['virus'] == v]['perc_resp_neg']).to_list()

fig, AX = pt_fcts.pie_chart_responsive_neurons_pos_neg(percentages_pos, percentages_neg, viruses, [])

firgurename = 'natimg_pie_charts_' + 'NDNF' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg')

#%% 4. Rasterplots with all reponsive ROIs sorted by mean 

def average_and_concatenate_sessions(arr_to_conc_1, arr_to_conc_2, viruses, states_names):
    concatenated = {}
    for v, state in product(viruses, states_names):
        key = f"{v}-{state}"
        concatenated[key] = np.array([], dtype=float).reshape(0, arr_to_conc_1[key][0].shape[2])
        for i in range(len(arr_to_conc_1[key])):
            a = np.mean(arr_to_conc_1[key][i], axis=0)
            concatenated[key] = np.concatenate((concatenated[key], a), axis=0)
        for i in range(len(arr_to_conc_2[key])):
            a = np.mean(arr_to_conc_2[key][i], axis=0)
            concatenated[key] = np.concatenate((concatenated[key], a), axis=0)
    return concatenated

dFoF_means = average_and_concatenate_sessions(dFoF_pos, dFoF_neg, viruses, states_names)

baselineCond = (ep.t>pos_stat_test_props['interval_pre'][0]) & (ep.t<pos_stat_test_props['interval_pre'][1])
response_window = (ep.t>pos_stat_test_props['interval_post'][0]) & (ep.t<pos_stat_test_props['interval_post'][1])

figurename = 'natimg_rastermap_res_dfof_bsl_sub_'+ 'NDNF' + '.svg'
fig, AX = pt_fcts.plot_rastermap(dFoF_means, ep, viruses, state_cond='all', 
                                 baselineSubtraction=True, baselineCond=baselineCond, 
                                 sort_fcts_options=dict(response_window=response_window))

fig.savefig(os.path.join(savepath_fig, figurename), transparent=True, format='svg')

figurename = 'natimg_rastermap_res_dfof_'+ 'NDNF' + '.svg'
fig, AX = pt_fcts.plot_rastermap(dFoF_means, ep, viruses, state_cond='all', 
                                 sort_fcts_options=dict(response_window=response_window))

fig.savefig(os.path.join(savepath_fig, figurename), transparent=True, format='svg')