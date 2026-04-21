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

# %%

# TO LOOP OVER NWB FILES WITH VISUAL STIMULUS --- DRIFITING GRATING ---  multisession
pname = "moving-dots"
neuron = 'NDNF-cond-CB1'  # 'PN-cond-NDNF-CB1' or 'NDNF-cond-CB1'

folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna',
                        neuron, 'NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol=pname)

# %% COMPUTATION

viruses = ['sgRosa', 'sgCnr1']
state_metric = 'speed' # 'speed' or 'pupil' or 'speed & pupil'
varied_parameter = 'speed'

ep_props = dict(quantities=['dFoF', 'running_speed'],
                prestim_duration=1.,
                dt_sampling=params.dt_sampling)

included_mice = None
means = None

run_means = None
pos_means, neg_means = None, None

# INITILIAZE DICTIONARIES TO STORE RESPONSES, PERCENTAGES AND BEHAVIORAL QUANTITIES
percentages = {v : [] for v in viruses}
pos_percentages = {v : [] for v in viruses}
neg_percentages = {v : [] for v in viruses}

# LOOP OVER SESSIONS

for i, filename in enumerate(DATASET['files']):
    
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    
    print(i+1, '--', filename, '--', data.nROIs)

    # determine virus        
    if 'sgRosa' in data.nwbfile.virus:
        virus = 'sgRosa'
    elif 'sgCnr1' in data.nwbfile.virus:
        virus = 'sgCnr1'
    else :
        raise ValueError("Virus not identified in session %s" % filename)\
    
    #date condition
    #if data.nwbfile.session_start_time.date() < datetime.date(2026, 4, 1) and virus == 'sgRosa':
    #if data.nwbfile.session_start_time.date() > datetime.date(2026, 4, 1):
    if True:

        if 'dFoF' in ep_props['quantities']:
            data.build_dFoF(**params.dFoF_options, verbose=True)
                
        if data.nROIs>0:

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
                                                    params.stat_test_props,
                                                    response_args=params.response_args,
                                                    response_significance_threshold=params.response_significance_threshold,
                                                    loop_over_cells=True,
                                                    repetition_keys=['repeat'],
                                                    verbose=False
                                                    )
                
                pos_evokedStats = ep.pre_post_statistics(\
                                                        params.pos_stat_test_props,
                                                        response_args=params.response_args,
                                                        response_significance_threshold=params.response_significance_threshold,
                                                        loop_over_cells=True,
                                                        repetition_keys=['repeat'],
                                                        verbose=False
                                                        )
                
                neg_evokedStats = ep.pre_post_statistics(\
                                                        params.neg_stat_test_props,
                                                        response_args=params.response_args,
                                                        response_significance_threshold=params.response_significance_threshold,
                                                        loop_over_cells=True,
                                                        repetition_keys=['repeat'],
                                                        verbose=False
                                                        )

                # find responsive ROIs for this contrast (from summary stats)
                percentages[virus].append(np.sum(evokedStats['significant'], axis=0)/evokedStats['significant'].shape[0]*100)
                pos_percentages[virus].append(np.sum(pos_evokedStats['significant'], axis=0)/pos_evokedStats['significant'].shape[0]*100)
                neg_percentages[virus].append(np.sum(neg_evokedStats['significant'], axis=0)/neg_evokedStats['significant'].shape[0]*100)

                ######### 3) Fill dictionnaries with responsive neurons data #########

                # INITILIAZE DICTIONARIES TO STORE RESPONSES AND BEHAVIORAL QUANTITIES
                if included_mice is None:
                    vparam_values = ep.varied_parameters[varied_parameter]
                    included_mice  = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}
                    means = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}
                    run_means = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}
                    pos_means = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}
                    neg_means = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}
                
                for k, vparam in enumerate(ep.varied_parameters[varied_parameter]):
                        
                    vparam_cond = (getattr(ep, varied_parameter)==vparam)

                    for state, state_filter in zip(states_names, states_filters):

                        if (np.sum(evokedStats['significant'][:, k], axis=0)>=params.NMIN_ROIS) and \
                                (np.sum(vparam_cond & state_filter) >= params.NMIN_EPISODES):
                            
                            print("cond: %s-%s-%s -> included %i ROIs and %i episodes" % (virus,state,round(vparam), np.sum(evokedStats['significant'][:, k]), np.sum(vparam_cond & state_filter)))
                            print("    -> %i ROIs out of %i ROIs are responsive" % (np.sum(evokedStats['significant'][:, k], axis=0), evokedStats['significant'].shape[0]))
                            
                            means[f"{virus}-{state}-{vparam}"].append(
                                ep.dFoF[vparam_cond & state_filter, :, :][:, evokedStats['significant'][:, k], :])
                            
                            pos_means[f"{virus}-{state}-{vparam}"].append(
                                ep.dFoF[vparam_cond & state_filter, :, :][:, pos_evokedStats['significant'][:, k], :])
                                                    
                            neg_means[f"{virus}-{state}-{vparam}"].append(
                                ep.dFoF[vparam_cond & state_filter, :, :][:, neg_evokedStats['significant'][:, k], :])
                            
                            included_mice[f"{virus}-{state}-{vparam}"].append(DATASET['subjects'][i])
                            
                        else:
                            print("cond: %s-%s-%s -> [XX] response not included (%i ROIs, %i eps)" % (virus,state,round(vparam), np.sum(evokedStats['significant'][:, k]), np.sum(vparam_cond & state_filter)))
        else :
            print("session %s has no ROIs, excluded from analysis" % filename)

    print('')
                
means = tools.remove_empty_sessions(means)
pos_means = tools.remove_empty_sessions(pos_means)
neg_means = tools.remove_empty_sessions(neg_means)

# now "means" is a list (over sessions) 
#    of responses of shape (episodes, responsiveROIs, time)

#%% Averaged dF/F0 over responsive ROIs and over episodes across virus and behavioral states (std over sessions)

fig, AX = pt_fcts.plot_average_response(ep.t, means, 
                                        viruses, states_names, vparam_values, varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS)

#%% Averaged dF/F0 with baseline substracted over responsive ROIs and over episodes across virus and behavioral states (std over sessions)
baselineCond = (ep.t>-0.1) & (ep.t<0)

fig, AX = pt_fcts.plot_average_response(ep.t, means, 
                                        viruses, states_names, vparam_values, varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS, 
                                        baselineSubtraction=True, baselineCond=baselineCond)

#%% Pie chart of responsive neurons

fig, AX = pt_fcts.pie_chart_responsive_neurons(percentages, viruses, vparam_values)

#%% POSITIVE RESPONSES: Averaged dF/F0 with baseline substracted over responsive ROIs and over episodes across virus and behavioral states (std over sessions)
baselineCond = (ep.t>-0.1) & (ep.t<0)

fig, AX = pt_fcts.plot_average_response(ep.t, pos_means, 
                                        viruses, states_names, vparam_values, varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS)

#%% NEGATIVE RESPONSES: Averaged dF/F0 with baseline substracted over responsive ROIs and over episodes across virus and behavioral states (std over sessions)
baselineCond = (ep.t>-0.1) & (ep.t<0)

fig, AX = pt_fcts.plot_average_response(ep.t, neg_means, 
                                        viruses, states_names, vparam_values, varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS)

#%% Pie chart of responsive neurons separated by positive and negative responses

fig, AX = pt_fcts.pie_chart_responsive_neurons_pos_neg(pos_percentages, neg_percentages, viruses, vparam_values)

#%% calculate z-scores
import pandas as pd

def zscores(baselines, traces):
    return np.array([(traces[i, :] - np.mean(baselines, axis=1)[i]) / np.std(baselines, axis=1)[i] for i in range(traces.shape[0])])
    
zscores_session = {f'{v}-{speed}' : [] for v, speed in product(viruses, vparam_values)}
zscores_session_pos = {f'{v}-{speed}' : [] for v, speed in product(viruses, vparam_values)}
zscores_session_neg = {f'{v}-{speed}' : [] for v, speed in product(viruses, vparam_values)}

nb_rois = {f'{v}-{speed}' : [] for v, speed in product(viruses, vparam_values)}
nb_rois_pos = {f'{v}-{speed}' : [] for v, speed in product(viruses, vparam_values)}
nb_rois_neg = {f'{v}-{speed}' : [] for v, speed in product(viruses, vparam_values)}

baselineCond = (ep.t<0)

for v in ['sgCnr1', 'sgRosa']:
    
    for s in vparam_values:

        key = f'{v}-{s}'

        session_responses = [np.mean(r, axis=0) for r in means[f'{v}-still-{s}']]
        session_responses_pos = [np.mean(r, axis=0) for r in pos_means[f'{v}-still-{s}']]
        session_responses_neg = [np.mean(r, axis=0) for r in neg_means[f'{v}-still-{s}']]
        
        for i in range(len(session_responses)):
            zscores_session_i = zscores(session_responses[i][:, baselineCond], session_responses[i])
            zscores_session[key].append(np.mean(zscores_session_i, axis=0))
            nb_rois[key].append(zscores_session_i.shape[0])

        for i in range(len(session_responses_pos)):
            zscores_session_i = zscores(session_responses_pos[i][:, baselineCond], session_responses_pos[i])
            zscores_session_pos[key].append(np.mean(zscores_session_i, axis=0))
            nb_rois_pos[key].append(zscores_session_i.shape[0])

        for i in range(len(session_responses_neg)):
            zscores_session_i = zscores(session_responses_neg[i][:, baselineCond], session_responses_neg[i])
            zscores_session_neg[key].append(np.mean(zscores_session_i, axis=0))
            nb_rois_neg[key].append(zscores_session_i.shape[0])

        nb_rois[key] = np.array(nb_rois[key])
        nb_rois_pos[key] = np.array(nb_rois_pos[key])
        nb_rois_neg[key] = np.array(nb_rois_neg[key])
        zscores_session[key] = np.array(zscores_session[key])
        zscores_session_pos[key] = np.array(zscores_session_pos[key])
        zscores_session_neg[key] = np.array(zscores_session_neg[key])

d = {'time': ep.t}
for v in viruses:
    for s in vparam_values:
        d.update({f'{v}-{s}_avg': np.mean(zscores_session[f'{v}-{s}'], axis=0)})
        d.update({f'{v}-{s}_sem': sem(zscores_session[f'{v}-{s}'], axis=0)})
        d.update({f'{v}-{s}_avg_pos': np.mean(zscores_session_pos[f'{v}-{s}'], axis=0)})
        d.update({f'{v}-{s}_sem_pos': sem(zscores_session_pos[f'{v}-{s}'], axis=0)})
        d.update({f'{v}-{s}_avg_neg': np.mean(zscores_session_neg[f'{v}-{s}'], axis=0)})
        d.update({f'{v}-{s}_sem_neg': sem(zscores_session_neg[f'{v}-{s}'], axis=0)})

df = pd.DataFrame(data=d)
savepath = None
#df.to_excel(os.path.join(savepath, filename))

#%%
color_virus = {'sgRosa' : 'grey', 
               'sgCnr1': "darkred"}

fig, AX = pt.figure()
speed = 60

for k, v in enumerate(['sgCnr1']):
    
    key = f'{v}-{speed}'
    print(key)

    nb_mice = np.unique(included_mice[f'{v}-still-{s}']).shape[0]

    pt.plot(ep.t, np.mean(zscores_session[key], axis=0), 
            sy=sem(zscores_session[key], axis=0), 
            color=color_virus[v], label=v, ax=AX)
    pt.annotate(AX, 'N=%i (%i mice, %i rois)' % (len(zscores_session[key]), 
                                                nb_mice,
                                                nb_rois[key].sum()) +k*'\n',
                                                
                color=color_virus[v], **dict(xy=(0.05,1), ha='left', fontsize=4))

#%%

fig, AX = pt.figure()
speed = 60

for k, v in enumerate(['sgCnr1']):
    
    key = f'{v}-{speed}'
    print(key)

    nb_mice = np.unique(included_mice[f'{v}-still-{s}']).shape[0]

    pt.plot(ep.t, np.mean(zscores_session_pos[key], axis=0), 
            sy=sem(zscores_session_pos[key], axis=0), 
            color=color_virus[v], label=v, ax=AX)
    pt.annotate(AX, 'N=%i (%i mice, %i rois)' % (len(zscores_session_pos[key]), 
                                                nb_mice,
                                                nb_rois_pos[key].sum()) +k*'\n',
                                                
                color=color_virus[v], **dict(xy=(0.05,1), ha='left', fontsize=4))
    
#%%

fig, AX = pt.figure()
speed = 60

for k, v in enumerate(['sgCnr1']):
    
    key = f'{v}-{speed}'
    print(key)

    nb_mice = np.unique(included_mice[f'{v}-still-{s}']).shape[0]

    pt.plot(ep.t, np.mean(zscores_session_neg[key], axis=0), 
            sy=sem(zscores_session_neg[key], axis=0), 
            color=color_virus[v], label=v, ax=AX)
    pt.annotate(AX, 'N=%i (%i mice, %i rois)' % (len(zscores_session_neg[key]), 
                                                nb_mice,
                                                nb_rois_neg[key].sum()) +k*'\n',
                                                
                color=color_virus[v], **dict(xy=(0.05,1), ha='left', fontsize=4))