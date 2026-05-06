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
neuron = 'PN_cond-NDNF-CB1_WT-vs-KD'  # 'PN_cond-NDNF-CB1_WT-vs-KD' or 'NDNF-cond-CB1'

folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna',
                        neuron, 'NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol=pname)

# %% COMPUTATION

viruses = ['sgRosa', 'sgCnr1']
state_metric = 'speed' # 'speed' or 'pupil' or 'speed & pupil'
varied_parameter = 'Image-ID'
COMPUTE_RELIABILITY = True

ep_props = dict(quantities=['dFoF', 'running_speed'],
                prestim_duration=1.5,
                dt_sampling=params.dt_sampling)

included_mice = None
means = None

run_means = None
pos_means, neg_means = None, None
pos_means_rel, neg_means_rel = None, None
means_rel, included_mice_rel = None, None

# INITILIAZE DICTIONARIES TO STORE RESPONSES, PERCENTAGES AND BEHAVIORAL QUANTITIES
percentages = {v : [] for v in viruses}
pos_percentages = {v : [] for v in viruses}
neg_percentages = {v : [] for v in viruses}
reliability = {v : [] for v in viruses}

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

                    means_rel = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}
                    pos_means_rel = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}
                    neg_means_rel = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}
                    included_mice_rel  = {f"{v}-{cond}-{vparam}" : [] for v, cond, vparam in product(viruses, states_names, vparam_values)}

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
            
                ######### 4) Run reliability #########
                if COMPUTE_RELIABILITY :
                    
                    ep.dFoF = ep.dFoF[:, evokedStats['significant'].sum(axis=1) >= 1, :]
                    summary_reliability = ep.reliability(
                                                response_args=params.response_args,
                                                response_significance_threshold=0.01,
                                                loop_over_cells=True,
                                                repetition_keys=['repeat'],
                                                verbose=False)
                    reliability[virus].append({key: summary_reliability[key] for key in ['r', 'pval']})

                    for k, vparam in enumerate(ep.varied_parameters[varied_parameter]):

                        vparam_cond = (getattr(ep, varied_parameter)==vparam)

                        for state, state_filter in zip(states_names, states_filters):

                            if (np.sum(summary_reliability['significant'][:, k], axis=0)>=params.NMIN_ROIS) and \
                                    (np.sum(vparam_cond & state_filter) >= params.NMIN_EPISODES):
                                
                                print("cond: %s-%s-%s -> included %i ROIs and %i episodes" % (virus,state,round(vparam), np.sum(summary_reliability['significant'][:, k]), np.sum(vparam_cond & state_filter)))
                                print("    -> %i ROIs out of %i ROIs are responsive" % (np.sum(summary_reliability['significant'][:, k], axis=0), summary_reliability['significant'].shape[0]))
                                
                                session_response = ep.dFoF[vparam_cond & state_filter][:, summary_reliability['significant'][:, k], :]
                                average_trial_centered = (np.mean(session_response, axis=0).T - np.mean(session_response[:, :, ep.t<0], axis=(0, 2))).T

                                means_rel[f"{virus}-{state}-{vparam}"].append(ep.dFoF[vparam_cond & state_filter][:, summary_reliability['significant'][:, k], :])
                                
                                pos = np.array([auc(np.arange(0, session_response.shape[2]), average_trial_centered[i]) for i in range(session_response.shape[1])]) >= 0

                                if np.sum(pos) >= params.NMIN_ROIS:
                                    pos_means_rel[f"{virus}-{state}-{vparam}"].append(
                                        ep.dFoF[vparam_cond & state_filter][:, summary_reliability['significant'][:, k], :][:, pos, :])

                                if np.sum(~pos) >= params.NMIN_ROIS:                 
                                    neg_means_rel[f"{virus}-{state}-{vparam}"].append(
                                        ep.dFoF[vparam_cond & state_filter][:, summary_reliability['significant'][:, k], :][:, ~pos, :])
                                
                                included_mice_rel[f"{virus}-{state}-{vparam}"].append(DATASET['subjects'][i])

                            else:
                                print("cond: %s-%s -> [XX] response not included (%i ROIs, %i eps)" % (virus,state, np.sum(evokedStats['significant']), np.sum(state_filter)))

            else :
                print("%s is not a valid varied parameter" % varied_parameter)

        else :
            print("session %s has no ROIs, excluded from analysis" % filename)

    print('')
                
means = tools.remove_empty_sessions(means)
pos_means = tools.remove_empty_sessions(pos_means)
neg_means = tools.remove_empty_sessions(neg_means)

means_rel = tools.remove_empty_sessions(means_rel)
pos_means_rel = tools.remove_empty_sessions(pos_means_rel)
neg_means_rel = tools.remove_empty_sessions(neg_means_rel)

# now "means" is a list (over sessions) 
#    of responses of shape (episodes, responsiveROIs, time)

#%% Averaged dF/F0 over responsive ROIs and over episodes across virus and behavioral states (std over sessions)

figurepath = '/Users/macbookair/work/Figures/'
firgurename = 'dfof_natimgID_beh_mod_'+ neuron + '.svg'

fig, AX = pt_fcts.plot_average_response(ep.t, means, 
                                        viruses, states_names, vparam_values, varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS)
#plt.savefig(os.path.join(figurepath+firgurename), transparent=True, format='svg')

#%% Averaged dF/F0 with baseline substracted over responsive ROIs and over episodes across virus and behavioral states (std over sessions)
baselineCond = (ep.t>-0.1) & (ep.t<0)

fig, AX = pt_fcts.plot_average_response(ep.t, means, 
                                        viruses, states_names, vparam_values, varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS, 
                                        baselineSubtraction=True, baselineCond=baselineCond)

#%% Pie chart of responsive neurons

firgurename = 'pie_natimgID_beh_mod_'+ neuron + '.svg'
fig, AX = pt_fcts.pie_chart_responsive_neurons(percentages, viruses, vparam_values)
#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')

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


#%%

pt_fcts.plot_dist_reliability(reliability, viruses, varied_parameter=vparam_values, vparam_name=varied_parameter, plot_type='hist')
pt_fcts.plot_dist_reliability(reliability, viruses, varied_parameter=vparam_values, vparam_name=varied_parameter)
pt_fcts.plot_dist_reliability(reliability, viruses, varied_parameter=vparam_values, vparam_name=varied_parameter, plot_type='hist', only_significant=False)
pt_fcts.plot_dist_reliability(reliability, viruses, varied_parameter=vparam_values, vparam_name=varied_parameter, only_significant=False)

#%%

fig, AX = pt_fcts.plot_average_response(ep.t, neg_means_rel, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice_rel, params.NMIN_SESSIONS)

pt.set_common_ylims(AX)

#%%

fig, AX = pt_fcts.plot_average_response(ep.t, neg_means, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS)

pt.set_common_ylims(AX)

#%%
baselineCond = (ep.t<0)
fig, AX = pt_fcts.plot_average_response(ep.t, neg_means_rel, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice_rel, params.NMIN_SESSIONS, 
                                        baselineSubtraction=True, baselineCond=baselineCond)
pt.set_common_ylims(AX)