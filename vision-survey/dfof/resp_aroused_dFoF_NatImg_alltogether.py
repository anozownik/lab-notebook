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
neuron = 'NDNF-cond-CB1'  # 'PN_cond-NDNF-CB1_WT-vs-KD' or 'NDNF-cond-CB1'

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
                                                    repetition_keys=[varied_parameter, 'repeat'],
                                                    verbose=False
                                                    )
                
                pos_evokedStats = ep.pre_post_statistics(\
                                                        params.pos_stat_test_props,
                                                        response_args=params.response_args,
                                                        response_significance_threshold=params.response_significance_threshold,
                                                        loop_over_cells=True,
                                                        repetition_keys=[varied_parameter, 'repeat'],
                                                        verbose=False
                                                        )
                
                neg_evokedStats = ep.pre_post_statistics(\
                                                        params.neg_stat_test_props,
                                                        response_args=params.response_args,
                                                        response_significance_threshold=params.response_significance_threshold,
                                                        loop_over_cells=True,
                                                        repetition_keys=[varied_parameter, 'repeat'],
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
                    included_mice  = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
                    means = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
                    run_means = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
                    pos_means = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
                    neg_means = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}

                    means_rel = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
                    pos_means_rel = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
                    neg_means_rel = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
                    included_mice_rel  = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}

                for state, state_filter in zip(states_names, states_filters):

                    if (np.sum(evokedStats['significant'], axis=0)>=params.NMIN_ROIS) and \
                            (np.sum(state_filter) >= params.NMIN_EPISODES):
                        
                        print("cond: %s-%s -> included %i ROIs and %i episodes" % (virus,state, np.sum(evokedStats['significant']), np.sum(state_filter)))
                        print("    -> %i ROIs out of %i ROIs are responsive" % (np.sum(evokedStats['significant'], axis=0), evokedStats['significant'].shape[0]))
                        
                        means[f"{virus}-{state}"].append(
                            ep.dFoF[state_filter][:, evokedStats['significant'], :])
                        
                        pos_means[f"{virus}-{state}"].append(
                            ep.dFoF[state_filter][:, pos_evokedStats['significant'], :])
                                                
                        neg_means[f"{virus}-{state}"].append(
                            ep.dFoF[state_filter][:, neg_evokedStats['significant'], :])
                        
                        included_mice[f"{virus}-{state}"].append(DATASET['subjects'][i])

                        run_means[f"{virus}-{state}"].append(ep.running_speed[state_filter, :])
                        
                    else:
                        print("cond: %s-%s -> [XX] response not included (%i ROIs, %i eps)" % (virus,state, np.sum(evokedStats['significant']), np.sum(state_filter)))

                ######### 4) Run reliability #########
                if COMPUTE_RELIABILITY :
                    
                    ep.dFoF = ep.dFoF[:, evokedStats['significant'], :]
                    summary_reliability = ep.reliability(
                                                response_args=params.response_args,
                                                response_significance_threshold=0.01,
                                                loop_over_cells=True,
                                                repetition_keys=[varied_parameter, 'repeat'],
                                                verbose=True)
                    reliability[virus].append({key: summary_reliability[key] for key in ['r', 'pval']})

                    for state, state_filter in zip(states_names, states_filters):

                        if (np.sum(summary_reliability['significant'], axis=0)>=params.NMIN_ROIS) and \
                                (np.sum(state_filter) >= params.NMIN_EPISODES):
                            
                            print("cond: %s-%s -> included %i ROIs and %i episodes" % (virus,state, np.sum(summary_reliability['significant']), np.sum(state_filter)))
                            print("    -> %i ROIs out of %i ROIs are responsive" % (np.sum(summary_reliability['significant'], axis=0), summary_reliability['significant'].shape[0]))
                            
                            session_response = ep.dFoF[state_filter][:, summary_reliability['significant'], :]
                            average_trial_centered = (np.mean(session_response, axis=0).T - np.mean(session_response[:, :, ep.t<0], axis=(0, 2))).T

                            means_rel[f"{virus}-{state}"].append(ep.dFoF[state_filter][:, summary_reliability['significant'], :])
                            
                            pos = np.array([auc(np.arange(0, session_response.shape[2]), average_trial_centered[i]) for i in range(session_response.shape[1])]) >= 0

                            if np.sum(pos) >= params.NMIN_ROIS:
                                pos_means_rel[f"{virus}-{state}"].append(
                                    ep.dFoF[state_filter][:, summary_reliability['significant'], :][:, pos, :])

                            if np.sum(~pos) >= params.NMIN_ROIS:                 
                                neg_means_rel[f"{virus}-{state}"].append(
                                    ep.dFoF[state_filter][:, summary_reliability['significant'], :][:, ~pos, :])
                            
                            included_mice_rel[f"{virus}-{state}"].append(DATASET['subjects'][i])

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

np.save(os.path.join(tempfile.tempdir, 'means.npy'), means)
np.save(os.path.join(tempfile.tempdir, 'run_means.npy'), run_means)
np.save(os.path.join(tempfile.tempdir, 'pos_means.npy'), pos_means)
np.save(os.path.join(tempfile.tempdir, 'neg_means.npy'), neg_means)
np.save(os.path.join(tempfile.tempdir, 'included_mice.npy'), included_mice)
np.save(os.path.join(tempfile.tempdir, 'percentages.npy'), percentages)
np.save(os.path.join(tempfile.tempdir, 'pos_percentages.npy'), pos_percentages)
np.save(os.path.join(tempfile.tempdir, 'neg_percentages.npy'), neg_percentages)
#np.save(os.path.join(tempfile.tempdir, 'reliability.npy'), reliability)

#%% Averaged dF/F0 over responsive ROIs and over episodes across virus and behavioral states (std over sessions)

firgurename = 'dfof_beh_mod_all_natimg'+ neuron + '.svg'
figurepath = '/Users/macbookair/work/Figures/Natural_images/'

fig, AX = pt_fcts.plot_average_response(ep.t, means, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS)
#plt.savefig(os.path.join(figurepath+firgurename), transparent=True, format='svg')

#%% Averaged dF/F0 with baseline substracted over responsive ROIs and over episodes across virus and behavioral states (std over sessions)
baselineCond = (ep.t>-0.1) & (ep.t<0)

fig, AX = pt_fcts.plot_average_response(ep.t, means, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS, 
                                        baselineSubtraction=True, baselineCond=baselineCond)

#%% Pie chart of responsive neurons

firgurename = 'pie_all_natimg'+ neuron + '.svg'
fig, AX = pt_fcts.pie_chart_responsive_neurons(percentages, viruses, [])
#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')

#%% POSITIVE RESPONSES: Averaged dF/F0 with baseline substracted over responsive ROIs and over episodes across virus and behavioral states (std over sessions)
baselineCond = (ep.t>-0.1) & (ep.t<0)

firgurename = 'dfof_beh_mod_onlypos_natimg'+ neuron + '.svg'
figurepath = '/Users/macbookair/work/Figures/Natural_images/'

fig, AX = pt_fcts.plot_average_response(ep.t, pos_means, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS, 
                                        baselineSubtraction=True, baselineCond=baselineCond)
#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')

#%% NEGATIVE RESPONSES: Averaged dF/F0 with baseline substracted over responsive ROIs and over episodes across virus and behavioral states (std over sessions)

firgurename = 'dfof_beh_mod_onlyneg_natimg_'+ neuron + '.svg'

fig, AX = pt_fcts.plot_average_response(ep.t, neg_means, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS)
#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')
#%% Pie chart of responsive neurons separated by positive and negative responses
firgurename = 'pie_resp_neg_natimg_'+ neuron + '.svg'

fig, AX = pt_fcts.pie_chart_responsive_neurons_pos_neg(pos_percentages, neg_percentages, viruses, [])
#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')

#%% Speed across virus and behavioral states

fig, AX = pt_fcts.plot_average_behavior(ep.t, run_means, viruses, states_names, params.NMIN_SESSIONS, ylabel='speed(cm/s)')


#%% Rastermap session k: averaged dF/F0 over responsive ROIs per episode across virus 
cond = 'all' #behavioral condition to plot
session_id = {'sgRosa' : 0, 'sgCnr1': 0} #session index to plot

fig, AX = pt_fcts.rastermap_session(session_id, ep, means, viruses, state_cond=cond)

#%% Rastermap session k: averaged dF/F0 (baseline substracted) over responsive ROIs per episode across virus
cond = 'all' #behavioral condition to plot
baselineCond = (ep.t>-1.9) & (ep.t<0)
session_id = {'sgRosa' : 0, 'sgCnr1': 0} #session index to plot

fig, AX = pt_fcts.rastermap_session(session_id, ep, means, viruses, state_cond=cond, baselineSubtraction=True, baselineCond=baselineCond)

#%%
color_virus = {'sgRosa' : 'grey', 
               'sgCnr1': "darkred"}

def plot_dist_reliability(reliability, viruses, plot_type='violin',
                          only_significant=True, significant_threshold=0.01, savepath=None):

    r_values = {virus: [] for virus in viruses}
    for v in viruses:

        for i in range(len(reliability[v])):

            if only_significant:
                significant = reliability[v][i]['pval'] <= significant_threshold
                r_values[v].append(reliability[v][i]['r'][significant])
            else :
                r_values[v].append(reliability[v][i]['r'])
        
        r_values[v] = np.concatenate(r_values[v])

    if plot_type == 'violin':

        fig, AX = plt.subplots(1, 1, figsize=(1.5, 1.5))

        for k, virus in enumerate(viruses):
            pt.violin(r_values[virus], x=k*1, color=color_virus[virus], ax=AX)

        pt.set_plot(AX, ['left'])

    elif plot_type == 'hist' or plot_type == 'histogram':

        fig, AX = pt.figure(axes=(len(viruses), 1))

        for k, virus in enumerate(viruses):
            ymax = np.histogram(r_values[virus], bins=20)[0].max()

            AX[k].hist(r_values[virus], bins=20, color=color_virus[virus])
            AX[k].vlines(np.mean(r_values[virus]), 0, ymax, color='black', linewidth=0.5, linestyle='dashed')
            AX[k].annotate('mean=%.2f' % np.mean(r_values[virus]), xy=(np.mean(r_values[virus])-0.01, ymax), 
                        ha='right', fontsize=4)
            AX[k].set_title(virus+ ' (n=%d)' % len(r_values[virus]))
            AX[k].set_xlim(round(np.min(r_values[virus])*10)*0.1 - 0.1, 1)
    
    else :
        raise ValueError('plot_type value not recognized, should be "violin" or "hist"')

    return fig, AX

plot_dist_reliability(reliability, viruses, plot_type='hist')
plot_dist_reliability(reliability, viruses)
plot_dist_reliability(reliability, viruses, plot_type='hist', only_significant=False)
plot_dist_reliability(reliability, viruses, only_significant=False)

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