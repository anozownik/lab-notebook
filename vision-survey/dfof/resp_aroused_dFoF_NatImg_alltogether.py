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
COMPUTE_RELIABILITY = False
stat_test_props = params.stat_test_props
stat_test_props['interval_post'] = [0.5, 1.5]

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
                                        None, params.NMIN_SESSIONS, 
                                        baselineSubtraction=True, baselineCond=baselineCond)
#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')

#%% NEGATIVE RESPONSES: Averaged dF/F0 with baseline substracted over responsive ROIs and over episodes across virus and behavioral states (std over sessions)

firgurename = 'dfof_beh_mod_onlyneg_natimg_'+ neuron + '.svg'

fig, AX = pt_fcts.plot_average_response(ep.t, neg_means, 
                                        viruses, states_names, [], varied_parameter, 
                                        None, params.NMIN_SESSIONS)
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

pt_fcts.plot_dist_reliability(reliability, viruses, plot_type='hist')
pt_fcts.plot_dist_reliability(reliability, viruses)
pt_fcts.plot_dist_reliability(reliability, viruses, plot_type='hist', only_significant=False)
pt_fcts.plot_dist_reliability(reliability, viruses, only_significant=False)

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


#%% CV factor formula 1
resp_win = (ep.t > 0.2) & (ep.t < 1.5)
CV_coef = {v: [] for v in viruses}
fano_factor = {v: [] for v in viruses}


for v in viruses:

    nb_neurons = 0
    print(f"Processing {v}  ...")

    for i in range(len(means[f"{v}-all"])):

        nb_neurons += means[f"{v}-all"][i].shape[1]

        mean_trial = means[f"{v}-all"][i][:, :, resp_win].mean(axis=2)
        cv = np.std(mean_trial, axis=0) / np.mean(mean_trial, axis=0)
        fano = np.var(mean_trial, axis=0) / np.mean(mean_trial, axis=0)

        CV_coef[f"{v}"].append(cv)
        fano_factor[f"{v}"].append(fano)
    CV_coef[f"{v}"] = np.concatenate(CV_coef[f"{v}"])
    fano_factor[f"{v}"] = np.concatenate(fano_factor[f"{v}"])
    print(nb_neurons)

#%% CV factor formula 2
resp_win = (ep.t > 0.2) & (ep.t < 1.5)
CV_coef = {v : [] for v in viruses}
fano_factor = {v : [] for v in viruses}

for v in viruses:

    nb_neurons = 0
    print(f"Processing {v} ...")

    for i in range(len(means[f"{v}-all"])):

        nb_neurons += means[f"{v}-all"][i].shape[1]

        trials = means[f"{v}-all"][i][:, :, resp_win]
        cv = np.std(trials, axis=0).mean(axis=1) / np.mean(trials, axis=0).mean(axis=1)
        fano = np.var(trials, axis=0).mean(axis=1) /  np.mean(trials, axis=0).mean(axis=1)

        CV_coef[v].append(cv)
        fano_factor[v].append(fano)
    CV_coef[v] = np.concatenate(CV_coef[v])
    fano_factor[v] = np.concatenate(fano_factor[v])
    print(nb_neurons)


#%% plot CV
color_virus = {'sgRosa' : 'grey', 
               'sgCnr1': "darkred"}
xdispersion = 0.5

fig, AX = pt.figure(axes=(1, 1))


for k, v in enumerate(viruses):
    pt.bar([np.mean(CV_coef[f"{v}"])], x=[k], sy=np.std(CV_coef[f"{v}"]), ax=AX, color=color_virus[v])
    """ pt.scatter(np.ones(len(CV_coef[f"{v}"]))*k+np.random.rand(len(CV_coef[f"{v}"]))*xdispersion-xdispersion/2, 
                CV_coef[f"{v}"], ax=AX[0], color='black', alpha=0.5) """

pt.set_plot(AX, ylabel='CV')
AX.set_xticks([0, 1], viruses, fontsize=7)

fig, AX = pt.figure(axes=(1, 1))


for k, v in enumerate(viruses):
    pt.violin(CV_coef[f"{v}"], x=k*1, ax=AX, color=color_virus[v])

pt.set_plot(AX, ylabel='CV')
AX.set_xticks([0, 1], viruses, fontsize=7)

#%% plot fano factor
color_virus = {'sgRosa' : 'grey', 
               'sgCnr1': "darkred"}
xdispersion = 0.5

fig, AX = pt.figure(axes=(1, 1))


for k, v in enumerate(viruses):
    pt.bar([np.mean(fano_factor[f"{v}"])], x=[k], sy=np.std(fano_factor[f"{v}"]), ax=AX, color=color_virus[v])
    """ pt.scatter(np.ones(len(fano_factor[f"{v}"]))*k+np.random.rand(len(fano_factor[f"{v}"]))*xdispersion-xdispersion/2, 
                fano_factor[f"{v}"], ax=AX, color='black', alpha=0.5) """

pt.set_plot(AX, ylabel='Fano Factor')
AX.set_xticks([0, 1], viruses, fontsize=7)

fig, AX = pt.figure(axes=(1, 1))


for k, v in enumerate(viruses):
    pt.violin(fano_factor[f"{v}"], x=k*1, ax=AX, color=color_virus[v])

pt.set_plot(AX, ylabel='Fano Factor')
AX.set_xticks([0, 1], viruses, fontsize=7)

# %% plot gain of arousal modulation by fitting a multiplicative factor to the still condition to explain the run condition
from scipy.optimize import minimize
from scipy import stats
from itertools import compress

# fitting window
t_cond = ep.t > -0.1

# baseline window
baselineCond = (ep.t > -0.1) & (ep.t < 0)

gains = {}

for virus in ['sgRosa', 'sgCnr1']:

    gains[virus] = {}

    gains[virus] = []

    run_sessions   = means[f'{virus}-run']
    still_sessions = means[f'{virus}-still']
    run_mice = np.array(included_mice[f'{virus}-run']) 
    still_mice = np.array(included_mice[f'{virus}-still'])

    for mouse in np.unique(run_mice):

        run_sessions_mouse = list(compress(means[f'{virus}-run'], run_mice==mouse))
        still_sessions_mouse = list(compress(means[f'{virus}-still'], still_mice==mouse))

        run_sessions_mouse_average = []
        still_sessions_mouse_average = []

        # loop over sessions
        for run_data, still_data in zip(run_sessions_mouse, still_sessions_mouse):

            # average across episodes
            run_trace = np.mean(run_data, axis=(0, 1))
            still_trace = np.mean(still_data, axis=(0, 1))

            run_sessions_mouse_average.append(run_trace)
            still_sessions_mouse_average.append(still_trace)

        run_sessions_mouse_average = np.array(run_sessions_mouse_average)
        still_sessions_mouse_average = np.array(still_sessions_mouse_average)
        
        run_trace = np.mean(run_sessions_mouse_average, axis=0)
        still_trace = np.mean(still_sessions_mouse_average, axis=0)
        
        # baseline subtraction
        run_trace -= np.mean(run_trace[baselineCond])
        still_trace -= np.mean(still_trace[baselineCond])

        # skip flat traces
        if np.std(still_trace[t_cond]) < 1e-10:
            continue

        # multiplicative gain fit
        def to_minimize(x):

            return np.sum(
                (x[0] * still_trace[t_cond]
                    - run_trace[t_cond])**2
            )

        res = minimize(to_minimize, [1.0])

        gains[virus].append(res.x[0])
# %%
fig, AX = pt.figure(axes=(1, 1))

rosa = gains['sgRosa']
cnr1 = gains['sgCnr1']

AX.bar([0,1],[np.mean(rosa), np.mean(cnr1)],yerr=[stats.sem(rosa), stats.sem(cnr1)],color=['grey','darkred'])
AX.scatter(np.ones(len(rosa))*0+np.random.rand(len(rosa))*0.2-0.1, rosa, color='black', alpha=0.5)
AX.scatter(np.ones(len(cnr1))*1+np.random.rand(len(cnr1))*0.2-0.1, cnr1, color='black', alpha=0.5)

AX.set_xticks([0,1])
AX.set_xticklabels(['sgRosa','sgCnr1'])
AX.set_ylabel('gain')
pt.annotate(AX, f'WT: N={len(rosa)} mice\nKD: N={len(cnr1)} mice', xy=(0.5, 1), xycoords='axes fraction', ha='right', fontsize=4)
plt.tight_layout()

#%%

from sklearn.linear_model import LinearRegression
color_virus = {'sgRosa' : 'grey', 
               'sgCnr1': "darkred"}


nb_neurons = {}
for v in viruses:
    nb_neurons[v] = [el.shape[1] for el in means[f"{v}-all"]]

model_wt = LinearRegression()
model_wt.fit(np.array(nb_neurons['sgRosa']).reshape(-1,1), percentages['sgRosa'])

model_kt = LinearRegression()
model_kt.fit(np.array(nb_neurons['sgCnr1']).reshape(-1,1), percentages['sgCnr1'])


fig, ax = pt.figure(figsize=(3.,3.))
for i, (virus, model) in enumerate(zip(viruses, [model_wt, model_kt])):
    ax.scatter(nb_neurons[virus], percentages[virus], color=color_virus[virus], alpha=0.5)

    """ x = np.linspace(np.min(nb_neurons[virus]), np.max(nb_neurons[virus]), 100)
    y =  model.coef_ * x + model.intercept_
    ax.plot(x, y, color='black') """

pt.set_plot(ax, xlabel='Nb neurons', ylabel='% of resp.')