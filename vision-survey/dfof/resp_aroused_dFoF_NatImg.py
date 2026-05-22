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
stat_test_props['interval_post'] = [1, 2.8]

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
                                                    stat_test_props,
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
                                        None, params.NMIN_SESSIONS)

#%% NEGATIVE RESPONSES: Averaged dF/F0 with baseline substracted over responsive ROIs and over episodes across virus and behavioral states (std over sessions)
baselineCond = (ep.t>-0.1) & (ep.t<0)

fig, AX = pt_fcts.plot_average_response(ep.t, neg_means, 
                                        viruses, states_names, vparam_values, varied_parameter, 
                                        None, params.NMIN_SESSIONS)

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

#%% CV factor
resp_win = (ep.t > 0.2) & (ep.t < 1.5)
CV_coef = {f"{v}-{vparam}" : [] for v, vparam in product(viruses, vparam_values)}
fano_factor = {f"{v}-{vparam}" : [] for v, vparam in product(viruses, vparam_values)}

for v in viruses:
    for vparam in vparam_values:

        nb_neurons = 0
        print(f"Processing {v} - {vparam} ...")

        for i in range(len(means[f"{v}-all-{vparam}"])):

            nb_neurons += means[f"{v}-all-{vparam}"][i].shape[1]

            mean_trial = means[f"{v}-all-{vparam}"][i][:, :, resp_win].mean(axis=2)
            cv = np.std(mean_trial, axis=0) / np.mean(mean_trial, axis=0)
            fano = np.var(mean_trial, axis=0) / np.mean(mean_trial, axis=0)

            CV_coef[f"{v}-{vparam}"].append(cv)
            fano_factor[f"{v}-{vparam}"].append(fano)
        CV_coef[f"{v}-{vparam}"] = np.concatenate(CV_coef[f"{v}-{vparam}"])
        fano_factor[f"{v}-{vparam}"] = np.concatenate(fano_factor[f"{v}-{vparam}"])
        print(nb_neurons)

#%% Gain CV coefficient, fano factor

from scipy.optimize import minimize
from scipy import stats

# fitting window
t_cond = ep.t > -0.1

# baseline window
baselineCond = (ep.t > -0.1) & (ep.t < 0)

resp_win = (ep.t > 0.2) & (ep.t < 1.5)
CV_coef = {f"{v}-{vparam}" : [] for v, vparam in product(viruses, vparam_values)}
fano_factor = {f"{v}-{vparam}" : [] for v, vparam in product(viruses, vparam_values)}
gains = {f"{v}-{vparam}" : [] for v, vparam in product(viruses, vparam_values)}


for v in viruses:
    for vparam in vparam_values:

        nb_neurons = 0
        print(f"Processing {v} - {vparam} ...")

        for i in range(len(means[f"{v}-all-{vparam}"])):

            nb_neurons += means[f"{v}-all-{vparam}"][i].shape[1]

            trials = means[f"{v}-all-{vparam}"][i][:, :, resp_win]
            cv = np.std(trials, axis=0).mean(axis=1) / np.mean(trials, axis=0).mean(axis=1)
            fano = np.var(trials, axis=0).mean(axis=1) /  np.mean(trials, axis=0).mean(axis=1)

            CV_coef[f"{v}-{vparam}"].append(cv)
            fano_factor[f"{v}-{vparam}"].append(fano)
        CV_coef[f"{v}-{vparam}"] = np.concatenate(CV_coef[f"{v}-{vparam}"])
        fano_factor[f"{v}-{vparam}"] = np.concatenate(fano_factor[f"{v}-{vparam}"])
        print(nb_neurons)

#%%
color_virus = {'sgRosa' : 'grey', 
               'sgCnr1': "darkred"}
xdispersion = 0.5

fig, AX = pt.figure(axes=(5, 1))

for i, vparam in enumerate(vparam_values):

    for k, v in enumerate(viruses):
        pt.bar([np.mean(CV_coef[f"{v}-{vparam}"])], x=[k], sy=np.std(CV_coef[f"{v}-{vparam}"]), ax=AX[i], color=color_virus[v])
        """ pt.scatter(np.ones(len(CV_coef[f"{v}-{vparam}"]))*k+np.random.rand(len(CV_coef[f"{v}-{vparam}"]))*xdispersion-xdispersion/2, 
                   CV_coef[f"{v}-{vparam}"], ax=AX[i], color='black', alpha=0.5) """

    pt.annotate(AX[i], '%s=%.0f ' % (varied_parameter, vparam), (0.5, 1.1), ha='center')
    pt.set_plot(AX[i], ylabel='CV' if i==0 else '')
    AX[i].set_xticks([0, 1], viruses, fontsize=7)

fig, AX = pt.figure(axes=(5, 1))

for i, vparam in enumerate(vparam_values):

    for k, v in enumerate(viruses):
        pt.violin(CV_coef[f"{v}-{vparam}"], x=k*1, ax=AX[i], color=color_virus[v])
    
    pt.annotate(AX[i], '%s=%.0f ' % (varied_parameter, vparam), (0.5, 1.1), ha='center')
    pt.set_plot(AX[i], ylabel='CV' if i==0 else '')
    AX[i].set_xticks([0, 1], viruses, fontsize=7)

#%%
color_virus = {'sgRosa' : 'grey', 
               'sgCnr1': "darkred"}
xdispersion = 0.5

fig, AX = pt.figure(axes=(5, 1))

for i, vparam in enumerate(vparam_values):

    for k, v in enumerate(viruses):
        pt.bar([np.mean(fano_factor[f"{v}-{vparam}"])], x=[k], sy=np.std(fano_factor[f"{v}-{vparam}"]), ax=AX[i], color=color_virus[v])
        """ pt.scatter(np.ones(len(fano_factor[f"{v}-{vparam}"]))*k+np.random.rand(len(fano_factor[f"{v}-{vparam}"]))*xdispersion-xdispersion/2, 
                   fano_factor[f"{v}-{vparam}"], ax=AX[i], color='black', alpha=0.5) """

    pt.annotate(AX[i], '%s=%.0f ' % (varied_parameter, vparam), (0.5, 1.1), ha='center')
    pt.set_plot(AX[i], ylabel='Fano Factor' if i==0 else '')
    AX[i].set_xticks([0, 1], viruses, fontsize=7)

fig, AX = pt.figure(axes=(5, 1))

for i, vparam in enumerate(vparam_values):

    for k, v in enumerate(viruses):
        pt.violin(fano_factor[f"{v}-{vparam}"], x=k*1, ax=AX[i], color=color_virus[v])
    
    pt.annotate(AX[i], '%s=%.0f ' % (varied_parameter, vparam), (0.5, 1.1), ha='center')
    pt.set_plot(AX[i], ylabel='Fano Factor' if i==0 else '')
    AX[i].set_xticks([0, 1], viruses, fontsize=7)



# %%
from scipy.optimize import minimize
from scipy import stats

# fitting window
t_cond = ep.t > -0.1

# baseline window
baselineCond = (ep.t > -0.1) & (ep.t < 0)

gains = {}

for virus in ['sgRosa', 'sgCnr1']:

    gains[virus] = {}

    for img_id in [1.,2.,3.,4.,5.]:

        gains[virus][img_id] = []

        run_sessions   = means[f'{virus}-run-{img_id}']
        still_sessions = means[f'{virus}-still-{img_id}']

        # loop over sessions
        for run_data, still_data in zip(run_sessions, still_sessions):

            # shape = (episodes, neurons, time)

            nROIs = min(run_data.shape[1],
                        still_data.shape[1])

            for roi in range(nROIs):

                # average across episodes
                run_trace = np.mean(run_data[:, roi, :], axis=0)
                still_trace = np.mean(still_data[:, roi, :], axis=0)

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

                gains[virus][img_id].append(res.x[0])

# %%
fig, ax = plt.subplots(1,5, figsize=(12,3))

for i, img_id in enumerate([1.,2.,3.,4.,5.]):

    rosa = gains['sgRosa'][img_id]
    cnr1 = gains['sgCnr1'][img_id]

    ax[i].bar(
        [0,1],
        [np.mean(rosa), np.mean(cnr1)],
        yerr=[stats.sem(rosa), stats.sem(cnr1)],
        color=['grey','darkred']
    )

    ax[i].set_xticks([0,1])
    ax[i].set_xticklabels(['sgRosa','sgCnr1'])
    ax[i].set_ylabel('gain')
    ax[i].set_title(f'Image {img_id}')

plt.tight_layout()

# %%
fig, ax = plt.subplots(1,5, figsize=(14,3), sharey=True)

all_values = []

# first pass -> gather all gains for common y-limits
for virus in ['sgRosa', 'sgCnr1']:
    for img_id in [1.,2.,3.,4.,5.]:
        all_values.extend(gains[virus][img_id])

ymin = np.min(all_values)
ymax = np.max(all_values)

for i, img_id in enumerate([1.,2.,3.,4.,5.]):

    rosa = np.array(gains['sgRosa'][img_id])
    cnr1 = np.array(gains['sgCnr1'][img_id])

    # bars
    ax[i].bar(
        0,
        np.mean(rosa),
        yerr=stats.sem(rosa),
        color='grey',
        alpha=0.6,
        width=0.6
    )

    ax[i].bar(
        1,
        np.mean(cnr1),
        yerr=stats.sem(cnr1),
        color='darkred',
        alpha=0.6,
        width=0.6
    )

    # scatter points with jitter
    jitter1 = np.random.normal(0, 0.05, size=len(rosa))
    jitter2 = np.random.normal(0, 0.05, size=len(cnr1))

    ax[i].scatter(
        np.zeros(len(rosa)) + jitter1,
        rosa,
        color='black',
        s=12,
        alpha=0.7
    )

    ax[i].scatter(
        np.ones(len(cnr1)) + jitter2,
        cnr1,
        color='black',
        s=12,
        alpha=0.7
    )

    ax[i].set_xticks([0,1])
    ax[i].set_xticklabels(['sgRosa','sgCnr1'])

    ax[i].set_title(f'Image {img_id}')

    ax[i].set_ylim([ymin-0.1, ymax+0.1])

    if i == 0:
        ax[i].set_ylabel('gain')

plt.tight_layout()
# %%
plt.hist(gains['sgRosa'][1.0], bins=50)
plt.xlabel('gain')
plt.ylabel('count')