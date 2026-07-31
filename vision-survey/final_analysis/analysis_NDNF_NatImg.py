# %%
import numpy as np
import pandas as pd 
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
from scipy.stats import skew

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
pos_stat_test_props['interval_post'] = params.interval_post_ndnf[pname]
pos_stat_test_props['sign'] = 'positive'

neg_stat_test_props = params.stat_test_props.copy()
neg_stat_test_props['interval_post'] = params.interval_post_ndnf[pname]
neg_stat_test_props['sign'] = 'negative'

dFoF_options = params.dFoF_options
dFoF_options['with_computed_neuropil_fact'] = False

ep_props = dict(quantities=['dFoF', 'Deconvolved'] ,
                #prestim_duration=1.5,
                dt_sampling=params.dt_sampling,
                verbose=False)

savepath_fig = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\NDNF\figures\natural_images"
savepath_excel = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\NDNF\excels"
savepath_data = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\NDNF\data\natural_images"
os.makedirs(savepath_fig, exist_ok=True)
os.makedirs(savepath_data, exist_ok=True)
os.makedirs(savepath_excel, exist_ok=True)

# %% COMPUTATION

# Initialization

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
    
    print(i+1, '--', filename, '--', data.nROIs)

    if 'dFoF' in ep_props['quantities']:
        data.build_dFoF(**dFoF_options, verbose=True)

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
                                            repetition_keys=['repeat'],
                                            verbose=False
                                            )
        
        neg_evokedStats = ep.pre_post_statistics(\
                                            neg_stat_test_props,
                                            response_args=params.response_args,
                                            response_significance_threshold=params.response_significance_threshold,
                                            loop_over_cells=True,
                                            repetition_keys=['repeat'],
                                            verbose=False
                                            )

        ######### 3) Get episodes' traces #########
        
        img_id = 3.
        k=np.argwhere(ep.varied_parameters[varied_parameter]==img_id)[0][0]
        
        # INITILIAZE DICTIONARIES TO STORE RESPONSES AND BEHAVIORAL QUANTITIES
        if included_mice_pos is None:

            included_mice_pos  = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            included_mice_neg = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            dFoF_pos = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            dFoF_neg = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            deconvolved_pos = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            deconvolved_neg = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
                                    
        vparam_cond = (getattr(ep, varied_parameter)==img_id)

        for state, state_filter in zip(states_names, states_filters):

            if (pos_evokedStats['significant'][:, k].size != 0) and \
                (np.sum(pos_evokedStats['significant'][:, k], axis=0)>=params.NMIN_ROIS) and \
                    (np.sum(vparam_cond & state_filter) >= params.NMIN_EPISODES):
                
                dFoF_pos[f"{virus}-{state}"].append(
                    ep.dFoF[vparam_cond & state_filter][:, pos_evokedStats['significant'][:, k], :])
                
                
                deconvolved_pos[f"{virus}-{state}"].append(
                    ep.Deconvolved[vparam_cond & state_filter][:, pos_evokedStats['significant'][:, k], :])
                
                included_mice_pos[f"{virus}-{state}"].append(DATASET['subjects'][i])

            else:
                print("cond: %s-%s -> [XX] response not included (%i ROIs, %i eps)" % 
                    (virus,state, np.sum(pos_evokedStats['significant'][:, k]), np.sum(vparam_cond & state_filter)))
            
            if (neg_evokedStats['significant'][:, k].size != 0) and \
                (np.sum(neg_evokedStats['significant'][:, k], axis=0)>=params.NMIN_ROIS) and \
                    (np.sum(vparam_cond & state_filter) >= params.NMIN_EPISODES):
                
                dFoF_neg[f"{virus}-{state}"].append(
                    ep.dFoF[vparam_cond & state_filter][:, neg_evokedStats['significant'][:, k], :])
                
                deconvolved_neg[f"{virus}-{state}"].append(
                    ep.Deconvolved[vparam_cond & state_filter][:, neg_evokedStats['significant'][:, k], :])
                
                included_mice_neg[f"{virus}-{state}"].append(DATASET['subjects'][i])

            else:
                print("cond: %s-%s -> [XX] response not included (%i ROIs, %i eps)" % 
                    (virus,state, np.sum(neg_evokedStats['significant'][:, k]), np.sum(vparam_cond & state_filter)))

        ###############################

    else :
        print("session %s has no ROIs, excluded from analysis" % filename)
    
    print('')

dFoF_pos = tools.remove_empty_sessions(dFoF_pos)
deconvolved_pos = tools.remove_empty_sessions(deconvolved_pos)
dFoF_neg = tools.remove_empty_sessions(dFoF_neg)
deconvolved_neg = tools.remove_empty_sessions(deconvolved_neg)

#%%
dFoF_pos_means = {f"{virus}-{state}": [] for virus, state in product(viruses, states_names)}
dFoF_neg_means = {f"{virus}-{state}": [] for virus, state in product(viruses, states_names)}
for key in dFoF_pos.keys():
    for i in range(len(dFoF_pos[key])):
        dFoF_pos_means[key].append(np.mean(dFoF_pos[key][i], axis=1))
    for i in range(len(dFoF_neg[key])):
        dFoF_neg_means[key].append(np.mean(dFoF_neg[key][i], axis=1))

dFoF_pos_means_trials = {}
l = []
for key in dFoF_pos_means.keys():
    print("pos", key, len(dFoF_pos_means[key]))
    for i in range(len(dFoF_pos_means[key])):
        l.append(dFoF_pos_means[key][i])
    dFoF_pos_means_trials[key] = np.array(l)

dFoF_neg_means_trials = {}
l = []
for key in dFoF_neg_means.keys():
    print("neg", key, len(dFoF_neg_means[key]))
    for i in range(len(dFoF_neg_means[key])):
        l.append(dFoF_neg_means[key][i])
    dFoF_neg_means_trials[key] = np.array(l)

#%%
from scipy.stats import sem

color_virus = {'sgRosa-all' : 'grey', 
               'sgCnr1-all': "darkred"}

fig, AX = plt.subplots(2, 10, figsize=(15, 3), sharex=True, sharey="row")

for j in range(10):
    for k, key in enumerate(dFoF_pos_means_trials.keys()):
        pt.plot(ep.t, np.mean(dFoF_pos_means_trials[key], axis=0)[j] - np.mean(dFoF_pos_means_trials[key], axis=0)[j][(ep.t<0) & (ep.t>-1)].mean(), 
                #sy=sem(dFoF_pos_means_trials[key], axis=0)[j],
                ax=AX[0][j], color=color_virus[key])
        if j == 9:
            nb_mice = len(np.unique(included_mice_pos[key]))
            pt.annotate(AX[0][9], 'N=%i (%i mice)' % (len(included_mice_pos[key]), nb_mice)+k*'\n',
                        color=color_virus[key], **dict(xy=(1.1,0.7), ha='left', fontsize=7))
    
    for k, key in enumerate(dFoF_pos_means_trials.keys()):
        pt.plot(ep.t, np.mean(dFoF_neg_means_trials[key], axis=0)[j] - np.mean(dFoF_neg_means_trials[key], axis=0)[j][(ep.t<0) & (ep.t>-1)].mean(), 
                #sy=sem(dFoF_neg_means_trials[key], axis=0)[j],
                ax=AX[1][j], color=color_virus[key])
        if j == 9:
            nb_mice = len(np.unique(included_mice_neg[key]))
            pt.annotate(AX[1][9], 'N=%i (%i mice)' % (len(included_mice_neg[key]), nb_mice)+k*'\n',
                        color=color_virus[key], **dict(xy=(1.1,0.7), ha='left', fontsize=7))

    AX[0][j].set_title("ep {}".format(j+1))

plt.show()

firgurename = 'natimg_response_over_ep_' + 'NDNF' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")

#%%
import random

fig, AX = plt.subplots(5, 10, figsize=(15, 6), sharex=True, sharey="row")

for j in range(10):
    for k, key in enumerate(dFoF_pos_means_trials.keys()):
        pt.plot(ep.t, np.mean(dFoF_pos_means_trials[key], axis=0)[j] - np.mean(dFoF_pos_means_trials[key], axis=0)[j][(ep.t<0) & (ep.t>-1)].mean(), 
                        sy=sem(dFoF_pos_means_trials[key], axis=0)[j],
                        ax=AX[0][j], color=color_virus[key])
        if j == 9:
            nb_mice = len(np.unique(included_mice_pos[key]))
            pt.annotate(AX[0][9], 'N=%i (%i mice)' % (len(included_mice_pos[key]), nb_mice)+k*'\n',
                        color=color_virus[key], **dict(xy=(1.1,0.7), ha='left', fontsize=7))
    AX[0][j].set_title("ep {}".format(j+1))
    AX[0][0].set_ylabel("Average")


for k, key in enumerate(dFoF_pos_means_trials.keys()):
    num_session = random.sample(range(0, len(dFoF_pos[key])-1), 3)
    for i in range(2):
        num_roi = random.randint(0, dFoF_pos[key][num_session[i]].shape[1]-1)
        for j in range(10):
            pt.plot(ep.t, dFoF_pos[key][num_session[i]][:, num_roi, :][j] - dFoF_pos[key][num_session[i]][:, num_roi, :][j][(ep.t<0) & (ep.t>-1)].mean(), 
            #pt.plot(ep.t, dFoF_pos[key][num_session[i]][:, num_roi, :][j], 
                    ax=AX[1+k*2+i][j], color=color_virus[key])
        AX[1+k*2+i][0].set_ylabel(f"Sess. {num_session[i]+1}, ROI #{num_roi+1}")

firgurename = 'natimg_pos_response_over_ep_' + 'NDNF' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")

#%%
import random

fig, AX = plt.subplots(5, 10, figsize=(15, 6), sharex=True, sharey="row")

for j in range(10):
    for k, key in enumerate(dFoF_neg_means_trials.keys()):
        pt.plot(ep.t, np.mean(dFoF_neg_means_trials[key], axis=0)[j] - np.mean(dFoF_neg_means_trials[key], axis=0)[j][(ep.t<0) & (ep.t>-1)].mean(), 
                        sy=sem(dFoF_neg_means_trials[key], axis=0)[j],
                        ax=AX[0][j], color=color_virus[key])
        if j == 9:
            nb_mice = len(np.unique(included_mice_neg[key]))
            pt.annotate(AX[0][9], 'N=%i (%i mice)' % (len(included_mice_neg[key]), nb_mice)+k*'\n',
                        color=color_virus[key], **dict(xy=(1.1,0.7), ha='left', fontsize=7))
    AX[0][j].set_title("ep {}".format(j+1))
    AX[0][0].set_ylabel("Average")


for k, key in enumerate(dFoF_neg_means_trials.keys()):
    num_session = random.sample(range(0, len(dFoF_neg[key])-1), 3)
    for i in range(2):
        num_roi = random.randint(0, dFoF_neg[key][num_session[i]].shape[1]-1)
        for j in range(10):
            pt.plot(ep.t, dFoF_neg[key][num_session[i]][:, num_roi, :][j] - dFoF_neg[key][num_session[i]][:, num_roi, :][j][(ep.t<0) & (ep.t>-1)].mean(), 
            #pt.plot(ep.t, dFoF_neg[key][num_session[i]][:, num_roi, :][j], 
                    ax=AX[1+k*2+i][j], color=color_virus[key])
        AX[1+k*2+i][0].set_ylabel(f"Sess. {num_session[i]+1}, ROI #{num_roi+1}")

firgurename = 'natimg_neg_response_over_ep_' + 'NDNF' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")