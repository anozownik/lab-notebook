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
                        'PN_cond-NDNF-CB1_WT-vs-KD', 'final_NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol=pname)

#%% SET PARAMETERS

state_metric = 'speed' # 'speed' or 'pupil' or 'speed & pupil'

stat_test_props = params.stat_test_props.copy()
stat_test_props['interval_post'] = params.interval_post_pyr[pname]
stat_test_props['sign'] = 'positive'

dFoF_options = params.dFoF_options
dFoF_options['with_computed_neuropil_fact'] = True

ep_props = dict(quantities=['dFoF', 'Deconvolved', 'running_speed'],
                #prestim_duration=1.5,
                dt_sampling=params.dt_sampling,
                verbose=False)

# decoder parameters
nSeed = 20
denoising = True
decoder_window = [[2, 3], [0.5, 1.5]]

savepath_fig = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\PN\figures\natural_images"
savepath_excel = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\PN\excels"
savepath_data = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\PN\data\natural_images"
os.makedirs(savepath_fig, exist_ok=True)
os.makedirs(savepath_data, exist_ok=True)
os.makedirs(savepath_excel, exist_ok=True)

# %% COMPUTATION

# Initialization
DATASET["virus"], DATASET['alpha'] = [], []
DATASET["nb_neurons"], DATASET["nb_responsive_neurons"], DATASET['mean_dFoF'] = [], [], []
for i in range(len(decoder_window)):
    key = f"accuracies_{decoder_window[i][0]}-{decoder_window[i][1]}"
    DATASET[key] = []

DATASET_ROIS = dict(mouse=[], virus=[], session=[], rois_id=[], 
                    mean_avg_trial=[], mean_avg_trial_win=[],
                    responsive=[], skewness=[])

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

        ######### 1) Define behavioral states #########
        states_names, states_filters = tools.define_trials_arousal_state(ep, cond=state_metric)

        ######### 2) identify visually-responsive cells #########

        evokedStats = ep.pre_post_statistics(\
                                            stat_test_props,
                                            response_args=params.response_args,
                                            response_significance_threshold=params.response_significance_threshold,
                                            loop_over_cells=True,
                                            repetition_keys=[varied_parameter, 'repeat'],
                                            verbose=False
                                            )

        # find responsive ROIs for this contrast (from summary stats)
        DATASET["nb_responsive_neurons"].append(np.sum(evokedStats['significant'], axis=0) if evokedStats['significant'].size != 0 else 0)
        print(data.nROIs, 'ROIs in this session,', np.sum(evokedStats['significant'], axis=0), 'responsive ROIs')

        window = (ep.t>stat_test_props['interval_post'][0]) & (ep.t<stat_test_props['interval_post'][1])
        if evokedStats['significant'].size != 0:
            DATASET['mean_dFoF'].append(np.mean(ep.dFoF[:, evokedStats['significant'], :][:, :, window]))
        else :
            DATASET['mean_dFoF'].append(np.nan)

        ######### 3) Get episodes' traces #########

        # INITILIAZE DICTIONARIES TO STORE RESPONSES AND BEHAVIORAL QUANTITIES
        if included_mice is None:

            included_mice  = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            dFoF = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            deconvolved = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}
            run_means = {f"{v}-{cond}" : [] for v, cond in product(viruses, states_names)}

        for state, state_filter in zip(states_names, states_filters):

            if (evokedStats['significant'].size != 0) and \
                (np.sum(evokedStats['significant'], axis=0)>=params.NMIN_ROIS) and \
                    (np.sum(state_filter) >= params.NMIN_EPISODES):
                
                print("cond: %s-%s -> included %i ROIs and %i episodes" % 
                      (virus, state, np.sum(evokedStats['significant']), np.sum(state_filter)))
                print("    -> %i ROIs out of %i ROIs are responsive" % 
                      (np.sum(evokedStats['significant'], axis=0), evokedStats['significant'].shape[0]))
                
                dFoF[f"{virus}-{state}"].append(
                    ep.dFoF[state_filter][:, evokedStats['significant'], :])
                
                deconvolved[f"{virus}-{state}"].append(
                    ep.Deconvolved[state_filter][:, evokedStats['significant'], :])
                
                included_mice[f"{virus}-{state}"].append(DATASET['subjects'][i])

                run_means[f"{virus}-{state}"].append(ep.running_speed[state_filter, :])
                
            else:
                print("cond: %s-%s -> [XX] response not included (%i ROIs, %i eps)" % 
                      (virus,state, np.sum(evokedStats['significant']), np.sum(state_filter)))
            
        ######### 4) Decoder #########

        from sklearn.neighbors import KNeighborsClassifier
        from sklearn.model_selection import train_test_split
        from sklearn.metrics import accuracy_score

        def build_X_y(ep, 
            averaging_window=[2.5, 3.5], evokedStats=None):

            if evokedStats is None:
                dFoF = ep.dFoF
            else :
                dFoF = ep.dFoF[:, evokedStats['significant'], :]

            averaging_window_cond = (ep.t> averaging_window[0])& (ep.t<averaging_window[1])
            X = np.zeros((dFoF.shape[0], dFoF.shape[1])) # list of matrice responses (Nrois, Ntrials) 

            y = np.zeros(dFoF.shape[0], dtype=int) # label of all trials 

            i=0
            for id in np.unique(getattr(ep, 'Image-ID')):
                pattern_cond = (getattr(ep, 'Image-ID')==id)
                X[i:i+np.sum(pattern_cond), :] =\
                    dFoF[pattern_cond,:,:][:,:,averaging_window_cond].mean(axis=2)
                y[i:i+np.sum(pattern_cond)] = id
                i+=np.sum(pattern_cond)

            return X, y

        def average_per_class(X, y):
            X_avg = []
            y_avg = []
            for label in np.unique(y):
                X_avg.append(X[y == label].mean(axis=0))
                y_avg.append(label)
            return np.array(X_avg), np.array(y_avg)

        for i in range(len(decoder_window)):
            key = f"accuracies_{decoder_window[i][0]}-{decoder_window[i][1]}"

            model = KNeighborsClassifier(n_neighbors=1)
            X, y = build_X_y(ep, decoder_window[i])
            #X, y = build_X_y(ep, decoder_window[i], evokedStats)

            accuracies = []
            for random_state in range(nSeed):
                X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                            test_size=0.5,
                                                            random_state=random_state,
                                                            stratify=y)
                if denoising:
                    X_train, y_train = average_per_class(X_train, y_train)
                
                model.fit(X_train, y_train)
                y_pred = model.predict(X_test)
                score=accuracy_score(y_test, y_pred)
                accuracies.append(score)
            
            DATASET[key].append(np.mean(accuracies))

        ######### 5) Store ROIs' data #########
        for roi_idx in range(len(evokedStats['significant'])):
            DATASET_ROIS['mouse'].append(DATASET['subjects'][i])
            DATASET_ROIS['virus'].append(virus)
            DATASET_ROIS['session'].append(os.path.basename(filename)[:-4])
            DATASET_ROIS['rois_id'].append(roi_idx)
            DATASET_ROIS['skewness'].append(skew(np.mean(ep.dFoF[:, roi_idx, :], axis=0)))
            roi_bsl = np.mean(ep.dFoF[:, roi_idx, :][:, (ep.t>-1) & (ep.t<0)])
            if evokedStats['significant'][roi_idx]:
                DATASET_ROIS['responsive'].append(1)
                DATASET_ROIS['mean_avg_trial'].append(np.mean(ep.dFoF[:, roi_idx, :][:, (ep.t>0)]) - roi_bsl)
                DATASET_ROIS['mean_avg_trial_win'].append(np.mean(ep.dFoF[:, roi_idx, :][:, window]) - roi_bsl)
            else:
                DATASET_ROIS['responsive'].append(0)
                DATASET_ROIS['mean_avg_trial'].append(np.mean(ep.dFoF[:, roi_idx, :][:, (ep.t>0)]) - roi_bsl)
                DATASET_ROIS['mean_avg_trial_win'].append(np.mean(ep.dFoF[:, roi_idx, :][:, window]) - roi_bsl)

        ###############################

    else :
        print("session %s has no ROIs, excluded from analysis" % filename)
        DATASET["nb_responsive_neurons"].append(0)
        DATASET['mean_dFoF'].append(np.nan)
    
    print('')
#%%

# BUILD DATAFRAME ROIS
df_rois = pd.DataFrame(DATASET_ROIS)

excel_filename = pname + '_summary_data_rois_' + 'PN' + '.xlsx'
df_rois.to_excel(os.path.join(savepath_excel, excel_filename))

df_rois

#%%           
dFoF = tools.remove_empty_sessions(dFoF)
deconvolved = tools.remove_empty_sessions(deconvolved)
run_means = tools.remove_empty_sessions(run_means)

# SAVE IN NPY FILES
np.save(os.path.join(savepath_data, 'dFoF_means.npy'), dFoF)
np.save(os.path.join(savepath_data, 'deconvolved_means.npy'), deconvolved)
np.save(os.path.join(savepath_data, 'run_means.npy'), run_means)
np.save(os.path.join(savepath_data, 'included_mice.npy'), included_mice)

# BUILD DATAFRAME

DATASET['session'] = np.array([os.path.basename(file)[:-4] for file in DATASET['files']])

df = pd.DataFrame(DATASET).drop(columns=['protocol_ids', 'protocols', 'files']).set_index('session')
df['perc_resp'] = df['nb_responsive_neurons'] / df['nb_neurons'] * 100

excel_filename = 'natimg_summary_data_' + 'PYR' + '.xlsx'
df.to_excel(os.path.join(savepath_excel, excel_filename))

df

#%% 1.a Averaged dF/F0 over responsive ROIs and over episodes across virus and behavioral states (std over sessions)

fig, AX = pt_fcts.plot_average_response(ep.t, dFoF, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS)

firgurename = 'natimg_average_res_dfof_'+ 'PYR' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")

#%% 1.b Averaged dF/F0 with baseline substracted over responsive ROIs and over episodes across virus and behavioral states (std over sessions)
baselineCond = (ep.t>stat_test_props['interval_pre'][0]) & (ep.t<stat_test_props['interval_pre'][1])

fig, AX = pt_fcts.plot_average_response(ep.t, dFoF, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS, 
                                        baselineSubtraction=True, baselineCond=baselineCond)

firgurename = 'natimg_average_res_dfof_bsl_sub_'+ 'PYR' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")

#%% 1.c Averaged deconvolved trace over responsive ROIs and over episodes across virus and behavioral states (std over sessions)

fig, AX = pt_fcts.plot_average_response(ep.t, deconvolved, 
                                        viruses, states_names, [], varied_parameter, 
                                        included_mice, params.NMIN_SESSIONS)
AX[0].set_ylabel('deconvolved')

firgurename = 'natimg_average_res_deconvolved_' + 'PYR' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")

#%% 2. Pie chart of responsive neurons
percentages = {}

for v in viruses:
    percentages[v] = (df[df['virus'] == v]['perc_resp']).to_list()

fig, AX = pt_fcts.pie_chart_responsive_neurons(percentages, viruses, [])

firgurename = 'natimg_pie_charts_' + 'PYR' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")

#%% 3. Speed across virus and behavioral states

fig, AX = pt_fcts.plot_average_behavior(ep.t, run_means, viruses, states_names, params.NMIN_SESSIONS, ylabel='speed(cm/s)')

firgurename = 'natimg_average_speed_' + 'PYR' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")

# %% plot gain of arousal modulation by fitting a multiplicative factor to the still condition to explain the run condition
from scipy.optimize import minimize
from scipy import stats
from itertools import compress

# fitting window
t_cond = ep.t > 0.

# baseline window
baselineCond = (ep.t>stat_test_props['interval_pre'][0]) & (ep.t<stat_test_props['interval_pre'][1])

gains = {}
mice = {}

for virus in ['sgRosa', 'sgCnr1']:

    gains[virus] = []
    mice[virus] = []

    run_sessions   = dFoF[f'{virus}-run']
    still_sessions = dFoF[f'{virus}-still']
    run_mice = np.array(included_mice[f'{virus}-run']) 
    still_mice = np.array(included_mice[f'{virus}-still'])

    for mouse in np.unique(run_mice):

        run_sessions_mouse = list(compress(dFoF[f'{virus}-run'], run_mice==mouse))
        still_sessions_mouse = list(compress(dFoF[f'{virus}-still'], still_mice==mouse))

        run_sessions_mouse_average = []
        still_sessions_mouse_average = []

        # loop over sessions
        for run_data, still_data in zip(run_sessions_mouse, still_sessions_mouse):

            # average across episodes
            run_trace = np.mean(run_data, axis=(0, 1))
            still_trace = np.mean(still_data, axis=(0, 1))

            run_sessions_mouse_average.append(run_trace)
            still_sessions_mouse_average.append(still_trace)
        
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
        mice[virus].append(mouse)

df2 = pd.DataFrame({'mouse': mice['sgRosa']+mice['sgCnr1'], 
              'virus': ['sgRosa']*len(mice['sgRosa'])+['sgCnr1']*len(mice['sgCnr1']), 
              'gain': gains['sgRosa']+gains['sgCnr1']})

with pd.ExcelWriter(os.path.join(savepath_excel, excel_filename)) as writer:
    df.to_excel(writer, sheet_name='sessions')
    df2.to_excel(writer, sheet_name='mice')

# %% 

rosa = gains['sgRosa']
cnr1 = gains['sgCnr1']

fig, AX = pt.figure(axes=(1, 1))

AX.bar([0,1],[np.mean(rosa), np.mean(cnr1)],yerr=[stats.sem(rosa), stats.sem(cnr1)], color=['grey','darkred'])
AX.scatter(np.ones(len(rosa))*0+np.random.rand(len(rosa))*0.2-0.1, rosa, color='black', alpha=0.5)
AX.scatter(np.ones(len(cnr1))*1+np.random.rand(len(cnr1))*0.2-0.1, cnr1, color='black', alpha=0.5)

AX.set_xticks([0,1])
AX.set_xticklabels(['sgRosa','sgCnr1'])
AX.set_ylabel('gain')
pt.annotate(AX, f'WT: N={len(rosa)} mice\nKD: N={len(cnr1)} mice', xy=(0.5, 1), xycoords='axes fraction', ha='right', fontsize=4)

firgurename = 'natimg_gain_' + 'PYR' + '.svg'
fig.savefig(os.path.join(savepath_fig, firgurename), transparent=True, format='svg', bbox_inches="tight")

# %% COMPUTATION Highly skewd vs Lowly skewed Average response

ep_props = dict(quantities=['dFoF'],
                #prestim_duration=1.5,
                dt_sampling=params.dt_sampling,
                verbose=False)

viruses = ['sgRosa', 'sgCnr1']

# INITILIAZE DICTIONARIES TO STORE RESPONSES AND BEHAVIORAL QUANTITIES
included_mice_highly_skewed  = {f"{v}" : [] for v in viruses}
included_mice_lowly_skewed  = {f"{v}" : [] for v in viruses}

trial_avg_res_highly_skewed = {f"{v}" : [] for v in viruses}
trial_avg_res_lowly_skewed = {f"{v}" : [] for v in viruses}

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

        # Build episode data
        ep = physion.analysis.episodes.build.EpisodeData(data, 
                                                        **ep_props,
                                                        protocol_name=pname)
        
        ######### 1) identify highly skewed neurons #########
        highly_skewed_rois = skew(data.dFoF, axis=1) >= 2.7

        ######### 2) identify visually-responsive cells #########

        evokedStats = ep.pre_post_statistics(\
                                            stat_test_props,
                                            response_args=params.response_args,
                                            response_significance_threshold=params.response_significance_threshold,
                                            loop_over_cells=True,
                                            repetition_keys=[varied_parameter, 'repeat'],
                                            verbose=False
                                            )

        ######### 3) Get episodes' traces #########

        if (evokedStats['significant'].size != 0) and \
            (np.sum(evokedStats['significant'], axis=0)>=params.NMIN_ROIS):

            if np.sum(highly_skewed_rois & evokedStats['significant']) > 0:
            
                trial_avg_res_highly_skewed_session = np.mean(ep.dFoF[:, highly_skewed_rois & evokedStats['significant'], :], axis=0)
                trial_avg_res_highly_skewed[virus].append(trial_avg_res_highly_skewed_session)
                included_mice_highly_skewed[f"{virus}"].append(DATASET['subjects'][i])
            
            if np.sum((~highly_skewed_rois) & evokedStats['significant']) > 0:
                trial_avg_res_lowly_skewed_session = np.mean(ep.dFoF[:, (~highly_skewed_rois) & (evokedStats['significant']), :], axis=0)
                trial_avg_res_lowly_skewed[virus].append(trial_avg_res_lowly_skewed_session)
                included_mice_lowly_skewed[f"{virus}"].append(DATASET['subjects'][i])
                
        else:
            print("cond: %s -> [XX] response not included (%i ROIs)" % 
                    (virus, np.sum(evokedStats['significant'])))

        ###############################

    else :
        print("session %s has no ROIs, excluded from analysis" % filename)
    
    print('')

trial_avg_res_highly_skewed['sgRosa-all'] = trial_avg_res_highly_skewed['sgRosa']
trial_avg_res_highly_skewed['sgCnr1-all'] = trial_avg_res_highly_skewed['sgCnr1']
del trial_avg_res_highly_skewed['sgRosa']
del trial_avg_res_highly_skewed['sgCnr1']

trial_avg_res_lowly_skewed['sgRosa-all'] = trial_avg_res_lowly_skewed['sgRosa']
trial_avg_res_lowly_skewed['sgCnr1-all'] = trial_avg_res_lowly_skewed['sgCnr1']
del trial_avg_res_lowly_skewed['sgRosa']
del trial_avg_res_lowly_skewed['sgCnr1']

included_mice_highly_skewed['sgRosa-all'] = included_mice_highly_skewed['sgRosa']
included_mice_highly_skewed['sgCnr1-all'] = included_mice_highly_skewed['sgCnr1']
del included_mice_highly_skewed['sgRosa']
del included_mice_highly_skewed['sgCnr1']

included_mice_lowly_skewed['sgRosa-all'] = included_mice_lowly_skewed['sgRosa']
included_mice_lowly_skewed['sgCnr1-all'] = included_mice_lowly_skewed['sgCnr1']
del included_mice_lowly_skewed['sgRosa']
del included_mice_lowly_skewed['sgCnr1']

#%%

color_virus = {'sgRosa' : 'grey', 
               'sgCnr1': "darkred"}

def plot_average_response(t, responses, 
                          viruses, states, varied_parameter=[],
                          vparam_name='',
                          included_mice=None,
                          nmin_sessions=1,
                          baselineSubtraction=False,
                          baselineCond=None,
                          annotation_props=dict(xy=(0.05,1), ha='left', fontsize=4),
                          savepath=None):
    
    if len(varied_parameter) == 0:
        fig, AX = plot_average_response_no_vparam(t, responses, 
                                                  viruses, states, 
                                                  included_mice, nmin_sessions, 
                                                  baselineSubtraction, baselineCond,
                                                  annotation_props)

    if savepath is not None:
        plt.savefig(os.path.join(savepath), transparent=True, format='svg')

    return fig, AX

def plot_average_response_no_vparam(t, responses, 
                                    viruses, states, 
                                    included_mice=None,
                                    nmin_sessions=1,
                                    baselineSubtraction=False,
                                    baselineCond=None,
                                    annotation_props=dict(xy=(0.05,1), ha='left', fontsize=4)):
    
    fig, AX = pt.figure(axes=(len(states), 1))

    if len(states) == 1:
        AX = np.reshape(AX, (1))

    for j, state in enumerate(states):
        for k, virus in enumerate(viruses):

            key = f'{virus}-{state}'

            session_responses = np.array([np.mean(r, axis=(0)) for r in responses[key]])
            
            if baselineSubtraction:
                if baselineCond is None:
                    baselineCond = (t<0)
                if len(session_responses.shape) == 2:
                    session_responses = session_responses - session_responses[:, baselineCond].mean(axis=1, keepdims=True)

            if np.shape(session_responses)[0] >= nmin_sessions :

                if np.shape(session_responses)[0] == 1 :
                    pt.plot(t, session_responses[0],
                            color=color_virus[virus], ax=AX[j])
                else :
                    pt.plot(t, np.mean(session_responses,axis=0),
                            sy=sem(session_responses,axis=0),
                            color=color_virus[virus], ax=AX[j])
            
            if len(responses[key]) > 0:

                nb_rois = np.column_stack([(r.shape[0]) for r in responses[key]])

                if included_mice is not None:
                    nb_mice = np.unique(included_mice[key]).shape[0]
                    pt.annotate(AX[j],
                                'N=%i (%i mice, %i rois)' % (len(responses[key]), 
                                                                        nb_mice,
                                                                        nb_rois.sum())
                                                                +k*'\n',
                                                                
                                color=color_virus[virus], **annotation_props)
                else :
                    pt.annotate(AX[j],
                                'N=%i (%i rois)' % (len(responses[key]), 
                                                            nb_rois.sum())
                                                            +k*'\n',
                                                                
                                color=color_virus[virus], **annotation_props)          
            else :
                pt.annotate(AX[j], 'N=0' +k*'\n', color=color_virus[virus], **annotation_props)
    
            pt.annotate(AX[j], state, (0.5, 1.3), ha='center') 
        
        pt.set_plot(AX[j], xlabel='time (s)', ylabel='$\\Delta$F/F' if j==0 else '')
    
    #pt.set_common_ylims(AX)

    return fig, AX

baselineCond = (ep.t>stat_test_props['interval_pre'][0]) & (ep.t<stat_test_props['interval_pre'][1])

fig, AX = plot_average_response(ep.t, trial_avg_res_highly_skewed, 
                                        viruses, ['all'], [], varied_parameter, 
                                        included_mice_highly_skewed, params.NMIN_SESSIONS, 
                                        baselineSubtraction=True, baselineCond=baselineCond)
AX[0].set_title('highly skewed', y=1.5)
AX[0].set_yticks([0., 0.2, 0.4, 0.6, 0.8])

                        
fig, AX = plot_average_response(ep.t, trial_avg_res_lowly_skewed, 
                                        viruses, ['all'], [], varied_parameter, 
                                        included_mice_lowly_skewed, params.NMIN_SESSIONS, 
                                        baselineSubtraction=True, baselineCond=baselineCond)
AX[0].set_title('lowly skewed', y=1.5)
AX[0].set_yticks([0., 0.2, 0.4, 0.6, 0.8])
