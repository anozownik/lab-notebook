import numpy as np
import os, sys
sys.path.append('../../physion/src') # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')
from scipy import stats
import matplotlib.pyplot as plt

    
# %%

# custom dFoF:

from physion.imaging.Calcium import compute_F0


dFoF_options = dict(\
    method_for_F0='sliding_minmax',
    #method_for_F0='sliding_percentile',
    #method_for_F0='percentile', #more strict
    sliding_window= 300,
    percentile=10,
    roi_to_neuropil_fluo_inclusion_factor=1.1,
    neuropil_correction_factor=0.7, 
    with_computed_neuropil_fact=True)




# TO LOOP OVER NWB FILES WITH VISUAL STIMULUS --- DRIFITING GRATING ---  multisession

folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna',
                        'PN_cond-NDNF-CB1_WT-vs-KD', '20260325', 'PNs', 'NWBs', '2026_march_and_april')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol='Natural-Images-4-repeats')



# STATISTICS PROPERTIES --- DRIFTING GRATINGS ---
stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.8],                                   
                       test='wilcoxon',
                       sign='both')

response_significance_threshold = 0.05

# PLOT PROPERTIES --- DRIFTING GRATINGS ---

plot_props = dict(column_key='Image-ID',
                  with_annotation=True,
                  Ybar=0.5, Ybar_label="0.5$\Delta$F/F",
                  Xbar=0.5, Xbar_label="0.5s",
                  figsize=(13,1.8))


# RESPONSE ARGUMENTS --- DRIFTING GRATINGS ---



response_args = dict(quantity='dFoF')

summary_stats = []

RUNNING_SPEED_THRESHOLD = 0.5
NMIN_ROIS = 1
NMIN_EPISODES = 1

means = {} # 
for virus in ['sgRosa', 'sgCnr1']:
      for cond in ['all', 'run', 'still']:
              for c, img_id in enumerate([1., 2., 3., 4., 5.]):
                means['%s-%s-%s' % (virus, cond, img_id)] = []
percentages = {}
for virus in ['sgRosa', 'sgCnr1']:
      
        for c, img_id in enumerate([1., 2., 3., 4., 5.]):
                percentages['%s-%s' % (virus, img_id)] = []

for i, filename in enumerate(DATASET['files']):
    
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    
    print(i+1, '--', filename, '--', data.nROIs)
    # print(data.protocols)

    data.build_dFoF(**dFoF_options, verbose=True)
    #data.build_pupil_diameter()
    #data.build_facemotion()
    data.build_running_speed()
    
    if data.nROIs>0:

        ep = physion.analysis.episodes.build.EpisodeData(data, 
                                                        quantities=['dFoF', 'running_speed'],
                                                        prestim_duration=1.5,
                                                        protocol_name='Natural-Images-4-repeats')
        
        if 'Image-ID' in ep.varied_parameters:
               
                # determine virus        
                if 'sgRosa' in data.nwbfile.virus:
                        virus = 'sgRosa'
                elif 'sgCnr1':
                        virus = 'sgCnr1'

                # 1) identify visually-responsive cells

                
                evokedStats = ep.pre_post_statistics(\
                                        stat_test_props,
                                        response_args=\
                                        dict(quantity='dFoF'),
                                        response_significance_threshold=response_significance_threshold,
                                        loop_over_cells=True,
                                        repetition_keys=['repeat']
                                        )
                
                # 2) split rest / run
                withinEpisode = (ep.t>0) & (ep.t<ep.time_duration[0])
                run = np.mean(ep.running_speed[:,withinEpisode], axis=1) > RUNNING_SPEED_THRESHOLD

                for img_id in ep.varied_parameters['Image-ID']:
                
                        # find responsive ROIs for this ImageID (from summary stats)
                        imageCond = evokedStats['Image-ID']==img_id
                        responsiveROIs = evokedStats['significant'][:,imageCond].flatten()
                        percentages['%s-%s' % (virus, img_id)].append(np.sum(responsiveROIs)/len(responsiveROIs)*100)
                        
                        
                               
                        image_cond = (getattr(ep, 'Image-ID')==img_id)

                #
                        #print("for session: %s" % filename)
                        for cond, filter in zip(['all', 'run', 'still'],
                                                [run|~run, run, ~run]):
                                
                                if (np.sum(responsiveROIs)>=NMIN_ROIS) and \
                                        (np.sum(image_cond & filter)>= NMIN_EPISODES):
                                        print("cond: %s-%s-%s -> included %i ROIs and %i episodes" % (virus,cond,img_id, np.sum(responsiveROIs), np.sum(image_cond & filter)))
                                        #print("cond: %s-%s-%s -> %i ROIs out of %i ROIs are responsive" % (cond,virus,img_id, np.sum(responsiveROIs), len(responsiveROIs)))
                                        
                                        means['%s-%s-%s' % (virus, cond, img_id)].append(
                                                ep.dFoF[image_cond & filter, :, :][:, responsiveROIs, :])
                                        
                                        
                                        
                                
                                        
                                else:
                                        print("cond: %s -> [XX] response not included (%i ROIs, %i eps)" % (cond, np.sum(responsiveROIs), np.sum(image_cond & filter)))
                

# now "means" is a list (over sessions) 
#    of responses of shape (episodes, responsiveROIs, time)

#%%
# 
if 'PN_cond-NDNF-CB1_WT-vs-KD' in folder:
       neuron= 'PN-cond-NDNF-CB1'
elif 'NDNF-cond-CB1_WT-vs-KD':
       neuron= 'NDNF-cond-CB1'


figurepath = '/Users/macbookair/work/Figures/'
firgurename = 'dfof_natimgID_beh_mod_'+ neuron + '.svg'


from scipy.stats import sem

fig, AX = pt.figure(axes=(3,5))

NMIN_SESSIONS = 1

for j, cond in enumerate(['all', 'run', 'still']):
    for i, img_id in enumerate([1., 2., 3., 4., 5.]):
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):

                session_responses = [np.mean(m,axis=(0,1))\
                        for m in means['%s-%s-%s' % (virus,cond,img_id)]]
                
                if len(session_responses)>=NMIN_SESSIONS:
                        pt.plot(ep.t, 
                                np.mean(session_responses,axis=0),
                                sy=sem(session_responses, axis=0),
                                color=color, ax=AX[i][j])
                        
                pt.annotate(AX[i][j],
                            'N=%i' % len(session_responses)+k*'\n',
                            (0.2,0.6), ha='right',
                            color=color, fontsize=6)
        if i==0:
             pt.annotate(AX[i][j], cond, (0.5, 1))
        if j==0:
             pt.annotate(AX[i][j], 'Image-ID= %s ' % img_id,
                         (0,1), ha='right')

        pt.set_plot(AX[i][j], 
                xlabel='time (s)' if i==2 else '',
                ylabel='$\\Delta$F/F' if j==0 else '')
pt.set_common_ylims(AX)  

#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')
#%%
baselineCond = (ep.t>-0.1) & (ep.t<0)

fig, AX = pt.figure(axes=(3,5))

NMIN_SESSIONS = 0

for j, cond in enumerate(['all', 'run', 'still']):
    for i, img_id in enumerate([1., 2., 3., 4., 5.]):
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):

                session_responses = [np.mean(m,axis=(0,1))\
                        for m in means['%s-%s-%s' % (virus,cond,img_id)]]
                
                #if len(session_responses)>=NMIN_SESSIONS:
                pt.plot(ep.t, 
                        np.mean(session_responses,axis=0)- (np.mean(session_responses,axis=0)[baselineCond].mean()),
                        sy=sem(session_responses, axis=0),
                        color=color, ax=AX[i][j])
                        
                pt.annotate(AX[i][j],
                            'N=%i' % len(session_responses)+k*'\n',
                            (0.2,0.6), ha='right',
                            color=color, fontsize=6)
        if i==0:
             pt.annotate(AX[i][j], cond, (0.5, 1))
        if j==0:
             pt.annotate(AX[i][j], 'Image-ID= %s ' % img_id,
                         (0,1), ha='right')

        pt.set_plot(AX[i][j], 
                xlabel='time (s)' if i==2 else '',
                ylabel='$\\Delta$F/F' if j==0 else '')
pt.set_common_ylims(AX)  
#%%


firgurename = 'pie_natimgID_beh_mod_'+ neuron + '.svg'
fig, AX = pt.figure(axes=(5,2))


for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['blue','darkred']):
        for i, img_id in zip(range(5), ep.varied_parameters['Image-ID']):
        
                perc_resp_ROI = np.mean(percentages['%s-%s' % (virus,img_id)],axis=0)
                rest = 100-perc_resp_ROI
                print(perc_resp_ROI,'%s-%s' % (virus,img_id))
                pt.pie(data=[perc_resp_ROI,rest],
                        COLORS=[color,'grey'],
                        ext_labels = ['resp.','non-resp.'],
                        
                        pie_labels = ['%.1f'%perc_resp_ROI,'%.1f'%rest],
                        title = '%s-%s' % (virus,img_id),
                        ax=AX[k][i])
                #pt.annotate(AX[k][i],)
                                
#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')
   


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

# %%
t_cond = (ep.t > 0.2) & (ep.t < 1.5)
# %%
resp_win = (ep.t > 0.2) & (ep.t < 1.5)

CVs = {}
Fanos = {}

for virus in ['sgRosa', 'sgCnr1']:

    CVs[virus] = {}
    Fanos[virus] = {}

    for img_id in [1.,2.,3.,4.,5.]:

        CVs[virus][img_id] = []
        Fanos[virus][img_id] = []

        sessions = means[f'{virus}-all-{img_id}']

        for data in sessions:

            # shape:
            # (trials, neurons, time)

            nROIs = data.shape[1]

            for roi in range(nROIs):

                # one scalar response per trial
                trial_responses = np.mean(
                    data[:, roi, :][:, resp_win],
                    axis=1
                )

                mu = np.mean(trial_responses)

                # avoid unstable ratios
                if np.abs(mu) < 1e-4:
                    continue

                sigma = np.std(trial_responses)

                cv = sigma / np.abs(mu)

                fano = np.var(trial_responses) / np.abs(mu)

                CVs[virus][img_id].append(cv)
                Fanos[virus][img_id].append(fano)
#%%
fig, ax = plt.subplots(1,5, figsize=(14,3), sharey=True)

all_vals = []

for virus in ['sgRosa','sgCnr1']:
    for img_id in [1.,2.,3.,4.,5.]:
        all_vals.extend(CVs[virus][img_id])

ymax = np.percentile(all_vals, 95)

for i, img_id in enumerate([1.,2.,3.,4.,5.]):

    rosa = np.array(CVs['sgRosa'][img_id])
    cnr1 = np.array(CVs['sgCnr1'][img_id])

    # bars
    ax[i].bar(
        0,
        np.mean(rosa),
        yerr=sem(rosa),
        color='grey',
        alpha=0.6
    )

    ax[i].bar(
        1,
        np.mean(cnr1),
        yerr=sem(cnr1),
        color='darkred',
        alpha=0.6
    )

    # scatter
    ax[i].scatter(
        np.random.normal(0,0.05,len(rosa)),
        rosa,
        s=10,
        color='black',
        alpha=0.5
    )

    ax[i].scatter(
        1 + np.random.normal(0,0.05,len(cnr1)),
        cnr1,
        s=10,
        color='black',
        alpha=0.5
    )

    ax[i].set_xticks([0,1])
    ax[i].set_xticklabels(['sgRosa','sgCnr1'])

    ax[i].set_title(f'Image {img_id}')

    ax[i].set_ylim([0, ymax])

    if i == 0:
        ax[i].set_ylabel('CV')

plt.tight_layout()
# %%
fig, ax = plt.subplots(1,5, figsize=(14,3), sharey=True)

all_vals = []

for virus in ['sgRosa','sgCnr1']:
    for img_id in [1.,2.,3.,4.,5.]:
        all_vals.extend(Fanos[virus][img_id])

ymax = np.percentile(all_vals, 95)

for i, img_id in enumerate([1.,2.,3.,4.,5.]):

    rosa = np.array(Fanos['sgRosa'][img_id])
    cnr1 = np.array(Fanos['sgCnr1'][img_id])

    ax[i].bar(
        0,
        np.mean(rosa),
        yerr=sem(rosa),
        color='grey',
        alpha=0.6
    )

    ax[i].bar(
        1,
        np.mean(cnr1),
        yerr=sem(cnr1),
        color='darkred',
        alpha=0.6
    )

    ax[i].scatter(
        np.random.normal(0,0.05,len(rosa)),
        rosa,
        s=10,
        color='black',
        alpha=0.5
    )

    ax[i].scatter(
        1 + np.random.normal(0,0.05,len(cnr1)),
        cnr1,
        s=10,
        color='black',
        alpha=0.5
    )

    ax[i].set_xticks([0,1])
    ax[i].set_xticklabels(['sgRosa','sgCnr1'])

    ax[i].set_title(f'Image {img_id}')

    ax[i].set_ylim([0, ymax])

    if i == 0:
        ax[i].set_ylabel('Fano factor')

plt.tight_layout()
# %%
