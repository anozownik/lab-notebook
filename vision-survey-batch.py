
import numpy as np
import os, sys
sys.path.append('physion/src') # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')
from scipy import stats

    
# %%

"""" VISION SURVEY SUBPROTOCOLS 
static-patch' 'looming-stim' 'Natural-Images-4-repeats'
'drifting-grating' 'drifting-surround' 'moving-dots' 'grey-10min'
'black-2min' 'quick-spatial-mapping']"""

# custom dFoF:

from physion.imaging.Calcium import compute_F0


dFoF_options = dict(\
    method_for_F0='sliding_minmax',
    #method_for_F0='sliding_percentile',
    #method_for_F0='percentile', #more strict
    sliding_window= 300,
    percentile=10,
    roi_to_neuropil_fluo_inclusion_factor=1.,
    neuropil_correction_factor=0.7, 
    with_computed_neuropil_fact=False)

# PLOT PROPERTIES --- DRIFTING GRATINGS ---
stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='ttest')


# TO LOOP OVER NWB FILES WITH VISUAL STIMULUS --- DRIFITING GRATING ---  multisession

folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna',
                        'NDNF_cond-CB1_WT-vs-KD', 'NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol='drifting-grating')



# STATISTICS PROPERTIES --- DRIFTING GRATINGS ---
stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='ttest')


# PLOT PROPERTIES --- DRIFTING GRATINGS ---

plot_props = dict(column_key='contrast',
                  with_annotation=True,
                  Ybar=0.5, Ybar_label="0.5$\Delta$F/F",
                  Xbar=0.5, Xbar_label="0.5s",
                  figsize=(9,1.8))


# RESPONSE ARGUMENTS --- DRIFTING GRATINGS ---



response_args = dict(quantity='dFoF')

summary_stats = []

RUNNING_SPEED_THRESHOLD = 0.1

PUPIL_THRESHOLD = 2.55

means = {} # 
for key in ['sgRosa', 'sgCnr1']:
      means['all-%s' % key] = [[] for i in range(3)]
      means['run-%s' % key] = [[] for i in range(3)]
      means['still-%s' % key] = [[] for i in range(3)]
      means['all_pupil-%s' % key] = [[] for j in range(3)]
      means['high_arousal-%s' % key] = [[] for j in range(3)]
      means['low_arousal-%s' % key] = [[] for j in range(3)]

for i, filename in enumerate(DATASET['files']):
    
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    print(i+1, '--', filename, '--', data.nROIs)
    print(data.protocols)

    data.build_dFoF(**dFoF_options, verbose=True)
    data.build_pupil_diameter()
    data.build_facemotion()
    data.build_running_speed()
    
    if 'sgRosa' in data.nwbfile.virus:
        color = 'grey'
        key = 'sgRosa'
    elif 'sgCnr1':
        color = 'darkred'
        key = 'sgCnr1'
        

    if data.nROIs>0:

        epGrating = physion.analysis.process_NWB.EpisodeData(data, 
                                                        quantities=['dFoF', 'running_speed'],
                                                        protocol_name='drifting-grating')
        
        withinEpisode = (epGrating.t>0) & (epGrating.t<epGrating.time_duration[0])
        run = np.mean(epGrating.running_speed[:,withinEpisode], axis=1) > RUNNING_SPEED_THRESHOLD
        arousal_state = np.mean(epGrating.pupil_diameter[:,withinEpisode], axis=1) > PUPIL_THRESHOLD
        
        if 'contrast' in epGrating.varied_parameters:
                fig, AX = physion.dataviz.episodes.trial_average.plot(epGrating,
                                                                quantity='dFoF', with_std=False, with_stat_test=True, stat_test_props=stat_test_props,
                                                                color=color,
                                                                roiIndices='all',
                                                                **plot_props)
                for i in range(3):
                        contrast_cond = epGrating.find_episode_cond(key='contrast', index=i)
                        means['all-%s' % key][i].append(epGrating.dFoF[contrast_cond,:,:].mean(axis=(0,1)))
                        if np.sum(contrast_cond & run)>=2:
                            means['run-%s' % key][i].append(epGrating.dFoF[contrast_cond & run,:,:].mean(axis=(0,1)))
                        if np.sum(contrast_cond & ~run)>=2:
                            means['still-%s' % key][i].append(epGrating.dFoF[contrast_cond & ~run,:,:].mean(axis=(0,1)))

                for j in range(3):
                        contrast_cond = epGrating.find_episode_cond(key='contrast', index=j)
                        means['all_pupil-%s' % key][j].append(epGrating.dFoF[contrast_cond,:,:].mean(axis=(0,1)))
                        if np.sum(contrast_cond & arousal_state)>=2:
                            means['high_arousal-%s' % key][j].append(epGrating.dFoF[contrast_cond & arousal_state,:,:].mean(axis=(0,1)))
                        if np.sum(contrast_cond & ~arousal_state)>=2:
                            means['low_arousal-%s' % key][j].append(epGrating.dFoF[contrast_cond & ~arousal_state,:,:].mean(axis=(0,1)))
                pt.show()

                result = epGrating.compute_summary_data( stat_test_props=stat_test_props,
                        response_args=response_args,
                             response_significance_threshold=0.01,
                             verbose=True)
                print(result)
                summary_stats.append(result)
        
    else:
           print(' !!!!!!  ', filename)


            


        

#%%
from scipy.stats import sem

fig, AX = pt.figure(axes=(3,3))

for j, cond in enumerate(['all', 'run', 'still']):
    for i, c in zip(range(3), epGrating.varied_parameters['contrast']):
        for k, key, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','r']):
            if len(means['%s-%s' % (cond, key)][i])>1:
                pt.plot(epGrating.t, 
                        np.mean(means['%s-%s' % (cond, key)][i], axis=0),
                        sy=sem(means['%s-%s' % (cond, key)][i], axis=0),
                        color=color, ax=AX[i][j])
                pt.annotate(AX[i][j],
                            'N=%i' % len(means['%s-%s' % (cond, key)][i])+k*'\n',
                            (0,0), #ha='right',
                            color=color, fontsize=6)
        if i==0:
             pt.annotate(AX[i][j], cond, (0.5, 1))
        if j==0:
             pt.annotate(AX[i][j], 'contrast=%.1f ' % c,
                         (0,1), ha='right')

        pt.set_plot(AX[i][j], 
                    xlabel='time (s)' if i==2 else '',
                    ylabel='$\\Delta$F/F' if j==0 else '')
pt.set_common_ylims(AX)


# %%

fig, AX = pt.figure(axes=(3,3))

for j, cond in enumerate(['all_pupil', 'high_arousal', 'low_arousal']):
    for i, c in zip(range(3), epGrating.varied_parameters['contrast']):
        for k, key, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','r']):
            if len(means['%s-%s' % (cond, key)][i])>1:
                pt.plot(epGrating.t, 
                        np.mean(means['%s-%s' % (cond, key)][i], axis=0),
                        sy=sem(means['%s-%s' % (cond, key)][i], axis=0),
                        color=color, ax=AX[i][j])
                pt.annotate(AX[i][j],
                            'N=%i' % len(means['%s-%s' % (cond, key)][i])+k*'\n',
                            (0,0), #ha='right',
                            color=color, fontsize=6)
        if i==0:
             pt.annotate(AX[i][j], cond, (0.5, 1))
        if j==0:
             pt.annotate(AX[i][j], 'contrast=%.1f ' % c,
                         (0,1), ha='right')

        pt.set_plot(AX[i][j], 
                    xlabel='time (s)' if i==2 else '',
                    ylabel='$\\Delta$F/F' if j==0 else '')
pt.set_common_ylims(AX)

# %%
