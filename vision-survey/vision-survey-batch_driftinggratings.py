
import numpy as np
import os, sys
sys.path.append('../physion/src') # add src code directory for physion
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
    with_computed_neuropil_fact=True)

# PLOT PROPERTIES --- DRIFTING GRATINGS ---
stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='ttest')


# TO LOOP OVER NWB FILES WITH VISUAL STIMULUS --- DRIFITING GRATING ---  multisession

folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna',
                        'PN_cond-NDNF-CB1_WT-vs-KD', 'NWBs')

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

means = {} # 
for key in ['sgRosa', 'sgCnr1']:
      means['all-%s' % key] = [[] for i in range(3)]
      means['run-%s' % key] = [[] for i in range(3)]
      means['still-%s' % key] = [[] for i in range(3)]

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
# 

from scipy.stats import sem

fig, AX = pt.figure(axes=(3,3))

for j, cond in enumerate(['all', 'run', 'still']):
    for i, c in zip(range(3), epGrating.varied_parameters['contrast']):
        for k, key, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['r','b']):
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


#%% 
# Sort by contrast and behavior
# 
from scipy.stats import sem

fig, AX = pt.figure(axes=(3,3))

for j, cond in enumerate(['all', 'run', 'still']):
    for i, c in zip(range(3), epGrating.varied_parameters['contrast']):
        for k, key, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['r','b']):
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

#%%


# TO LOOP OVER NWB FILES WITH VISUAL STIMULUS --- NATURAL-IMAGES ---  multisession
DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol='Natural-Images-4-repeats')


means = [[] for i in range(5)] # array over condition
for i, filename in enumerate(DATASET['files']):
    
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    print(i+1, '--', filename, '--', data.nROIs)
    print(data.protocols)

    data.build_dFoF(**dFoF_options, verbose=True)
    #data.build_rawFluo()#(**dFoF_options, verbose=True)

    if data.nROIs>0:

        epNatImg = physion.analysis.process_NWB.EpisodeData(data, 
                                                        quantities=['dFoF'],
                                                        protocol_name='Natural-Images-4-repeats')
        
        if 'Image-ID' in epNatImg.varied_parameters:
                fig, AX = physion.dataviz.episodes.trial_average.plot(epNatImg,
                                                                quantity='dFoF', with_std=False,
                                                                roiIndices='all',
                                                                **plot_props)
                for i in range(5):
                        image_cond = epNatImg.find_episode_cond(key='Image-ID', index=i)
                        means[i].append(epNatImg.dFoF[image_cond,:,:].mean(axis=(0,1)))     
    else:
           print(' !!!!!!  ', filename)


#%%
"""#build z score per episode - CURRENTLY ON PAUSE
bsl = epGrating.dFoF[:, :, 0:1998].mean(axis=(0,1)) 


for i in range(epGrating.dFoF.nROIs):
        
        pt.plot(epGrating.dFoF[i, :, :])

        
        
from scipy import stats

fig, AX = pt.figure(axes=(3,1))
for i, mean in enumerate(means):
        print(np.mean(mean))
        pt.plot(epGrating.t, np.mean(mean, axis=0), 
                sy=stats.sem(mean, axis=0), ax=AX[i])
        #AX[i].axis('off')
#pt.draw_bar_scales(AX[0], Xbar=1, Xbar_label='1s', Ybar=0.1, Ybar_label='0.1dFoF')#

pt.set_common_ylims(AX)
"""




#%%
# First, plot responses per ROI of three contrasts dFoF
"""
rois = epGrating.dFoF

for i in data.nROIs:
        
        
        fig, ax = pt.figure()
        cond1 = epGrating.find_episode_cond(key='contrast', value= 0.2)
        
        
        pt.plot(epGrating.t, epGrating.dFoF[cond1,i,:].mean(axis=1), ax=ax)
        pt.set_plot(ax, 
            xlabel='time from stim. (s)',
            ylabel='dFoF',
            title= "mean dFoF of ROI:"+str(i))
        

#%%
pt.plot(epGrating.t,epGrating.dFoF[6,19,:])
pt.set_plot(ax, 
xlabel='time from stim. (s)',
ylabel='dFoF',
title= "mean dFoF of ROI:")




#%%

# First, plot responses per ROI of three contrasts dFoF

fig, ax = pt.figure()
cond1 = epGrating.find_episode_cond(key='contrast', value=0.2)


pt.plot(epGrating.t, epGrating.dFoF[cond1,i-1,:].mean(axis=0), ax=ax)
pt.set_plot(ax, 
        xlabel='time from stim. (s)',
        ylabel='dFoF',
        title= "mean dFoF of ROI:"+str(i))

#sy=np.std(epGrating.dFoF.mean(axis=0)
#%% 
#separate by arousal state


#build_pupildiameter


#%%
epNatImg = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF',
                                                                'running_speed'],
                                                    protocol_name='Natural-Images-4-repeats')


roiIndex = 1
cond1 = epNatImg.find_episode_cond(key='Image-ID', value=1)
fig, ax = pt.figure()
pt.plot(epNatImg.t, 
        epNatImg.dFoF[cond1,roiIndex,:].mean(axis=0), ax=ax)
pt.set_plot(ax, 
            xlabel='time from stim. (s)',
            ylabel='dFoF')

#sy=np.std(epGrating.dFoF.mean(axis=0),
# %%
fig, ax = pt.figure()
cond = data.pupil_diameter<0.9
data.pupil_diameter[cond] = 1
pt.plot(data.pupil_diameter, ax=ax)
#pt.set_plot(ax, xlim=(1400,1600))

# %%
epGrating = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF',
                                                                'running_speed'],
                                                    protocol_name='drifting-grating')"""

# %%
"""print(epGrating.varied_parameters)

# %%
from scipy import stats
fig, ax = pt.figure()
pt.plot(epGrating.t, 
        epGrating.running_speed.mean(axis=0),
        sy=stats.sem(epGrating.running_speed, axis=0), ax=ax)
pt.set_plot(ax, 
            xlabel='time from stim. (s)',
            ylabel='run')
# %%
roiIndex = 10
from scipy import stats
fig, ax = pt.figure()
pt.plot(epGrating.t, 
        # epGrating.dFoF[:,roiIndex,:].mean(axis=0),
        epGrating.dFoF.mean(axis=(0,1)),
        sy=np.std(epGrating.dFoF.mean(axis=0), axis=0))
pt.set_plot(ax, 
            xlabel='time from stim. (s)',
            ylabel='dFoF')

# %%
fig, ax = pt.figure()
cond1 = epGrating.find_episode_cond(key='contrast', value=0.2)
pt.plot(epGrating.dFoF[cond1, :, :].mean(axis=(0,1)), 
        color='grey', ax=ax)
cond2 = epGrating.find_episode_cond(key='contrast', value=1.0)
pt.plot(epGrating.dFoF[cond2, :, :].mean(axis=(0,1)), 
        ax=ax)

# %%
EPISODES = []
for p, protocol in enumerate(data.protocols):
    print(p, protocol)
    EPISODES.append(\
        physion.analysis.process_NWB.EpisodeData(data, 
                                    quantities=['dFoF',
                                                'running_speed'],
                                    protocol_name=protocol))


# %%
fig, AX = pt.figure(axes=(1, len(data.protocols)), hspace=2.)
for p, protocol in enumerate(data.protocols):
    pt.plot(EPISODES[p].t, EPISODES[p].dFoF.mean(axis=(0,1)),
            ax=AX[p])
    pt.set_plot(AX[p], 
                xlabel='time (s)', ylabel='dFoF',
                title=data.protocols[p])


# %%
EPISODES[8].find_episode_cond(key='x-center', 0)"""

# %%
x