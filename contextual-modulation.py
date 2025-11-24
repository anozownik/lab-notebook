# %%
import numpy as np
import os, sys
sys.path.append('physion/src') # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')
from scipy import stats


# %%
filename = os.path.join(os.path.expanduser('~'), 'work', 'DATA', 
                        'NDNF', 'NWBs',
                        '2025_10_06-17-01-40.nwb'
                        )

#2025_02_26-11-53-54.nwb '2025_05_06-18-14-23 '2025_03_19-17-14-07.nwb' '2025_04_09-16-07-54.nwb' '2025_05_06-16-28-19.nwb' '2025_06_25-12-14-38.nwb'
#'2025_03_18-14-45-26.nwb' '2025_03_19-15-42-35.nwb' '2025_03_05-15-32-48.nwb' '2025_04_01-18-04-01.nwb' '2025_06_25-15-19-24.nwb'

data = physion.analysis.read_NWB.Data(filename,
                                     verbose=False)

# custom dFoF:

from physion.imaging.Calcium import compute_F0


dFoF_options = dict(\
    method_for_F0='hamming',
    sliding_window= 0.5,
    percentile=10,
    roi_to_neuropil_fluo_inclusion_factor=1.15,
    neuropil_correction_factor=0.7, 
    with_computed_neuropil_fact=False)

# %%
# we first perform the dFoF determination with the above params
#    (this restrict the available ROIs in the future)
data.build_dFoF(**dFoF_options, verbose=True)
#print(data.correctedFluo)
valid = data.valid_roiIndices
rejected = [i for i in range(data.Fluorescence.data.shape[1]) if (i not in valid)]
dFoF_options = dict()
#data.build_pupil_diameter()

# %%


data.build_Zscore_dFoF(verbose =True)




# %%
coord_map = {
            (np.float64(-36.0), np.float64(-23.0)): (2, 0),
            (np.float64(-36.0), np.float64(0.0)):   (2, 1),
            (np.float64(-36.0), np.float64(23.0)):  (2, 2),
            (np.float64(0.0), np.float64(-23.0)):   (1, 0),
            (np.float64(0.0), np.float64(0.0)):     (1, 1),
            (np.float64(0.0), np.float64(23.0)):    (1, 2),
            (np.float64(36.0), np.float64(-23.0)):  (0, 0),
            (np.float64(36.0), np.float64(0.0)):    (0, 1),
            (np.float64(36.0), np.float64(23.0)):   (0, 2),
        }



#%%
epGrating = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF', 'Zscore_dFoF',
                                                                'running_speed'],
                                                    protocol_name='drifting-grating')
#%%

folder = os.path.join(os.path.expanduser('~'), 'work', 'DATA', 
                        'NDNF', 'NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol='drifting-grating')

#%%

plot_props = dict(column_key='contrast',
                  with_annotation=True, figsize=(9,1.8))
#%%
for filename in DATASET['files']:
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    print(filename)
    epGrating = physion.analysis.process_NWB.EpisodeData(data, 
                                                quantities=['dFoF',
                                                                'running_speed'],
                                                protocol_name='drifting-grating')
    fig, AX = physion.dataviz.episodes.trial_average.plot(epGrating,
                                                quantity='dFoF', with_std=False,
                                                roiIndices='all',
                                                **plot_props)
    

                                        



#%%
"""
for filename in DATASET['files']:
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    ep = physion.analysis.process_NWB.EpisodeData(data, 
                                              prestim_duration=4,
                                              quantities=['pupil_diameter'])
    fig, axs = pt.figure(axes_extents=[[[1,1],[1,1]]])

    for ax, val, c in zip(axs, [3,6], ['grey','white']):
        cond = ep.find_episode_cond(key='Image-ID', value=val)
        pt.plot(ep.t, 
                ep.pupil_diameter[cond, :].mean(axis=(0)),
                sy=stats.sem(ep.pupil_diameter[cond, :], axis=(0)),
                color=c, ax=ax)
        pt.set_plot(ax=ax,
            xlabel='time from stim. (s)',
            ylabel='pupil diam. (mm)')

    pt.show()
"""

#%%
plot_props = dict(column_key='contrast',
                  with_annotation=True,
                  figsize=(9,1.8))


#%%

fig, AX = physion.dataviz.episodes.trial_average.plot(epGrating,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)


#%%

for i in epGrating.data.valid_roiIndices:
        fig, AX = physion.dataviz.episodes.trial_average.plot(epGrating, with_std=False,
                                                      roiIndex=i,
                                                      **plot_props)



#%%

plot_props = dict(column_key='Image-ID',
                  with_annotation=True,
                  figsize=(9,1.8))

epNatImg = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF', 'Zscore_dFoF',
                                                                'running_speed'],
                                                    protocol_name='Natural-Images-4-repeats')

fig, AX = physion.dataviz.episodes.trial_average.plot(epNatImg,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)


#%%

for i in epNatImg.data.valid_roiIndices:
        fig, AX = physion.dataviz.episodes.trial_average.plot(epNatImg, with_std=False,
                                                      roiIndex=i,
                                                      **plot_props)
        
#%%
""""static-patch' 'looming-stim' 'Natural-Images-4-repeats'
 'drifting-grating' 'drifting-surround' 'moving-dots' 'grey-10min'
 'black-2min' 'quick-spatial-mapping']"""

plot_props = dict(column_key='angle',
                  with_annotation=True,
                  figsize=(9,1.8))

epStaticPatch = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF', 'Zscore_dFoF',
                                                                'running_speed'],
                                                    protocol_name='static-patch')

fig, AX = physion.dataviz.episodes.trial_average.plot(epStaticPatch,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)


#%%

for i in epStaticPatch.data.valid_roiIndices:
        fig, AX = physion.dataviz.episodes.trial_average.plot(epStaticPatch, with_std=False,
                                                      roiIndex=i,
                                                      **plot_props)
        

#%%
plot_props = dict(
                  with_annotation=True,
                  figsize=(1,1))

epLooming = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF', 'Zscore_dFoF',
                                                                'running_speed'],
                                                    protocol_name='looming-stim')

fig, AX = physion.dataviz.episodes.trial_average.plot(epLooming,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)


#%%

for i in epLooming.data.valid_roiIndices:
        fig, AX = physion.dataviz.episodes.trial_average.plot(epLooming, with_std=False,
                                                      roiIndex=i,
                                                      **plot_props)

#%%

plot_props = dict(
                  with_annotation=True,
                  figsize=(1,1))

epInverseRF = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF', 'Zscore_dFoF',
                                                                'running_speed'],
                                                    protocol_name='drifting-surround')

fig, AX = physion.dataviz.episodes.trial_average.plot(epInverseRF,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)


#%%

for i in epLooming.data.valid_roiIndices:
        fig, AX = physion.dataviz.episodes.trial_average.plot(epLooming, with_std=False,
                                                      roiIndex=i,
                                                      **plot_props)

#%%

plot_props = dict(
                  with_annotation=True,
                  figsize=(1,1))

epMovingDots = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF', 'Zscore_dFoF',
                                                                'running_speed'],
                                                    protocol_name='moving-dots')

fig, AX = physion.dataviz.episodes.trial_average.plot(epMovingDots,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)


#%%

for i in epLooming.data.valid_roiIndices:
        fig, AX = physion.dataviz.episodes.trial_average.plot(epMovingDots, with_std=False,
                                                      roiIndex=i,
                                                      **plot_props)
        

#%%
plot_props = dict(
                  with_annotation=True,
                  figsize=(1,1))

epBlack2Min = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF', 'Zscore_dFoF',
                                                                'running_speed'],
                                                    protocol_name='black-2min')

fig, AX = physion.dataviz.episodes.trial_average.plot(epBlack2Min,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)


#%%

for i in epLooming.data.valid_roiIndices:
        fig, AX = physion.dataviz.episodes.trial_average.plot(epBlack2Min, with_std=False,
                                                      roiIndex=i,
                                                      **plot_props)






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
