# %%
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


# TO CHECK A SINGLE SESSION

filename = os.path.join(os.path.expanduser('~'), 'work', 'DATA', 
                        'PN', 'NWBs',
                        '2025_09_17-15-14-34.nwb'
                        )



data = physion.analysis.read_NWB.Data(filename,
                                     verbose=False)

# custom dFoF:

from physion.imaging.Calcium import compute_F0


dFoF_options = dict(\
    method_for_F0='sliding_minmax',
    #method_for_F0='sliding_percentile',
    #method_for_F0='percentile', #more strict
    sliding_window= 60,
    percentile=10,
    roi_to_neuropil_fluo_inclusion_factor=1.,
    neuropil_correction_factor=0.7, 
    with_computed_neuropil_fact=True)

# %%
# we first perform the dFoF determination with the above params
#    (this restrict the available ROIs in the future)
data.build_dFoF(**dFoF_options, verbose=True)
#print(data.correctedFluo)
valid = data.valid_roiIndices
rejected = [i for i in range(data.Fluorescence.data.shape[1]) if (i not in valid)]


data.build_pupil_diameter()
data.build_facemotion()
data.build_running_speed()


# %%



# BUILDING --- DRIFTING GRATINGS ---- EPISODE 
epGrating = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF',
                                                                'running_speed', 'pupil_diameter', 'facemotion'],
                                                    protocol_name='drifting-grating')


# %%
# PLOT PROPERTIES --- DRIFTING GRATINGS ---

plot_props = dict(column_key='contrast',
                  with_annotation=True,
                  Ybar=0.5, Ybar_label="0.5$\Delta$F/F",
                  figsize=(9,1.8))
#%%


#RESPONSES OF --- DRIFTING GRATINGS --- SINGLE ROIs single session 

for i in range(epGrating.data.nROIs):
        roiIndex=i,
        fig, AX = physion.dataviz.episodes.trial_average.plot(epGrating, with_std=False,
                                                        roiIndex=roiIndex,
                                                        **plot_props)
        pt.show()



#%%

#RESPONSES OF --- DRIFTING GRATINGS --- AVERAGE OF ONE SESSION  single session
print(3)
fig, AX = physion.dataviz.episodes.trial_average.plot(epGrating,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)

#%%








#%%


# BUILDING --- NATURAL IMAGES ---- EPISODE 

epNatImg = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF'],
                                                    protocol_name='Natural-Images-4-repeats')

#%%
# PLOT PROPERTIES --- NATURAL IMAGES ---
plot_props = dict(column_key='Image-ID',
                  with_annotation=True,
                  figsize=(9,1.8))




#%%


#RESPONSES OF ---NATURAL IMAGES --- SINGLE ROIs - single session 

for i in range(epNatImg.data.nROIs):
        roiIndex=i,
        fig, AX = physion.dataviz.episodes.trial_average.plot(epNatImg, with_std=False,
                                                        roiIndex=roiIndex,
                                                        **plot_props)
        pt.show()



#%%

#RESPONSES OF --- NATURAL IMAGES --- AVERAGE OF ONE SESSION - single session

fig, AX = physion.dataviz.episodes.trial_average.plot(epNatImg,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)



                                        
# %%

#RESPONSES OF --- NATURAL IMAGES --- AVERAGE OF ONE SESSION - single session with SEM

fig, AX = pt.figure(axes=(5,1))
for i, mean in enumerate(means):
        print(np.mean(mean))
        pt.plot(epNatImg.t, np.mean(mean, axis=0), 
                sy=stats.sem(mean, axis=0), ax=AX[i])
        #AX[i].axis('off')
#pt.draw_bar_scales(AX[0], Xbar=1, Xbar_label='1s', Ybar=0.1, Ybar_label='0.1dFoF')#

pt.set_common_ylims(AX)
        

#%%




#%%

# BUILDING --- LOOMING STIMULUS ---- EPISODE 

epLooming = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF',
                                                                'running_speed', 'pupil_diameter', 'facemotion'],
                                                    protocol_name='looming-stim')

#%%

# PLOT PROPERTIES --- LOOMING STIMULUS ---

plot_props = dict(
                  with_annotation=True,
                  figsize=(1,1))
#%%

# RESPONSES OF --- LOOMING STIMULUS --- SINGLE ROIs - single session 

for i in range(epLooming.data.nROIs):
        roiIndex=i,
        fig, AX = physion.dataviz.episodes.trial_average.plot(epLooming, with_std=False,
                                                      roiIndex=roiIndex,
                                                      **plot_props)


#%%

#RESPONSES OF --- LOOMING STIMULUS --- AVERAGE OF ONE SESSION of Calcium, pupildiameter and running speed - single session


fig, AX = physion.dataviz.episodes.trial_average.plot(epLooming,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)
pt.show()



fig, AX = physion.dataviz.episodes.trial_average.plot(epLooming,
                                                      quantity='pupil_diameter', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)
pt.show()


fig, AX = physion.dataviz.episodes.trial_average.plot(epLooming,
                                                      quantity='running_speed', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)
pt.show()


fig, AX = physion.dataviz.episodes.trial_average.plot(epLooming,
                                                      quantity='facemotion', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)
pt.show()



#%%




#%%

# BUILDING --- STATIC PATCH ---- EPISODE 

epStaticPatch = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF',
                                                                'running_speed'],
                                                    protocol_name='static-patch')

#%%
# PLOT PROPERTIES --- STATIC PATCH ---

plot_props = dict(column_key='angle',
                  with_annotation=True,
                  figsize=(9,1.8))

#%%

# RESPONSES OF --- STATIC PATCH --- SINGLE ROIs - single session 

for i in range(epStaticPatch.data.nROIs):
        roiIndex=i,
        fig, AX = physion.dataviz.episodes.trial_average.plot(epStaticPatch, with_std=False,
                                                      roiIndex=i,
                                                      **plot_props)

#%%
#RESPONSES OF --- STATIC PATCH --- AVERAGE OF ONE SESSION

fig, AX = physion.dataviz.episodes.trial_average.plot(epStaticPatch,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)

#%%

#%%


# BUILDING --- INVERSE RECEPTIVE FIELD ---- EPISODE 

epInverseRF = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF', 
                                                                'running_speed'],
                                                    protocol_name='drifting-surround')

#%%

# PLOT PROPERTIES --- INVERSE RECEPTIVE FIELD ---

plot_props = dict(
                  with_annotation=True,
                  figsize=(1,1))


#%%


# RESPONSES OF --- INVERSE RECEPTIVE FIELD --- SINGLE ROIs - single session 

for i in range(epInverseRF.data.nROIs):
        roiIndex = i
        fig, AX = physion.dataviz.episodes.trial_average.plot(epInverseRF, with_std=False,
                                                      roiIndex=i,
                                                      **plot_props)


#%%
#RESPONSES OF --- INVERSE RECEPTIVE FIELD --- AVERAGE OF ONE SESSION

fig, AX = physion.dataviz.episodes.trial_average.plot(epInverseRF,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)


#%%



# BUILDING --- MOVING DOTS ---- EPISODE 

epMovingDots = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF', 
                                                                'running_speed'],
                                                    protocol_name='moving-dots')

#%%

# PLOT PROPERTIES --- MOVING DOTS ---

plot_props = dict(
                  with_annotation=True,
                  figsize=(1,1))

#%%

# RESPONSES OF --- MOVING DOTS --- SINGLE ROIs - single session 


for i in range(epMovingDots.data.nROIs):
        roiIndex=i,
        fig, AX = physion.dataviz.episodes.trial_average.plot(epMovingDots, with_std=False,
                                                      roiIndex=roiIndex,
                                                      **plot_props)


#%%
#RESPONSES OF --- INVERSE RECEPTIVE FIELD --- AVERAGE OF ONE SESSION

fig, AX = physion.dataviz.episodes.trial_average.plot(epMovingDots,
                                                      quantity='dFoF', with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)


#%%



#%%
plot_props = dict(
                  with_annotation=True,
                  figsize=(1,1))

epBlack2Min = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['dFoF', 
                                                                'running_speed'],
                                                    protocol_name='black-2min')

fig, AX = physion.dataviz.episodes.trial_average.plot(epBlack2Min,
                                                      quantity='dFoF', 
                                                                 with_std=False,
                                                      roiIndices='all',
                                                      **plot_props)


#%%

for i in range(epBlack2Min.data.nROIs):
        roiIndex=i,
        fig, AX = physion.dataviz.episodes.trial_average.plot(epBlack2Min, with_std=False,
                                                      roiIndex=roiIndex,
                                                      **plot_props)
        
#%%


#%%

# TO LOOP OVER NWB FILES WITH VISUAL STIMULUS --- DRIFITING GRATING ---  multisession

folder = os.path.join(os.path.expanduser('~'), 'work', 'DATA', 
                        'PN', 'NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol='drifting-grating')


means = [[] for i in range(3)] # array over contrasts
for i, filename in enumerate(DATASET['files']):
    
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    print(i+1, '--', filename, '--', data.nROIs)
    print(data.protocols)

    data.build_dFoF(**dFoF_options, verbose=True)
    data.build_pupil_diameter()
    data.build_facemotion()
    data.build_running_speed()

    if data.nROIs>0:

        epGrating = physion.analysis.process_NWB.EpisodeData(data, 
                                                        quantities=['dFoF'],
                                                                        
                                                        protocol_name='drifting-grating')
        
        if 'contrast' in epGrating.varied_parameters:
                fig, AX = physion.dataviz.episodes.trial_average.plot(epGrating,
                                                                quantity='dFoF', with_std=False,
                                                                roiIndices='all',
                                                                **plot_props)
                for i in range(3):
                        contrast_cond = epGrating.find_episode_cond(key='contrast', index=i)
                        means[i].append(epGrating.dFoF[contrast_cond,:,:].mean(axis=(0,1)))
                pt.show()     
    else:
           print(' !!!!!!  ', filename)


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
