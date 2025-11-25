# %%
import numpy as np
import matplotlib.pyplot as plt
import os, sys
sys.path.append('physion/src') # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('white')
from scipy import stats


folder = os.path.join(os.path.expanduser('~'), 'work', 'DATA', 
                        'Asahi', 'NWBs')
filename = os.path.join(folder, '2025_09_16-15-29-43.nwb')

data = physion.analysis.read_NWB.Data(filename,
                                     verbose=False)

# custom dFoF:
#dFoF_options = dict()
data.build_pupil_diameter()



# %%
fig, ax = pt.figure()
cond = data.pupil_diameter<0.9 # everything below 0.9 pupil diameter is not considered
data.pupil_diameter[cond] = 1
pt.plot( data.pupil_diameter, ax=ax)
pt.set_plot(ax, xlim=(1400,1600))

# %%
ep = physion.analysis.process_NWB.EpisodeData(data, 
                                              prestim_duration=4,
                                              quantities=['pupil_diameter'])

# %%
fig, ax = pt.figure()
cond1 = ep.find_episode_cond(key='Image-ID', value=3)
pt.plot(ep.t, ep.pupil_diameter[cond1, :].mean(axis=(0)),
        color='grey', ax=ax)
#cond2 = epGrating.find_episode_cond(key='contrast', value=1.0)
#pt.plot(epGrating.dFoF[cond2, :, :].mean(axis=(0,1)), 
 #       ax=ax)
pt.set_plot(ax, 
            xlabel='time from stim. (s)',
            ylabel='pupil diam. (mm)')




# %%
#start of function, goals are a) to take correct NWBs, sorted by protocol. b) build episodes as above c) plot/subplot the average pupildiameter across subjects
def func(filename):
    pass
    return fig



# %%


# BIG LOOPS, LESS ELEGANT
# selecting the NWBs by the correct protocol and putting them in an object or class
DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol='Asahi-unorthodox-SIZE-vars')



# now i have a collection of NWBs + attributes, I want to  build a parameter from them, e.g. pupil diameter
for filename in DATASET['files']:
    data = physion.analysis.read_NWB.Data(filename,
                                     verbose=False)
    data.build_pupil_diameter()
    #fig, ax = pt.figure()
    #pt.plot(data.pupil_diameter, ax=ax)


# Building an episode (ep) that considers the pupildiameter of the data and 4 seconds before    

    ep = physion.analysis.process_NWB.EpisodeData(data, 
                                              prestim_duration=4,
                                              quantities=['pupil_diameter'])
    
# Defining figure structure: list of axes determines the sizes of subplots
    fig, axs = pt.figure(axes_extents=[[[1,1],[1,1],[1,1]]])

    
    cond_medium = ep.find_episode_cond(key='Image-ID', value=3)
    pt.plot(ep.t, ep.pupil_diameter[cond_medium, :].mean(axis=(0)), 
    color='grey', ax=axs[0])
    pt.set_plot(ax=axs[0],
        xlabel='time from stim. (s)',
        ylabel='pupil diam. (mm)')    

    
    
    cond_small = ep.find_episode_cond(key='Image-ID', value=6)
    pt.plot(ep.t, ep.pupil_diameter[cond_small, :].mean(axis=(0)), 
    color='white', ax=axs[1])
    pt.set_plot(ax=axs[1], 
        xlabel='time from stim. (s)',
        ylabel='pupil diam. (mm)')
    

    cond_big = ep.find_episode_cond(key='Image-ID', value=9)
    pt.plot(ep.t, ep.pupil_diameter[cond_big, :].mean(axis=(0)), 
    color='blue', ax=axs[2])
    pt.set_plot(ax=axs[2], 
        xlabel='time from stim. (s)',
        ylabel='pupil diam. (mm)')
    pt.show()

#%%
DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol='Asahi-unorthodox-SIZE-vars')

# FOR "SIZE VARIATION" ASAHI UNORTHODOX

# now i have a collection of NWBs + attributes (?), I want to  build a parameter from them, e.g. pupil diameter
for filename in DATASET['files']:
    data = physion.analysis.read_NWB.Data(filename,
                                     verbose=False)
    data.build_pupil_diameter()
    #fig, ax = pt.figure()
    #pt.plot(data.pupil_diameter, ax=ax)


    

    ep = physion.analysis.process_NWB.EpisodeData(data, 
                                              prestim_duration=4,
                                              quantities=['pupil_diameter'])
    

    fig, axs = pt.figure(axes_extents=[[[1,1],[1,1],[1,1]]])

    for ax, val, c in zip(axs, [3,6,9], ['grey','white','blue']):
        cond = ep.find_episode_cond(key='Image-ID', value=val)
        pt.plot(ep.t, ep.pupil_diameter[cond, :].mean(axis=(0)), 
        color=c, ax=ax)
        pt.set_plot(ax=ax,
            xlabel='time from stim. (s)',
            ylabel='pupil diam. (mm)')

    pt.show()



                                                   




#%%


# FOR "NORMAL" ASAHI UNORTHODOX

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol='Asahi-unorthodox')

for filename in DATASET['files']:
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    data.build_pupil_diameter()
    fig, ax = pt.figure(ax_scale=(3,3))
    pt.plot(data.pupil_diameter, ax=ax)


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

# %%
# FOR "NORMAL" ASAHI


DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol='Asahi-classical')



for filename in DATASET['files']:
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    data.build_pupil_diameter()
    fig, ax = pt.figure(ax_scale=(3,3))
    pt.plot(data.pupil_diameter, ax=ax)
    pt.show()

    ep = physion.analysis.process_NWB.EpisodeData(data, 
                                              prestim_duration=3,
                                              quantities=['pupil_diameter'])
    fig, axs = pt.figure(axes_extents=[[[1,1],[1,1]]])

    for ax, val, c in zip(axs, [1,2], ['grey','white']):
        cond = ep.find_episode_cond(key='Image-ID', value=val)
        pt.plot(ep.t, 
                ep.pupil_diameter[cond, :].mean(axis=(0)),
                sy=stats.sem(ep.pupil_diameter[cond, :], axis=(0)),
                color=c, ax=ax)
        pt.set_plot(ax=ax,
            xlabel='time from stim. (s)',
            ylabel='pupil diam. (mm)')

    pt.show()
# %%
