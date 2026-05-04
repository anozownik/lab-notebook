import numpy as np
import os, sys
sys.path.append('./physion/src') # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')
from scipy import stats

#%%

#filename = os.path.join(os.path.expanduser('~'),  'DATA', 'Adrianna', 
#                        'PN_cond-NDNF-CB1_WT-vs-KD', 'spontaneous', 'NWBs',
#                        '2025_07_30-16-35-55.nwb'
                        
 


#data = physion.analysis.read_NWB.Data(filename,
                                     #verbose=False)


#%%

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


#%%

#data.build_dFoF(**dFoF_options, verbose=True)

#ep = physion.analysis.episodes.build.EpisodeData(data,quantities=['dFoF'], protocol_name ='imaging during spontaneous behavior')

# %%
#Matrix = np.corrcoef(data.dFoF)
# %%

# filter the matrix for unique pairs
#triu = np.array(np.triu(np.ones(Matrix.shape), k=1),
#                dtype=bool)
# 
#Npairs = np.sum(triu.flatten()) # number of unique pairs

# %%
#Matrix = np.triu(Matrix, k=1)
#pt.matrix(Matrix, vmin=0, vmax=1)

# %%
"""
fig, ax = pt.figure()
ax.hist(Matrix[triu])
pt.set_plot(ax, xlabel='corr. coef.', ylabel='# pairs',
            title= '%.2f $\pm$ %.2f ' % (\
                    np.mean(Matrix[triu]),
                    stats.sem(Matrix[triu])),
                    xlim=[-0.5,1])"""

# %%

import matplotlib.pyplot as plt
from scipy.signal import savgol_filter 
from scipy.stats import zscore
folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna', 
                         'PN_cond-NDNF-CB1_WT-vs-KD','spontaneous', 'NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder)


plot_props = dict(with_annotation=True,
                  Ybar=0.5, Ybar_label="0.5$\Delta$F/F",
                  Xbar=0.5, Xbar_label="0.5s",
                  figsize=(9,1.8))



for i, filename in enumerate(DATASET['files']):
    
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    print(i+1, '--', filename, '--', data.nROIs,data.nwbfile.virus)
    print(data.protocols)

    data.build_dFoF(**dFoF_options, verbose=True)
    data.build_running_speed()
    data.build_pupil_diameter()
    

        

    if data.nROIs>0:
        
        """
        ep = physion.analysis.episodes.build.EpisodeData(data, 
                                                        quantities=['dFoF', 'running_speed'],
                                                        protocol_name='imaging during spontaneous behavior')
        
        
        
       
        fig, AX = physion.dataviz.episodes.trial_average.plot(ep,
                                                                quantity='dFoF', with_std=False, 
                                                                color=color,
                                                                roiIndices='all',
                                                                **plot_props)
        
        fig, AX = physion.dataviz.episodes.trial_average.plot(ep,
                                                                quantity='running_speed', with_std=False, 
                                                                color=color,
                                                                
                                                                **plot_props)
        """
        if 'sgRosa' in data.nwbfile.subject.genotype:
            color = 'grey'
            key = 'sgRosa'
            genotype='WT'
        elif 'sgCnr1':
            color = 'darkred'
            key = 'sgCnr1'
            genotype='KD'

    figurepath = '/Users/macbookair/work/Figures/spontaneous/'
    # build figure with traces aligned to dfof time stamps 
    time_dfof =[]
    for i in range(len(data.t_dFoF)):
        time_dfof.append(data.t_dFoF[i]-data.t_dFoF[0])
    
    figurename= data.filename+'_'+genotype+'.svg'
    fig,axs = plt.subplots(3,1, figsize= (12, 8),sharex=True)
    Y = savgol_filter(np.mean(data.dFoF,axis=0), 50, 2)
    Y2 = savgol_filter(data.running_speed, 100, 3)

    pupil_zscore = zscore(data.pupil_diameter)
    Y3 = savgol_filter(pupil_zscore, 100, 3)
    axs[0].set_ylim(bottom=-1.25, top=4)
    axs[0].plot(np.array(data.t_pupil),Y3)
    axs[0].set_ylabel('pupil z-score')
    
    """axs[1].set_ylim(bottom=-0.25, top=4)
    axs[1].plot(data.t_running_speed[300:(300+len(data.t_dFoF))], Y2[300:(300+len(data.t_dFoF))])
    axs[1].set_ylabel('running speed [cm/s]')"""

    axs[1].set_ylim(bottom=-0.25, top=10)
    axs[1].plot(data.t_running_speed, Y2)
    axs[1].set_ylabel('running speed [cm/s]')

    time_dfof =[]
    for i in range(len(data.t_dFoF)):
        time_dfof.append(data.t_dFoF[i]-data.t_dFoF[0])

    axs[2].set_ylim(bottom=0, top=5)
    axs[2].plot(data.t_dFoF, Y, color=color)
    axs[2].set_ylabel('$\\Delta$F/F')
    axs[2].set_xlabel('Time [s]')
    
    #plt.savefig(os.path.join(figurepath+figurename), format='svg')
    plt.show()
    
        
    Matrix = np.corrcoef(data.dFoF)


    # filter the matrix for unique pairs
    triu = np.array(np.triu(np.ones(Matrix.shape), k=1),
                    dtype=bool)
    
    Npairs = np.sum(triu.flatten()) # number of unique pairs

    bar_legend_args={'label':'corr. coeff.',
                                        'ticks':[-1,0,1],
                                        'ticks_labels':['-1', '0', '1']}
    Matrix = np.triu(Matrix, k=1)
    pt.matrix(Matrix, vmin=-1, vmax=1, colormap='seismic',bar_legend_args=bar_legend_args)



    fig, ax = pt.figure()
    ax.hist(Matrix[triu])
    pt.set_plot(ax, xlabel='corr. coef.', ylabel='# pairs',
                title= '%s %.2f $\pm$ %.2f ' % (\
                        genotype,
                        np.mean(Matrix[triu]),
                        stats.sem(Matrix[triu])),
                        
                        xlim=[-0.5,1])

# %%
