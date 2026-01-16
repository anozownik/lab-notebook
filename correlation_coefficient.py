import numpy as np
import os, sys
sys.path.append('./physion/src') # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')
from scipy import stats

#%%

filename = os.path.join(os.path.expanduser('~'),  'DATA', 'Adrianna', 
                        'PN_cond-NDNF-CB1_WT-vs-KD', 'spontaneous', 'NWBs',
                        '2025_07_30-16-35-55.nwb'
                        )
 


data = physion.analysis.read_NWB.Data(filename,
                                     verbose=False)


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
    with_computed_neuropil_fact=True)


#%%

data.build_dFoF(**dFoF_options, verbose=True)

ep = physion.analysis.episodes.build.EpisodeData(data,quantities=['dFoF'], protocol_name ='imaging during spontaneous behavior')

# %%
Matrix = np.corrcoef(data.dFoF)
# %%

# filter the matrix for unique pairs
triu = np.array(np.triu(np.ones(Matrix.shape), k=1),
                dtype=bool)
# 
Npairs = np.sum(triu.flatten()) # number of unique pairs

# %%
Matrix = np.triu(Matrix, k=1)
pt.matrix(Matrix, vmin=0, vmax=1)

# %%

fig, ax = pt.figure()
ax.hist(Matrix[triu])
pt.set_plot(ax, xlabel='corr. coef.', ylabel='# pairs',
            title= '%.2f $\pm$ %.2f ' % (\
                    np.mean(Matrix[triu]),
                    stats.sem(Matrix[triu])),
                    xlim=[-0.5,1])

# %%


folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna', 
                        'PN_cond-NDNF-CB1_WT-vs-KD', 'spontaneous','NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder)


plot_props = dict(with_annotation=True,
                  Ybar=0.5, Ybar_label="0.5$\Delta$F/F",
                  Xbar=0.5, Xbar_label="0.5s",
                  figsize=(9,1.8))



for i, filename in enumerate(DATASET['files']):
    
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    print(i+1, '--', filename, '--', data.nROIs)
    print(data.protocols)

    data.build_dFoF(**dFoF_options, verbose=True)
    #data.build_running_speed()
    
    if 'sgRosa' in data.nwbfile.virus:
        color = 'grey'
        key = 'sgRosa'
    elif 'sgCnr1':
        color = 'darkred'
        key = 'sgCnr1'
        

    if data.nROIs>0:
        """

        ep = physion.analysis.process_NWB.EpisodeData(data, 
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
                                                                
                                                                **plot_props)"""
        
        Matrix = np.corrcoef(data.dFoF)


        # filter the matrix for unique pairs
        triu = np.array(np.triu(np.ones(Matrix.shape), k=1),
                        dtype=bool)
        
        Npairs = np.sum(triu.flatten()) # number of unique pairs


        Matrix = np.triu(Matrix, k=1)
        pt.matrix(Matrix, vmin=0, vmax=1)



        fig, ax = pt.figure()
        ax.hist(Matrix[triu])
        pt.set_plot(ax, xlabel='corr. coef.', ylabel='# pairs',
                    title= '%s %.2f $\pm$ %.2f ' % (\
                            filename,
                            np.mean(Matrix[triu]),
                            stats.sem(Matrix[triu])),
                            
                            xlim=[-0.5,1])
# %%
