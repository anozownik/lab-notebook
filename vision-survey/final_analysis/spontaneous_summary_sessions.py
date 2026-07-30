#%%
import numpy as np
import pandas as pd
import os, sys
sys.path.append('./physion/src') # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')
from scipy import stats

sys.path.append('./physion/src') # add src code directory for physion
sys.path.append('./utils') 

import params

import matplotlib.pyplot as plt
from scipy.signal import savgol_filter 
from scipy.stats import zscore

#%%
folder = os.path.join(os.path.expanduser('~'),  'DATA', 'Adrianna', 
                        #'NDNF-cond-CB1', 'spontaneous', 'NWBs')
                        'PN_cond-NDNF-CB1_WT-vs-KD', 'spontaneous', 'NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol=None)


#%%
plot = False
DATASET['virus'], DATASET['mean_corr_coef'] = [], []

plot_props = dict(with_annotation=True,
                  Ybar=0.5, Ybar_label="0.5$\Delta$F/F",
                  Xbar=0.5, Xbar_label="0.5s",
                  figsize=(9,1.8))


for i, filename in enumerate(DATASET['files']):
    
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    print(i+1, '--', filename, '--', data.nROIs,data.nwbfile.virus)
    print(data.protocols)

    data.build_dFoF(**params.dFoF_options, verbose=True)
    data.build_running_speed()
    data.build_pupil_diameter()

    if data.nROIs>0:
        
        if 'sgRosa' in data.nwbfile.subject.genotype:
            color = 'grey'
            key = 'sgRosa'
            genotype='WT'
        elif 'sgCnr1':
            color = 'darkred'
            key = 'sgCnr1'
            genotype='KD'
    
    DATASET['virus'].append(key)
    
    # Correlation matrix
    corr_matrix = np.corrcoef(data.dFoF) # Pearson product-moment correlation coefficient
    corr_triu = np.triu(corr_matrix, k=1) # upper triangular matrix
    triu_mask = np.array(np.triu(np.ones(corr_matrix.shape), k=1), dtype=bool)

    DATASET['mean_corr_coef'].append(np.mean(corr_triu[triu_mask]))

    if plot :
        
        bar_legend_args={'label':'corr. coeff.',
                        'ticks':[-1,0,1],
                        'ticks_labels':['-1', '0', '1'],
                        }

        pt.matrix(corr_triu, vmin=-1, vmax=1, colormap='seismic', bar_legend_args=bar_legend_args)

        fig, ax = pt.figure()
        ax.hist(corr_triu[triu_mask])
        pt.set_plot(ax, xlabel='corr. coef.', ylabel='# pairs',
                    title= '%s %.2f $\pm$ %.2f ' % (\
                            genotype,
                            np.mean(corr_triu[triu_mask]),
                            stats.sem(corr_triu[triu_mask])),
                            
                            xlim=[-0.5,1])

# %%
# BUILD DATAFRAME

if 'PN' in folder:
    savepath_excel = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\PN\excels"
    os.makedirs(savepath_excel, exist_ok=True)
    excel_filename = 'spontaneous_summary_data_' + 'PN' + '.xlsx'
else :
    savepath_excel = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\NDNF\excels"
    os.makedirs(savepath_excel, exist_ok=True)
    excel_filename = 'spontaneous_summary_data_' + 'NDNF' + '.xlsx'

DATASET['session'] = np.array([os.path.basename(file)[:-4] for file in DATASET['files']])

df = pd.DataFrame(DATASET).drop(columns=['protocol_ids', 'protocols', 'files']).set_index('session')

df.to_excel(os.path.join(savepath_excel, excel_filename))

df