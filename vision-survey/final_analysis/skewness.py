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
from scipy.stats import skew

# %%

# TO LOOP OVER NWB FILES WITH VISUAL STIMULUS --- DRIFITING GRATING ---  multisession

folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna',
                        #'NDNF-cond-CB1', 'NWBs')
                        'PN_cond-NDNF-CB1_WT-vs-KD', 'final_NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                                            for_protocol=None)


#%%
DATASET_ROIS = dict(mouse=[], virus=[], session=[], rois_id=[], mean_dFoF=[], skewness=[])

dFoF_options = params.dFoF_options
dFoF_options['with_computed_neuropil_fact'] = False

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

    data.build_dFoF(**dFoF_options, verbose=True)

    if data.nROIs>0:
        
        for roi_idx in range(data.dFoF.shape[0]):
            DATASET_ROIS['mouse'].append(DATASET['subjects'][i])
            DATASET_ROIS['virus'].append(virus)
            DATASET_ROIS['session'].append(os.path.basename(filename)[:-4])
            DATASET_ROIS['rois_id'].append(roi_idx)
            DATASET_ROIS['skewness'].append(skew(data.dFoF[roi_idx]))
            DATASET_ROIS['mean_dFoF'].append(np.mean(data.dFoF[roi_idx]))
    
    print('')

#%%
# BUILD DATAFRAME ROIS

if 'PN' in folder:
    savepath_excel = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\PN\excels"
    os.makedirs(savepath_excel, exist_ok=True)
    excel_filename = 'summary_data_rois_' + 'PN' + '.xlsx'
else :
    savepath_excel = r"Y:\raw-imaging\Adrianna\experiments\analysis\Adrianna\NDNF\excels"
    os.makedirs(savepath_excel, exist_ok=True)
    excel_filename = 'summary_data_rois_' + 'NDNF' + '.xlsx'



df_rois = pd.DataFrame(DATASET_ROIS)
df_rois.to_excel(os.path.join(savepath_excel, excel_filename))

df_rois
