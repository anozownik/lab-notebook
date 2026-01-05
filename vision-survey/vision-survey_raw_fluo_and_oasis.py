import numpy as np
import os, sys
sys.path.append('../physion/src') # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')
from scipy import stats


#%%
# PLOT PROPERTIES --- DRIFTING GRATINGS ---
stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='ttest')

#%%

filename = os.path.join(os.path.expanduser('~'),  'DATA', 'Adrianna',
                        'PN_cond-NDNF-CB1_WT-vs-KD', 'NWBs',
                        '2025_10_06-15-30-53.nwb'
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
    with_computed_neuropil_fact=True,
    with_correctedFluo_and_F0=True)

# we first perform the dFoF determination with the above params
#    (this restrict the available ROIs in the future)
data.build_dFoF(**dFoF_options, verbose=True)
#print(data.correctedFluo)
valid = data.valid_roiIndices
rejected = [i for i in range(data.Fluorescence.data.shape[1]) if (i not in valid)]

data.build_rawFluo()
data.build_pupil_diameter()
data.build_facemotion()
data.build_running_speed()
data.build_neuropil()
data.build_Deconvolved()
#%% 
# to mark virus for later 
if 'sgRosa' in data.nwbfile.virus:
       color = 'grey'

elif 'sgCnr1' in data.nwbfile.virus:
       color = 'darkred'

# BUILDING --- DRIFTING GRATINGS ---- EPISODE 
epGrating = physion.analysis.process_NWB.EpisodeData(data, 
                                                    quantities=['rawFluo', 'correctedFluo','Deconvolved','dFoF', 'neuropil','running_speed', 'pupil_diameter', 'facemotion'],
                                                    protocol_name='drifting-grating')

#%%
plot_props = dict(column_key='contrast',
                  with_annotation=True,
                  with_axis = True,
                  #Ybar=50, Ybar_label=" 50a.u",
                  #Xbar=0.5, Xbar_label="0.5s",
                  figsize=(9,1.8))# %%






stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='ttest')

for i in range(epGrating.data.nROIs):
        
        roiIndex=i,
        fig, AX = physion.dataviz.episodes.trial_average.plot(epGrating, quantity = 'Deconvolved',with_std=False,
                                                        roiIndex=roiIndex,with_stat_test=True, stat_test_props=stat_test_props,
                                                        **plot_props)
        pt.show()
# %%
pt.plot(epGrating.t,epGrating.Deconvolved[0,0,:])
# %%
