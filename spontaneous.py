import numpy as np
import os, sys
sys.path.append('./physion/src') # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')
from scipy import stats




#%%

filename = os.path.join(os.path.expanduser('~'),  'DATA', 'Adrianna',
                        'spontaneous', 'NWBs',
                        '2025_04_02-15-11-39.nwb'
                        )
 

data = physion.analysis.read_NWB.Data(filename,
                                     verbose=False)


# %%
from physion.imaging.Calcium import compute_F0


dFoF_options = dict(\
    method_for_F0='sliding_minmax',
    #method_for_F0='sliding_percentile',
    #method_for_F0='percentile', #more strict
    sliding_window= 300,
    percentile=10,
    roi_to_neuropil_fluo_inclusion_factor=1.,
    neuropil_correction_factor=0.7, 
    with_computed_neuropil_fact=False,
    with_correctedFluo_and_F0=False)

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
# %%
