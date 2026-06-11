# %%
import numpy as np
import matplotlib.pyplot as plt
import os, sys
sys.path += ['./physion/src']

import physion.utils.plot_tools as pt
from physion.analysis.read_NWB import Data,\
      scan_folder_for_NWBfiles
pt.set_style('dark')
from scipy import stats

datafolder = os.path.expanduser(\
                '~/DATA/Adrianna/NDNF/NWBs')

DATASET = scan_folder_for_NWBfiles(datafolder,
                                   for_protocol='vision-survey')
# %%

filename = os.path.join(
    datafolder, '2025_07_31-17-44-26.nwb')
data = Data(filename)
data.build_dFoF()

# %%
from physion.analysis.episodes.build import EpisodeData
ep = EpisodeData(data, 
                 quantities = ['dFoF'],
                 protocol_name='drifting-grating')
# %%

significant = {'positive':[], 'negative':[]}
responses = {'positive':[], 'negative':[]}
for i in range(ep.dFoF.shape[1]):

    for sign in ['positive', 'negative']:
        stat = ep.stat_test_for_evoked_responses(
            response_args=dict(roiIndex=i),
            interval_pre=[-1,0],
            interval_post=[0.5,1.5],
            sign=sign
        )
        if stat.significant(threshold=0.05):
            # fig, ax = pt.figure()
            # ax.plot(ep.t, ep.dFoF[:,i,:].mean(axis=0))
            significant[sign].append(i)
            responses[sign].append(\
                    ep.dFoF[:,i,:].mean(axis=0))
        


# %%
from physion.dataviz.imaging import show_CaImaging_FOV

for sign in ['positive', 'negative']:

    fig, ax = pt.figure(ax_scale=(1.5,1.5))
    for i, resp in zip(significant[sign], responses[sign]):
        ax.plot(ep.t, resp-resp[ep.t<0].mean(), label='roi #%i' % (i+1))
    ax.legend(loc=(1.,0))
    fig= show_CaImaging_FOV(data, NL=4,
                            with_annotation=False,
                            with_ROI_annotation=True,
                            roiIndex=significant[sign])
    fig[1].set_title(sign+' ROIs')
# %%
