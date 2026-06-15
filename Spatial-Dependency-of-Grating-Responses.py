# %%
import numpy as np
import matplotlib.pyplot as plt
import os, sys
sys.path += ['./physion/src']

import physion.utils.plot_tools as pt
from physion.analysis.read_NWB import Data,\
      scan_folder_for_NWBfiles
from physion.analysis.episodes.build import EpisodeData
from physion.dataviz.imaging import show_CaImaging_FOV

pt.set_style('dark')
from scipy import stats

datafolder = os.path.expanduser(\
                '~/DATA/Adrianna/NDNF/NWBs')

# %%
def compute(filename):
    data = Data(filename)
    data.build_dFoF()
    ep = EpisodeData(data, 
                 quantities = ['dFoF'],
                 protocol_name='drifting-grating')
    significant = {'positive':[], 'negative':[]}
    responses = {'positive':[], 'negative':[]}
    for i in range(ep.dFoF.shape[1]):

        for sign in ['positive', 'negative']:
            stat = ep.stat_test_for_evoked_responses(
                response_args=dict(roiIndex=i),
                interval_pre=[-1,0],
                interval_post=[1,2],
                sign=sign
            )
            if stat.significant(threshold=0.05):
                # fig, ax = pt.figure()
                # ax.plot(ep.t, ep.dFoF[:,i,:].mean(axis=0))
                significant[sign].append(i)
                responses[sign].append(\
                        ep.dFoF[:,i,:].mean(axis=0))
    return significant, responses, data, ep
        


# %%

filename = os.path.join(
    datafolder, '2025_07_31-17-44-26.nwb')

def analyze(filename):

    significant, responses, data, ep = compute(filename)

    fig, AX = pt.figure((4,1), ax_scale=(1,1.5), top=1.5)
    for ax1, ax2, sign in zip(AX[1::2], AX[::2],
                              ['positive', 'negative']):

        for i, resp in zip(significant[sign], responses[sign]):
            ax1.plot(ep.t, resp-resp[ep.t<0].mean())
        show_CaImaging_FOV(data, NL=4,
                        with_annotation=False,
                        with_ROI_annotation=True,
                        roiIndex=significant[sign],
                        ax=ax2)
        ax2.set_title('%i %s ROIs' %\
            (len(significant[sign]), sign))
    ratio = len(significant['negative'])/\
        (len(significant['negative'])+len(significant['positive']))
    fig.suptitle( '%s, %s : ratio = %.1f ' % (\
        data.metadata['subject_ID'],
        os.path.basename(filename), ratio))
    return ratio

analyze(filename)

# %%
DATASET = scan_folder_for_NWBfiles(datafolder)
results = {'file':[], 'ratio':[]}
for f in DATASET['files']:

    try:
        ratio = analyze(f)
        results['ratio'].append(ratio)
        results['file'].append(os.path.basename(f))
    except BaseException as be:
        pass

# %%

fig, ax = pt.figure()
ax.hist(results['ratio'])
pt.set_plot(ax, xlabel='neg.-to-pos. resp.\nratio',
            ylabel='count')
