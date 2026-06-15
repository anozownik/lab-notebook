# %%
import numpy as np
import matplotlib.pyplot as plt
import os, sys
sys.path += ['./physion/src']

import physion.utils.plot_tools as pt
from physion.analysis.read_NWB import Data,\
      scan_folder_for_NWBfiles
from physion.analysis.episodes.build import EpisodeData

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
    return significant, responses, data
        


# %%

filename = os.path.join(
    datafolder, '2025_07_31-17-44-26.nwb')


from physion.dataviz.imaging import show_CaImaging_FOV

def analyze(filename):

    significant, responses, data = compute(filename)

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
    score = len(significant['negative'])/\
        (len(significant['negative'])+len(significant['positive']))
    fig.suptitle( '%s, %s : score = %.1f ' % (\
        data.metadata['subject_ID'],
        os.path.basename(filename), score))
    return score
# %%
# datafolder = os.path.expanduser(\
#                 '~/DATA/Adrianna/PN/NWBs')

# %%
DATASET = scan_folder_for_NWBfiles(datafolder)
results = {'file':[], 'score':[]}
for f in DATASET['files']:

    try:
        score = analyze(f)
        results['score'].append(score)
        results['file'].append(os.path.basename(f))
    except BaseException as be:
        pass

# %%

# pt.plt.hist(results['score'])

# %%
for i, f in enumerate(os.listdir(datafolder)[:2]):
    if '.nwb' in f:
        print(i, ')', f)
        data = Data(os.path.join(datafolder, f))
        print(data.df_name)

# %%

# DATASET = scan_folder_for_NWBfiles(datafolder,
#                                    for_protocol='vision-survey')
