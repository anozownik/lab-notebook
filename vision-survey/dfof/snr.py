# %%
import numpy as np
import os, sys, tempfile

sys.path.append('../../physion/src') # add src code directory for physion
sys.path.append('../../utils') 

import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')

from physion.analysis.episodes.build import EpisodeData
from physion.dataviz.episodes.trial_average import plot

import plots_fcts as pt_fcts
import params
import tools

from scipy import stats
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import sem
from itertools import product
import datetime
from sklearn.metrics import auc

# %%
neuron = 'PN_cond-NDNF-CB1_WT-vs-KD'  # 'PN_cond-NDNF-CB1_WT-vs-KD' or 'NDNF-cond-CB1'

folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna',
                        neuron, 'NWBs')

nwb = "2025_09_17-16-04-30.nwb"
filename = os.path.join(folder, nwb)

data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)

print(data.protocols)

dFoF_options = dict(\
    method_for_F0='sliding_minmax',
    sliding_window=300,
    roi_to_neuropil_fluo_inclusion_factor=1.6,
    roi_to_neuropil_fluo_inclusion_factor_metric='std',
    neuropil_correction_factor=0.7, 
    with_computed_neuropil_fact=True)

data.build_dFoF(**dFoF_options, verbose=False)

#%%

pnames = [pname for pname in data.protocols if ('grey' in pname) or ('black' in pname)]
print(pnames)
snr_stim = {}

for pname in pnames:
    ep = EpisodeData(data, quantities=['dFoF'], 
                    protocol_name=pname,
                    with_visual_stim=True)

    print(ep.varied_parameters)

    sigma_noise = np.std(ep.dFoF, axis=(0,2))

    snr_stim.update({pname : ( np.percentile(ep.dFoF, 95, axis=(0,2)) - np.percentile(ep.dFoF, 5, axis=(0,2)) ) / sigma_noise})
    
    snr = snr_stim[pname]
    suite2p_rois = data.valid_roiIndices_suite2p[np.where(snr < 2.5)]
    fig, ax = plt.subplots(len(np.where(snr < 2.5)[0]), 1, figsize=(5, 5), sharex=True)
    for i, roi in enumerate(np.where(snr < 2.5)[0]):
        pt.plot(ep.t, ep.dFoF[0, roi, :], ax=ax[i])
        ax[i].set_ylabel(f"ROI: {suite2p_rois[i]}")

    fig.suptitle(f"{nwb}, {pname}")

#%%
for key in snr_stim.keys():
    
    snr = snr_stim[key]

    fig, AX = plt.subplots(1, 1, figsize=(2, 2))
    hist = AX.hist(snr)
    AX.vlines(2.5, 0, hist[0].max(), color='r', linestyles='--')
    AX.set_xlabel('SNR')
    AX.set_title(f"{nwb}, {key}")

    print(data.valid_roiIndices_suite2p[np.where(snr < 2.5)])

#%%
snr = np.zeros(snr_stim[pnames[0]].shape[0])
for i, key in enumerate(snr_stim.keys()):
    snr += snr_stim[key]
snr /= len(snr_stim.keys())

fig, AX = plt.subplots(1, 1, figsize=(2, 2))
hist = AX.hist(snr)
AX.vlines(2.5, 0, hist[0].max(), color='r', linestyles='--')
AX.set_xlabel('SNR')
AX.set_title(f"{nwb}")

print(data.valid_roiIndices_suite2p[np.where(snr < 2.5)])

#%%
suite2p_rois = data.valid_roiIndices_suite2p[np.where(snr < 2.5)]

for pname in pnames:
    ep = EpisodeData(data, quantities=['dFoF'], 
                    protocol_name=pname,
                    with_visual_stim=True)
    
    fig, ax = plt.subplots(len(np.where(snr < 2.5)[0]), 1, figsize=(5, 5), sharex=True)
    for i, roi in enumerate(np.where(snr < 2.5)[0]):
        pt.plot(ep.t, ep.dFoF[0, roi, :], ax=ax[i])
        ax[i].set_ylabel(f"ROI: {suite2p_rois[i]}")

    fig.suptitle(f"{nwb}, {pname}")
