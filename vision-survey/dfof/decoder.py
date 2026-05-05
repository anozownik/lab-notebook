# %% [markdown]
# # NN decoder for neural activity patterns
# this notebook is copied from yzerlaut Decoding-Responses notebook !

# %%
import os, sys
import numpy as np
from sklearn import linear_model, model_selection, svm
#sys.path += ['/../../../physion/src'] # add src code directory for physion
sys.path.append('../../physion/src')
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
import physion
import physion.utils.plot_tools as pt
from physion.analysis.read_NWB import Data
from physion.analysis.episodes.build import EpisodeData
pt.set_style('ticks')


# %% [markdown]
# ## Load data

# %%
# load a datafile
filename = os.path.join(os.path.expanduser('~'),
                        'DATA', 'Adrianna', 'PN_cond-NDNF-CB1_WT-vs-KD','NWBs',
                        '2025_09_17-15-14-34.nwb')
                        
# example for wt: late'2025_09_17-15-14-34.nwb' mouse 056; wt early 2026_03_12-17-33-32.nwb mouse 006sG; '2026_03_19-16-34-03.nwb' 5 weeks mouse 006sF, 2026_03_09-15-22-59.nwb 4 weeks mouse 006sF
# example for kd: late '2025_09_16-11-10-39.nwb' mouse 075; 8 weeks '2025_09_17-17-56-58.nwb' mouse 1, 3 weeks "2026_03_12-16-39-07.nwb" mouse 004sE, 8 weeks 2026_04_09-14-47-28.nwb mouse 004sE

from physion.imaging.Calcium import compute_F0


dFoF_options = dict(\
    method_for_F0='sliding_minmax',
    #method_for_F0='sliding_percentile',
    #method_for_F0='percentile', #more strict
    sliding_window= 300,
    percentile=10,
    roi_to_neuropil_fluo_inclusion_factor=1.1,
    neuropil_correction_factor=0.7, 
    with_computed_neuropil_fact=False)

data = Data(filename)
data.build_dFoF(**dFoF_options, verbose=False)
ep= EpisodeData(data, protocol_name='Natural-Images-4-repeats',
                quantities =['dFoF'])


# %% [markdown]

# ## 


# --------------------------------------------- #
# Transforming to matrix and list  for sklearn X, y respectively
# --------------------------------------------- #

def build_X_y(ep, averaging_window=[0, 2]):

    averaging_window_cond = (ep.t> averaging_window[0])& (ep.t<averaging_window[1])
    X = np.zeros((ep.dFoF.shape[0], ep.dFoF.shape[1])) # list of matrice responses (Nrois, Ntrials) 

    y = np.zeros(ep.dFoF.shape[0], dtype=int) # label of all trials 

    i=0
    for id in np.unique(getattr(ep, 'Image-ID')):
        pattern_cond = (getattr(ep, 'Image-ID')==id)
        X[i:i+np.sum(pattern_cond), :] =\
            ep.dFoF[pattern_cond,:,:][:,:,averaging_window_cond].mean(axis=2)
        y[i:i+np.sum(pattern_cond)] = id
        i+=np.sum(pattern_cond)

    return X, y
X, y = build_X_y(ep, averaging_window_cond)
#%%
# normalization of input data
normed=True
if normed:
     from sklearn.preprocessing import StandardScaler
     scaler = StandardScaler()
     X = scaler.fit_transform(X,y)

#%%

# %%
fig, AX = pt.figure(axes=(1, len(np.unique(y))), ax_scale=(2,.6))
for id, ax in zip(np.unique(y), AX):

        ax.bar(range(data.nROIs),
        X[y==id].mean(axis=0),
        yerr = X[y==id].std(axis=0),
        color=None)

        pt.set_plot(ax,
        ylabel='$\\Delta$F/F',
        xlabel='' if ax!=AX[-1] else 'ROIs',
        xticks_labels=[] if ax!=AX[-1] else None)
        pt.annotate(ax, 'image #%i' % id, (0,1), va='top')
        pt.annotate(ax, '(%i trials)' % np.sum(y==id),
        (1,1), va='top', ha='right', fontsize='small')
pt.set_common_ylims(AX)
# %%
# split data into training and test sets - stratified strategy to have same amount of a given image in the train and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                test_size=0.5,
                                                random_state=42,
                                                stratify=y)

# check that the split is balanced
pd.value_counts(y_test)

# %%
# possibility to denoise training set by averaging
USE_DENOISING = False

def average_per_class(X, y):
    X_avg = []
    y_avg = []
    for label in np.unique(y):
        X_avg.append(X[y == label].mean(axis=0))
        y_avg.append(label)
    return np.array(X_avg), np.array(y_avg)


if USE_DENOISING:
    X_train, y_train = average_per_class(X_train, y_train)

"""if denoising:
        X_train = pd.DataFrame(
               np.array([X_train[y_train==id].mean(axis=0)\
                        for id in np.unique(y_train)]),\
                        columns=['ROI %i' % i for i in range(data.nROIs)])
        y_train = pd.DataFrame(\
                 np.array([id for id in np.unique(y_train)]),\
                 columns=['Image-ID'])  """      
# %%
fig, AX = pt.figure(axes=(1, len(np.unique(y)) + 1), 
                    ax_scale=(2, .6))

# --- Test set (single trials) ---
for x in X_test:
    AX[0].plot(range(data.nROIs), x, lw=0.1)

pt.annotate(AX[0],
            f'Test set: ({len(X_test)} single trials)',
            (0, 1), va='top')

# plot mean response across ALL test data
AX[0].plot(range(data.nROIs),
           X_test.mean(axis=0),
           color='k', lw=2)


# --- Training set (averaged per class) ---
for i, (label, ax) in enumerate(zip(np.unique(y), AX[1:])):

    class_data = X_train[y_train == label]

    mean_response = class_data.mean(axis=0)
    std_response = class_data.std(axis=0)

    ax.bar(range(data.nROIs),
           mean_response,
           yerr=std_response)

    pt.annotate(ax,
                f'Training Set: image #{label}',
                (0, 1), va='top')

    pt.annotate(ax,
                f'({len(class_data)} trials)',
                (1, 1), va='top', ha='right', fontsize='small')


# --- formatting ---
for ax in AX:
    pt.set_plot(ax,
                ylabel='$\\Delta$F/F',
                xlabel='' if ax != AX[-1] else 'ROIs',
                xticks_labels=[] if ax != AX[-1] else None)

"""fig, AX = pt.figure(axes=(1, len(np.unique(y))+1), 
                    ax_scale=(2,.6))

for x in np.array(X_test):
    AX[0].plot(range(data.nROIs), x, lw=0.1, color=None)
pt.annotate(AX[0], 'Test set: (%i single trials)' % len(X_test), 
            (0,1), va='top')
AX[0].plot(range(data.nROIs), X[y==id].mean(axis=0), color=None)



for id, ax in zip(np.unique(y), AX[1:]):

            
                ax.bar(range(data.nROIs), 
                        X_train.iloc[id,:],
                        yerr=X_train.iloc[id,:].std(axis=0), 
                        color=None)
                pt.annotate(ax, 'Training Set: image #%i' % id, (0,1), va='top')
                pt.annotate(ax, '(%i samples)' % np.sum(y_train['Image-ID']==id), 
                                (1,1), va='top', ha='right', fontsize='small')

for ax in AX:
    pt.set_plot(ax, 
                ylabel='$\\Delta$F/F',
                xlabel='' if ax!=AX[-1] else 'ROIs',
                xticks_labels=[] if ax!=AX[-1] else None)"""
# %%


# Decoding with a linear NN classifier


import matplotlib.pyplot as plt
from sklearn.neighbors import NearestCentroid


model = svm.SVC(kernel='linear') #NearestCentroid() #KNeighborsClassifier(n_neighbors=2)


model.fit(X_train, y_train)
y_pred = model.predict(X_test)

from sklearn.metrics import accuracy_score

score=accuracy_score(y_test, y_pred)
chance= 1./len(np.unique(y))

fig, ax = pt.figure(ax_scale=(0.8,1.))
ax.bar([0], [score], label='SVM, score=%.2f' % score)
ax.plot([-1, 1], [chance, chance], ':', label='chance level, %.2f' % chance)
ax.legend(loc=(1.,0.0), frameon=False)
pt.set_plot(ax, ylabel='accuracy', xticks_labels=[])


# %%

# %%
nSeed = 20
denoising = True
averaging_window = [0, 2] # seconds, interval to average to get single activation level per neuron

model = svm.SVC(kernel='linear') #NearestCentroid() #KNeighborsClassifier(n_neighbors=2)
model = KNeighborsClassifier(n_neighbors=1)
folder = os.path.join(os.path.expanduser('~'),
                      'DATA', 'Adrianna', 'PN_cond-NDNF-CB1_WT-vs-KD', '20260325','PNs','NWBs', '2026_march_and_april')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol='Natural-Images-4-repeats')


# %%
DATASET["virus"], DATASET["accuracies"] = [], []

for filename in DATASET['files']: # loop over files if you want to test multiple sessions
    data = Data(filename)
    data.build_dFoF(**dFoF_options, verbose=False)
    ep= EpisodeData(data, protocol_name='Natural-Images-4-repeats',
                    quantities =['dFoF'], prestim_duration=0)  
    X, y = build_X_y(ep)
    accuracies = []
    for random_state in range(nSeed):
        X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    test_size=0.5,
                                                    random_state=random_state,
                                                    stratify=y)
        if denoising:
            X_train, y_train = average_per_class(X_train, y_train)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        score=accuracy_score(y_test, y_pred)
        accuracies.append(score)
    DATASET["accuracies"].append(np.mean(accuracies))
    DATASET["virus"].append(data.nwbfile.virus)
# %%

DATASET["accuracies"] =np.array(DATASET["accuracies"])

cond_wt = np.array(DATASET["virus"]) == 'syn-Gcamp6s + sgRosa'
cond_kd = np.array(DATASET["virus"]) == 'syn-Gcamp6s + sgCnr1'
    

fig, ax = pt.figure(ax_scale=(0.8,1.))
ax.bar([0], [DATASET["accuracies"][cond_wt].mean()], yerr=[DATASET["accuracies"][cond_wt].std()], label='wt,score=%.2f' % DATASET["accuracies"][cond_wt].mean())
ax.bar([1], [DATASET["accuracies"][cond_kd].mean()], yerr=[DATASET["accuracies"][cond_kd].std()], label='kd, score=%.2f' % DATASET["accuracies"][cond_kd].mean())
ax.scatter(np.zeros(sum(cond_wt)), DATASET["accuracies"][cond_wt], color='grey')
ax.scatter(np.ones(sum(cond_kd)), DATASET["accuracies"][cond_kd], color='grey')
ax.plot([-1, 1], [chance, chance], ':', label='chance level, %.2f' % chance)
ax.legend(loc=(1.,0.0), frameon=False)
pt.set_plot(ax, ylabel='accuracy', xticks_labels=[])







# %%
