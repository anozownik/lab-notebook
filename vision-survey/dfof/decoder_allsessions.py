# %% [markdown]
# # NN decoder for neural activity patterns
# this notebook is copied from yzerlaut Decoding-Responses notebook !

# %%
import os, sys
import numpy as np
from sklearn import linear_model, model_selection, svm
#sys.path += ['/../../../physion/src'] # add src code directory for physion
sys.path.append('../../physion/src')
sys.path.append('../../utils') 

import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
import physion
import physion.utils.plot_tools as pt
from physion.analysis.read_NWB import Data
from physion.analysis.episodes.build import EpisodeData
pt.set_style('ticks')
import datetime
from sklearn.metrics import accuracy_score
import params

# %% [markdown]
# ## Load data


# %%

def build_X_y(ep, averaging_window=[2.5, 3.5], evokedStats=None):

    if evokedStats is None:
        dFoF = ep.dFoF
    else :
        dFoF = ep.dFoF[:, evokedStats['significant'], :]

    averaging_window_cond = (ep.t> averaging_window[0])& (ep.t<averaging_window[1])
    X = np.zeros((dFoF.shape[0], dFoF.shape[1])) # list of matrice responses (Nrois, Ntrials) 

    y = np.zeros(dFoF.shape[0], dtype=int) # label of all trials 

    i=0
    for id in np.unique(getattr(ep, 'Image-ID')):
        pattern_cond = (getattr(ep, 'Image-ID')==id)
        X[i:i+np.sum(pattern_cond), :] =\
            dFoF[pattern_cond,:,:][:,:,averaging_window_cond].mean(axis=2)
        y[i:i+np.sum(pattern_cond)] = id
        i+=np.sum(pattern_cond)

    return X, y

def average_per_class(X, y):
    X_avg = []
    y_avg = []
    for label in np.unique(y):
        X_avg.append(X[y == label].mean(axis=0))
        y_avg.append(label)
    return np.array(X_avg), np.array(y_avg)


# %%
nSeed = 20
denoising = True
averaging_window = [2, 3] # seconds, interval to average to get single activation level per neuron


model = svm.SVC(kernel='linear') #NearestCentroid() #KNeighborsClassifier(n_neighbors=2)
model = KNeighborsClassifier(n_neighbors=1)
folder = os.path.join(os.path.expanduser('~'),
                      'DATA', 'Adrianna', 'PN_cond-NDNF-CB1_WT-vs-KD', 'NWBs')

pname = 'Natural-Images-4-repeats'
DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol=pname)


# %%
DATASET["virus"], DATASET["accuracies"], DATASET["nb_neurons"], DATASET["nb_responsive_neurons"] = [], [], [], []
DATASET["alpha"] = []

for filename in DATASET['files']: # loop over files if you want to test multiple sessions
    data = Data(filename)

    #if data.nwbfile.session_start_time.date() > datetime.date(2026, 4, 1):
    if True:
        
        data.build_dFoF(**params.dFoF_options, verbose=False)
        DATASET["alpha"].append(data.neuropil_correction_factor)
        ep= EpisodeData(data, protocol_name=pname,
                        quantities =['dFoF'], prestim_duration=0)
        
        evokedStats = ep.pre_post_statistics(\
                                params.stat_test_props,
                                response_args=params.response_args,
                                response_significance_threshold=params.response_significance_threshold,
                                loop_over_cells=True,
                                repetition_keys=['Image-ID', 'repeat'],
                                verbose=False
                                ) 
        DATASET["nb_responsive_neurons"].append(np.sum(evokedStats['significant']))
        
        DATASET["nb_neurons"].append(ep.dFoF.shape[1])

        X, y = build_X_y(ep, averaging_window)
        #X, y = build_X_y(ep, averaging_window, evokedStats)
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

    else :
        DATASET["accuracies"].append(np.nan)
        DATASET["virus"].append(data.nwbfile.virus)

# %% Accuracies per session

DATASET["accuracies"] = np.array(DATASET["accuracies"])
chance= 1./len(np.unique(y))
cond_wt = np.array(DATASET["virus"]) == 'syn-Gcamp6s + sgRosa'
cond_kd = np.array(DATASET["virus"]) == 'syn-Gcamp6s + sgCnr1'

fig, ax = pt.figure(ax_scale=(0.8,1.))
ax.bar([0], [np.nanmean(DATASET["accuracies"][cond_wt])], yerr=[np.nanstd(DATASET["accuracies"][cond_wt])], label='wt,score=%.2f' % np.nanmean(DATASET["accuracies"][cond_wt]))
ax.bar([1], [np.nanmean(DATASET["accuracies"][cond_kd])], yerr=[np.nanstd(DATASET["accuracies"][cond_kd])], label='kd, score=%.2f' % np.nanmean(DATASET["accuracies"][cond_kd]))
ax.scatter(np.zeros(sum(cond_wt)), DATASET["accuracies"][cond_wt], color='grey')
ax.scatter(np.ones(sum(cond_kd)), DATASET["accuracies"][cond_kd], color='grey')
ax.plot([-1, 1], [chance, chance], ':', label='chance level, %.2f' % chance)
ax.legend(loc=(1.,0.0), frameon=False)
pt.set_plot(ax, ylabel='accuracy', xticks_labels=[])

# %% Accuracies averaged per subject
DATASET["accuracies"] = np.array(DATASET["accuracies"])
chance= 1./len(np.unique(y))

data = { key : [] for key in ["subjects", "virus", "accuracies"]}

for subject in np.unique(DATASET["subjects"]):
    subject_cond = np.array(DATASET["subjects"]) == subject
    data["subjects"].append(subject)
    data["virus"].append(np.array(DATASET["virus"])[subject_cond][0])
    data["accuracies"].append(np.mean(DATASET["accuracies"][subject_cond]))

data["accuracies"] = np.array(data["accuracies"])
cond_wt = np.array(data["virus"]) == 'syn-Gcamp6s + sgRosa'
cond_kd = np.array(data["virus"]) == 'syn-Gcamp6s + sgCnr1'
    

fig, ax = pt.figure(ax_scale=(0.8,1.))
ax.bar([0], [np.nanmean(data["accuracies"][cond_wt])], yerr=[np.nanstd(data["accuracies"][cond_wt])], label='wt,score=%.2f' % np.nanmean(data["accuracies"][cond_wt]))
ax.bar([1], [np.nanmean(data["accuracies"][cond_kd])], yerr=[np.nanstd(data["accuracies"][cond_kd])], label='kd, score=%.2f' % np.nanmean(data["accuracies"][cond_kd]))
ax.scatter(np.zeros(sum(cond_wt)), data["accuracies"][cond_wt], color='grey')
ax.scatter(np.ones(sum(cond_kd)), data["accuracies"][cond_kd], color='grey')
ax.plot([-1, 1], [chance, chance], ':', label='chance level, %.2f' % chance)
ax.legend(loc=(1.,0.0), frameon=False)
pt.set_plot(ax, ylabel='accuracy', xticks_labels=[])

#%%
DATASET['session'] = np.array([os.path.basename(file)[:-4] for file in DATASET['files']])
df = pd.DataFrame(DATASET).drop(columns=['protocol_ids', 'protocols', 'files', 'dates', 'ages']).set_index('session')
df['perc_resp'] = df['nb_responsive_neurons'] / df['nb_neurons']
df = df.drop(columns=['nb_neurons'])

(df[df['virus'] == 'syn-Gcamp6s + sgRosa']).drop(columns=['virus']).sort_values(by=['subjects', 'session'])

#(df[df['virus'] == 'syn-Gcamp6s + sgCnr1']).drop(columns=['virus']).sort_values(by=['subjects', 'session'])

#%%
DATASET['session'] = np.array([os.path.basename(file)[:-4] for file in DATASET['files']])
if 'nb_responsive_neurons' in DATASET.keys():
    del DATASET['nb_responsive_neurons']

df = pd.DataFrame(DATASET).drop(columns=['protocol_ids', 'protocols','files', 'dates']).set_index('session')

df.sort_values(by=['virus', 'subjects', 'session'])

path = r"Y:\raw-imaging\Adrianna\experiments\analysis\decoder_results_final_sessions.xlsx"
#df.to_excel(path, index=True)

#%%

from sklearn.linear_model import LinearRegression
viruses = np.unique(DATASET['virus'])
color_virus = {'syn-Gcamp6s + sgRosa' : 'grey', 
               'syn-Gcamp6s + sgCnr1': "darkred"}

DATASET['nb_neurons'] = np.array(DATASET['nb_neurons'])
DATASET['virus'] = np.array(DATASET['virus'])

model = LinearRegression()
model.fit(DATASET['nb_neurons'].reshape(-1, 1), DATASET['accuracies'])

model_wt = LinearRegression()
model_wt.fit(DATASET['nb_neurons'][DATASET['virus'] == 'syn-Gcamp6s + sgRosa'].reshape(-1,1), DATASET['accuracies'][DATASET['virus'] == 'syn-Gcamp6s + sgRosa'])

model_kt = LinearRegression()
model_kt.fit(DATASET['nb_neurons'][DATASET['virus'] == 'syn-Gcamp6s + sgCnr1'].reshape(-1,1), DATASET['accuracies'][DATASET['virus'] == 'syn-Gcamp6s + sgCnr1'])

fig, ax = pt.figure(figsize=(3.,3.))

for i, virus in enumerate(viruses):
    ax.scatter(DATASET['nb_neurons'][DATASET['virus'] == virus], DATASET['accuracies'][DATASET['virus'] == virus], 
               color=color_virus[virus], alpha=0.5, label=virus)
pt.set_plot(ax, xlabel='Nb neurons', ylabel='accuracy')

x = np.linspace(np.min(DATASET['nb_neurons']), np.max(DATASET['nb_neurons']), 100)
y =  model.coef_ * x + model.intercept_
ax.plot(x, y, color='black')
pt.annotate(ax, 'y=%.3f *x + %.2f' % (model.coef_, model.intercept_), xy=(0.7, 1.1), ha='center')

""" x = np.linspace(np.min(DATASET['nb_neurons']), np.max(DATASET['nb_neurons']), 100)

y =  model_wt.coef_ * x + model_wt.intercept_
ax.plot(x, y, color=color_virus['syn-Gcamp6s + sgRosa'])
y =  model_kt.coef_ * x + model_kt.intercept_
ax.plot(x, y, color=color_virus['syn-Gcamp6s + sgCnr1'])
pt.annotate(ax, 'y_wt=%.3f *x + %.2f\ny_kt=%.3f *x + %.2f' % (model_wt.coef_, model_wt.intercept_, model_kt.coef_, model_kt.intercept_), 
            xy=(0.7, 1.1), ha='center', fontsize=5) """

#%%

from sklearn.linear_model import LinearRegression
viruses = np.unique(DATASET['virus'])
color_virus = {'syn-Gcamp6s + sgRosa' : 'grey', 
               'syn-Gcamp6s + sgCnr1': "darkred"}

DATASET['alpha'] = np.array(DATASET['alpha'])
DATASET['virus'] = np.array(DATASET['virus'])

model = LinearRegression()
model.fit(DATASET['alpha'].reshape(-1, 1), DATASET['accuracies'])

model_wt = LinearRegression()
model_wt.fit(DATASET['alpha'][DATASET['virus'] == 'syn-Gcamp6s + sgRosa'].reshape(-1,1), DATASET['accuracies'][DATASET['virus'] == 'syn-Gcamp6s + sgRosa'])

model_kt = LinearRegression()
model_kt.fit(DATASET['alpha'][DATASET['virus'] == 'syn-Gcamp6s + sgCnr1'].reshape(-1,1), DATASET['accuracies'][DATASET['virus'] == 'syn-Gcamp6s + sgCnr1'])

fig, ax = pt.figure(figsize=(3.,3.))

for i, virus in enumerate(viruses):
    ax.scatter(DATASET['alpha'][DATASET['virus'] == virus], DATASET['accuracies'][DATASET['virus'] == virus], 
               color=color_virus[virus], alpha=0.5, label=virus)
pt.set_plot(ax, xlabel='alphas', ylabel='accuracy')

""" x = np.linspace(np.min(DATASET['alpha']), np.max(DATASET['alpha']), 100)
y =  model.coef_ * x + model.intercept_
ax.plot(x, y, color='black')
pt.annotate(ax, 'y=%.3f *x + %.2f' % (model.coef_, model.intercept_), xy=(0.7, 1.1), ha='center') """


x = np.linspace(np.min(DATASET['alpha']), np.max(DATASET['alpha']), 100)

y =  model_wt.coef_ * x + model_wt.intercept_
ax.plot(x, y, color=color_virus['syn-Gcamp6s + sgRosa'])
y =  model_kt.coef_ * x + model_kt.intercept_
ax.plot(x, y, color=color_virus['syn-Gcamp6s + sgCnr1'])
pt.annotate(ax, 'y_wt=%.3f *x + %.2f\ny_kt=%.3f *x + %.2f' % (model_wt.coef_, model_wt.intercept_, model_kt.coef_, model_kt.intercept_), 
            xy=(0.7, 1.1), ha='center', fontsize=5)

#%%
from sklearn.linear_model import LinearRegression
viruses = np.unique(DATASET['virus'])
color_virus = {'syn-Gcamp6s + sgRosa' : 'grey', 
               'syn-Gcamp6s + sgCnr1': "darkred"}

DATASET['alpha'] = np.array(DATASET['alpha'])
DATASET['virus'] = np.array(DATASET['virus'])
DATASET['recording_order'] = np.arange(len(DATASET['virus']))

model = LinearRegression()
model.fit(DATASET['recording_order'].reshape(-1, 1), DATASET['alpha'])

model_wt = LinearRegression()
model_wt.fit(DATASET['recording_order'][DATASET['virus'] == 'syn-Gcamp6s + sgRosa'].reshape(-1,1), DATASET['alpha'][DATASET['virus'] == 'syn-Gcamp6s + sgRosa'])

model_kt = LinearRegression()
model_kt.fit(DATASET['recording_order'][DATASET['virus'] == 'syn-Gcamp6s + sgCnr1'].reshape(-1,1), DATASET['alpha'][DATASET['virus'] == 'syn-Gcamp6s + sgCnr1'])

fig, ax = pt.figure(figsize=(3.,3.))

for i, virus in enumerate(viruses):
    ax.scatter(DATASET['recording_order'][DATASET['virus'] == virus], DATASET['alpha'][DATASET['virus'] == virus], 
               color=color_virus[virus], alpha=0.5, label=virus)
pt.set_plot(ax, xlabel='recording order', ylabel='alphas')

""" x = np.linspace(np.min(DATASET['recording_order']), np.max(DATASET['recording_order']), 100)
y =  model.coef_ * x + model.intercept_
ax.plot(x, y, color='black')
pt.annotate(ax, 'y=%.3f *x + %.2f' % (model.coef_, model.intercept_), xy=(0.7, 1.1), ha='center') """


x = np.linspace(np.min(DATASET['recording_order']), np.max(DATASET['recording_order']), 100)

y =  model_wt.coef_ * x + model_wt.intercept_
ax.plot(x, y, color=color_virus['syn-Gcamp6s + sgRosa'])
y =  model_kt.coef_ * x + model_kt.intercept_
ax.plot(x, y, color=color_virus['syn-Gcamp6s + sgCnr1'])
pt.annotate(ax, 'y_wt=%.3f *x + %.2f\ny_kt=%.3f *x + %.2f' % (model_wt.coef_, model_wt.intercept_, model_kt.coef_, model_kt.intercept_), 
            xy=(0.7, 1.1), ha='center', fontsize=5)