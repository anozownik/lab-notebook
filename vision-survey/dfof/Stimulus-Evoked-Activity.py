# %% [markdown]
# # Trial Averaging

# %%
import os, sys
import numpy as np
sys.path += ['../src'] # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')


# %% [markdown]
# ## Load data

# %%
# load a datafile
filename = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna',
                        'PN_cond-NDNF-CB1_WT-vs-KD', '20260325','PNs','NWBs', '2025', '2025_09_16-11-10-39.nwb') #\\iss\bacci\raw-imaging\Adrianna\experiments\PYR\2026_04_14\12-24-38\2026_04_14_12-24-38_output_1\2026_04_14_12-24-38_1_Natural-Images-4-repeats

data = physion.analysis.read_NWB.Data(filename,
                                     verbose=False)

# %% [markdown]
# ## Build episodes (stimulus-aligned)


from physion.imaging.Calcium import compute_F0


dFoF_options = dict(\
    method_for_F0='sliding_minmax',
    #method_for_F0='sliding_percentile',
    #method_for_F0='percentile', #more strict
    sliding_window= 300,
    percentile=10,
    roi_to_neuropil_fluo_inclusion_factor=1.1,
    neuropil_correction_factor=0.7, 
    with_computed_neuropil_fact=True)

data.build_dFoF(**dFoF_options, verbose=True)
# %%

# find protocol of full-field gratings
p_name = [p for p in data.protocols if 'Natural-Images-4-repeats' in p][0]
episodes = physion.analysis.episodes.build.EpisodeData(data, 
                                                       quantities=['dFoF'],
                                                       prestim_duration=1.5,
                                                       protocol_name=p_name)

# %% [markdown]
# ## Plot properties

# %%
# plot over varying angles
plot_props = dict(column_key='Image-ID',
                  with_annotation=True,
                  with_std = False,
                  color='r',
                  with_stat_test=True,
                  figsize=(9,1.8))

# %% [markdown]
# ## Pupil variations

# %%
#fig, AX = physion.dataviz.episodes.trial_average.plot(episodes,
#                                                      quantity='pupil_diameter',
#                                                    #   with_std=False,
#                                                      **plot_props)
# %% [markdown]
# ## Average over all ROIs 

# %%
fig, AX = physion.dataviz.episodes.trial_average.plot(episodes,
                                                      quantity='dFoF',
                                                      roiIndex=range(data.nROIs),
                                                      **plot_props)

# %% [markdown]
# ## Single ROIs 

# %%


for i in range(data.nROIs):
        fig, AX = physion.dataviz.episodes.trial_average.plot(episodes,
                                                      roiIndex=i,
                                                      **plot_props)
        


# %%
fig, AX = physion.dataviz.episodes.trial_average.plot(episodes,
                                                      roiIndex=2,
                                                      **plot_props)
#%%


                


import matplotlib.pyplot as plt

#%%
baselineCond = (episodes.t > -1.9) & (episodes.t < 0)
fig, ax = plt.subplots(episodes.data.nROIs, 5, figsize=(15, 3 * episodes.data.nROIs), sharey=True)

for i in range(episodes.data.nROIs):
    for j, img_id in enumerate([1., 2., 3., 4., 5.]):
        
        image_cond = (getattr(episodes, 'Image-ID') == img_id)
        
        data = episodes.dFoF[image_cond, i, :]
        
        mean_trace = np.mean(data, axis=0)
        bsl_mean = mean_trace[baselineCond].mean()
        bsl_substracted = mean_trace - bsl_mean
        
        ax[i, j].plot(episodes.t, bsl_substracted)
        ax[i, j].set_title(f'ROI {i}, Img {img_id}')

plt.tight_layout()
plt.show()

# %%

# statisical test of evoked responses:
stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='wilcoxon',                                            
                       sign='positive')

summary = episodes.pre_post_statistics(\
                stat_test_props,
                # repetition_keys=['repeat', 'angle', 'contrast'],
                response_args=dict(quantity='dFoF',
                                   roiIndex=0),
                response_significance_threshold=0.01,
                verbose=True)

for key in summary:
    print('- %s : %s' % (key, summary[key]))

# %%
# now LOOPING over cells
summary = episodes.pre_post_statistics(\
                stat_test_props,
                response_args=dict(quantity='dFoF'),
                response_significance_threshold=0.01,
                loop_over_cells=True,
                repetition_keys=['repeat'],
                verbose=False)


# %%
summary = episodes.reliability(
        response_args=dict(quantity='dFoF'),
        stat_test_props=dict(n_samples=500, seed=2),
        loop_over_cells=True,
        verbose=False,
)

# %%
# visualizing the computation of reliability:

from physion.analysis.episodes.trial_statistics \
        import run_reliability_test
summary = run_reliability_test(episodes,
                episodes.find_episode_cond(key=['contrast', 'angle'],
                                           value=[0.5, 0]),
                dict(quantity='dFoF', roiIndex=1),
                stat_test_props=dict(n_samples=100, seed=2),
                response_significance_threshold=0.05,
                return_samples=True)

# %%
fig, [ax0, ax]= pt.figure(axes=(2,1), wspace=3.)
pt.plot(episodes.t, 
        np.mean(summary['real'], axis=0), 
        sy=np.std(summary['real'], axis=0), 
        ax=ax0, color='tab:blue', label='real')
pt.plot(episodes.t, 
        np.mean(summary['shuffled'], axis=0), 
        sy=np.std(summary['shuffled'], axis=0), 
        ax=ax0, color='tab:grey', label='shuffled')
pt.set_plot(ax0, 
            xlabel='time (s)',
            ylabel='$\\Delta$F/F')
ax.hist(summary['null_corr_list'], bins=np.linspace(-1,1,20),
        label='Null correlations', color='tab:grey', density=True)
ax.hist(summary['corr_list'], bins=np.linspace(-1,1,20),
        label='True correlations', color='tab:blue', density=True)
ax.axvline(summary['r'], color='green' if summary['significant'] else 'red', linestyle='--', label='Reliability r=%.2f' % summary['r'])
# ax.axvline(perc_threshold, color='black', linestyle='--', label='%.0fth percentile of null dist=%.2f' %(percentile, perc_threshold))
pt.set_plot(ax, xlabel='Correlation coefficient', ylabel='Count')
ax.set_title('r=%(r).2f, p%(pval).1e' % summary)
ax0.legend(loc=(1.,.2), frameon=False)
ax.legend(loc=(1.,.2), frameon=False)


# %%
