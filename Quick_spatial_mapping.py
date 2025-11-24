# %% [markdown]
# # Quick spatial Mapping
# ### Description quick spatial mapping
# 9 static patches at 9 positions x=(-36,0,36) and y=(-26,0,26)


# %%
# load packages:
import os, sys
sys.path += ['/Users/macbookair/work/lab-notebook/physion/src'] # add src code directory for physion
from physion.analysis.read_NWB import Data, scan_folder_for_NWBfiles
from physion.analysis.process_NWB import EpisodeData
from physion.utils  import plot_tools as pt
import numpy as np

import itertools
from physion.dataviz.episodes.trial_average import plot as plot_trial_average
import matplotlib.pyplot as plt



# %% [markdown]
# ## Load data

# %%
datafolder = os.path.join(os.path.expanduser('~'), 'work', 'DATA', 'NDNF','NWBs')
SESSIONS = scan_folder_for_NWBfiles(datafolder)
SESSIONS['nwbfiles'] = [os.path.basename(f) for f in SESSIONS['files']]

#%%
dFoF_options = {'roi_to_neuropil_fluo_inclusion_factor' : 1.0, # ratio to discard ROIs with weak fluo compared to neuropil
                 'method_for_F0' : 'sliding_percentile', # either 'minimum', 'percentile', 'sliding_minimum', or 'sliding_percentile'
                 'sliding_window' : 300. , # seconds (used only if METHOD= 'sliding_minimum' | 'sliding_percentile')
                 'percentile' : 10. , # for baseline (used only if METHOD= 'percentile' | 'sliding_percentile')
                 'neuropil_correction_factor' : 0.8 }# fraction of neuropil substracted to fluorescence

coord_map = {
            (np.float64(-36.0), np.float64(-23.0)): (2, 0),
            (np.float64(-36.0), np.float64(0.0)):   (2, 1),
            (np.float64(-36.0), np.float64(23.0)):  (2, 2),
            (np.float64(0.0), np.float64(-23.0)):   (1, 0),
            (np.float64(0.0), np.float64(0.0)):     (1, 1),
            (np.float64(0.0), np.float64(23.0)):    (1, 2),
            (np.float64(36.0), np.float64(-23.0)):  (0, 0),
            (np.float64(36.0), np.float64(0.0)):    (0, 1),
            (np.float64(36.0), np.float64(23.0)):   (0, 2),
        }

#%%
filename = SESSIONS['files'][0]
data = Data(filename,
                verbose=False)

data.build_dFoF(**dFoF_options, verbose=False)
data.build_pupil_diameter()
data.build_running_speed()

quantities = ['dFoF']
protocol = "quick-spatial-mapping"
ep = EpisodeData(data, 
                quantities = quantities, 
                protocol_name = protocol, 
                verbose=False)

varied_keys = [k for k in ep.varied_parameters.keys() if k!='repeat']
varied_values = [ep.varied_parameters[k] for k in varied_keys]
print(varied_keys)
print(varied_values)



#%% plot quick spatial mapping
from scipy import stats
def plot_qsm(index, coord_map, diffs):
    filename = SESSIONS['files'][index]
    data = Data(filename,
                verbose=False)

    data.build_dFoF(**dFoF_options, verbose=False)
    data.build_pupil_diameter()
    data.build_running_speed()

    quantities = ['dFoF']
    protocol = "quick-spatial-mapping"
    ep = EpisodeData(data, 
                    quantities = quantities, 
                    protocol_name = protocol, 
                    verbose=False)
    
    fig, AX = plt.subplots(3, 3, figsize = (10,10))  # 3x3 grid

    for values in itertools.product(*varied_values):

        stim_cond = ep.find_episode_cond(key=varied_keys, value=values)
        response = ep.get_response2D(quantity='dFoF',
                                    episode_cond=stim_cond)
        response = ep.dFoF[stim_cond,:,:].mean(axis=1)
        
        i, j = coord_map[values]
        pt.plot(ep.t, response.mean(axis=0), sy=stats.sem(response,axis=0), ax=AX[i][j])
        AX[i][j].axvspan(0, 1, color="grey", alpha=0.3)
        

    #quantification to rank 
    summary = ep.compute_summary_data(stat_test_props={})
    center_cond = (summary['x-center']==0) & (summary['y-center']==0)
    center_value = summary['value'][center_cond]
    center_resp = summary['significant'][center_cond]
    fig.text(0.4, 0.26, f"Center responsive ={center_resp}", ha="center", fontsize=12)
    around_value = summary['value'][~center_cond].mean()
    diff_value = center_value-around_value
    fig.text(0.4, 0.24, f"Center ={center_value[0]:.3f}", ha="center", fontsize=12)
    fig.text(0.4, 0.22, f"Around ={around_value:.3f}", ha="center", fontsize=12)
    fig.text(0.4, 0.20, f"Difference Center - Around ={diff_value[0]:.3f}", ha="center", fontsize=12)
    diffs.append(diff_value)
    #outdir = os.path.join(os.path.expanduser('~'),"Output_expe","In_Vivo", "quick-spatial_map")
    #base = os.path.basename(filename)   # get file name only
    #name, _ = os.path.splitext(base)    # drop extension
    #save_path = os.path.join(outdir, f"{name}_qsm.png")
    
    #fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return 0

diffs = []
plot_qsm(5, coord_map, diffs)

# %%

summary = ep.compute_summary_data(stat_test_props={})
center_cond = (summary['x-center']==0) & (summary['y-center']==0)
center_value = summary['value'][center_cond]
around_value = summary['value'][~center_cond].mean()
diff_value = center_value-around_value

center_resp = summary['significant'][center_cond]

print("center value : ",center_value)
print("center significant : ", center_resp)
print("mean value around : ",around_value)
print("difference center - around : ",diff_value)
#%%
plot_qsm(index=0, coord_map=coord_map, diffs=diffs)
#%%
print(SESSIONS['files'])
print(SESSIONS['files'][0][-1])

#%% [markdown]
# ## Plot for all files
#%%
diffs = []
for index in range(len(SESSIONS['files'])):
    plot_qsm(index=index, coord_map=coord_map, diffs=diffs)


# %%
print(diffs)

plt.scatter(np.arange(len(diffs)), diffs)

# %% [markdown]
# ### plot per animal
#%%

# Example: diffs and matching animal IDs
animal_ids = ["1", "2", "3", "2", "3", 
              "4", "5", "6", "7", "8",
              "9", "10", "7", "2", "3",
              "4", "5", "9", "10", "10",
              "7", "5", "5", "2"]

# Find unique animals
unique_animals = np.unique(animal_ids)

means = []
for animal in unique_animals:
    animal_diffs = [d for d, a in zip(diffs, animal_ids) if a == animal]
    means.append(np.mean(animal_diffs))

plt.figure(figsize=(8,5))

# Bar plot of mean per animal
plt.bar(range(len(unique_animals)), means, alpha=0.6, color="skyblue", label="Mean per animal")

# Overlay scatter of individual points
for i, animal in enumerate(unique_animals):
    animal_diffs = [d for d, a in zip(diffs, animal_ids) if a == animal]
    # Add some jitter so dots don’t overlap
    x = np.random.normal(i, 0.05, size=len(animal_diffs))
    plt.scatter(x, animal_diffs, color="black", alpha=0.8)
    # Color points: red if < 0.05, black otherwise
    colors = ["red" if val > 0.03 else "black" for val in animal_diffs]
    plt.scatter(x, animal_diffs, color=colors, alpha=0.8)

# Add horizontal line at 0.05
plt.axhline(0.03, color="gray", linestyle="--", linewidth=1)

plt.xticks(range(len(unique_animals)), unique_animals)
plt.xlabel("Animal")
plt.ylabel("Difference")
plt.title("Differences per animal")


def responsive_center():
    files = []
    
    return files