import numpy as np
import os, sys
sys.path.append('../physion/src') # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')
from scipy import stats
import matplotlib.pyplot as plt
from scipy.optimize import minimize

    
# %%

# custom dFoF:

from physion.imaging.Calcium import compute_F0


dFoF_options = dict(\
    method_for_F0='sliding_minmax',
    #method_for_F0='sliding_percentile',
    #method_for_F0='percentile', #more strict
    sliding_window= 300,
    percentile=10,
    roi_to_neuropil_fluo_inclusion_factor=1.,
    neuropil_correction_factor=0.7, 
    with_computed_neuropil_fact=True)




# TO LOOP OVER NWB FILES WITH VISUAL STIMULUS --- DRIFITING GRATING ---  multisession

folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna',
                        'PN_cond-NDNF-CB1_WT-vs-KD', 'NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol='Natural-Images-4-repeats')



# STATISTICS PROPERTIES --- DRIFTING GRATINGS ---
stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='ttest',
                       sign='both')

pos_stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='ttest',
                       sign ='positive')
response_significance_threshold =0.05

neg_stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='ttest',
                       sign ='negative')

#response_significance_threshold = 0.05

# PLOT PROPERTIES --- DRIFTING GRATINGS ---

plot_props = dict(column_key='Image-ID',
                  with_annotation=True,
                  Ybar=0.5, Ybar_label="0.5$\Delta$F/F",
                  Xbar=0.5, Xbar_label="0.5s",
                  figsize=(13,1.8))


# RESPONSE ARGUMENTS --- DRIFTING GRATINGS ---



response_args = dict(quantity='dFoF')

summary_stats = []

RUNNING_SPEED_THRESHOLD = 0.5
NMIN_ROIS = 3
NMIN_EPISODES = 2

means = {} # 
for virus in ['sgRosa', 'sgCnr1']:
      for cond in ['all', 'aroused', 'still']:
              
                means['%s-%s' % (virus, cond)] = []


run_means = {} # 
for virus in ['sgRosa', 'sgCnr1']:
      for cond in ['all', 'aroused', 'still']:
              
               run_means['%s-%s' % (virus, cond)] = []

pupil_means = {} # 
for virus in ['sgRosa', 'sgCnr1']:
      for cond in ['all', 'aroused', 'still']:
              
               pupil_means['%s-%s' % (virus, cond)] = []



pos_means = {} # 
for virus in ['sgRosa', 'sgCnr1']:
      for cond in ['all', 'aroused', 'still']:
              pos_means['%s-%s' % (virus, cond)] = []

neg_means = {} # 
for virus in ['sgRosa', 'sgCnr1']:
      for cond in ['all', 'aroused', 'still']:
              
                neg_means['%s-%s' % (virus, cond)] = []



percentages = {}
for virus in ['sgRosa', 'sgCnr1']:
        percentages['%s' % virus] = []

pos_percentages = {}
for virus in ['sgRosa', 'sgCnr1']:
                pos_percentages['%s' % virus] = []


neg_percentages = {}
for virus in ['sgRosa', 'sgCnr1']:
      neg_percentages['%s' % virus] = []


for i,filename in enumerate(DATASET['files']):
    
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    
    print(i+1,'--', filename, '--', data.nROIs)
    # print(data.protocols)

    data.build_dFoF(**dFoF_options, verbose=True)
    data.build_pupil_diameter()
    data.build_facemotion()
    data.build_running_speed()
    
    if data.nROIs>0 and hasattr(data, 'pupil_diameter'):

        ep = physion.analysis.episodes.build.EpisodeData(data, 
                                                        quantities=['dFoF', 'running_speed','pupil_diameter'],
                                                        protocol_name='Natural-Images-4-repeats')
        
        if 'Image-ID' in ep.varied_parameters:
               
                # determine virus        
                if 'sgRosa' in data.nwbfile.virus:
                        virus = 'sgRosa'
                elif 'sgCnr1':
                        virus = 'sgCnr1'

                # 1) identify visually-responsive cells
                evokedStats = ep.pre_post_statistics_over_cells(\
                                                        stat_test_props,
                                                        response_args=\
                                                        dict(quantity='dFoF'),
                                                        response_significance_threshold=response_significance_threshold,
                                                        repetition_keys=['Image-ID','repeat']
                                                        )
                pos_evokedStats = ep.pre_post_statistics_over_cells(\
                                                                pos_stat_test_props,
                                                                response_args=\
                                                                dict(quantity='dFoF'),
                                                                response_significance_threshold=response_significance_threshold,
                                                                repetition_keys=['Image-ID','repeat']
                                                                )

                neg_evokedStats = ep.pre_post_statistics_over_cells(\
                                                                neg_stat_test_props,
                                                                response_args=\
                                                                dict(quantity='dFoF'),
                                                                response_significance_threshold=response_significance_threshold,
                                                                repetition_keys=['Image-ID','repeat']
                                                                )

                
                # 2) split rest / run
                withinEpisode = (ep.t>0) & (ep.t<ep.time_duration[0]) # 2 conditions: episode must have values over zero and ???

                preEpisode = 0
                Ep_run_speed = ep.running_speed[:,withinEpisode].mean(axis=1)
                run = np.mean(ep.running_speed[:,withinEpisode], axis=1) > RUNNING_SPEED_THRESHOLD
                Ep_pupil_size = ep.pupil_diameter[:,withinEpisode].mean(axis=1)
                pupil_bins = np.linspace(Ep_pupil_size.min(), Ep_pupil_size.max(), 15)
                # binning the data according to pupil level for analysis:
                

                bins = np.digitize(Ep_pupil_size, pupil_bins)
                speed_binned, sb_std = np.zeros(len(pupil_bins)), np.zeros(len(pupil_bins))

                for b in np.unique(bins):
                        speed_binned[b-1] = np.mean(Ep_run_speed[bins==b])
                        sb_std[b-1] = np.std(Ep_run_speed[bins==b])
                
                def func(t, X):
                        """ threshold-linear function """
                        return np.array([X[1]*(tt-X[0]) if tt>X[0] else 0 for tt in t])
    
                def to_minimize(X):
                        return np.sum((speed_binned-func(pupil_bins, X))**2)
                
                res = minimize(to_minimize,[pupil_bins.mean(), 1])
                
                pupil_threshold = res.x[0]

                pupil = np.mean(ep.pupil_diameter[:,withinEpisode], axis=1) > pupil_threshold

                aroused = pupil
                
                responsiveROIs = evokedStats['significant'].flatten()
                pos_responsiveROIs = pos_evokedStats['significant'][:,:].flatten()
                neg_responsiveROIs = neg_evokedStats['significant'][:,:].flatten()
                responsive = np.sum(responsiveROIs)/len(responsiveROIs)*100
                percentages['%s' % virus].append(responsive)
                pos_percentages['%s' % virus].append(np.sum(pos_responsiveROIs)/len(pos_responsiveROIs)*100)
                neg_percentages['%s' % virus].append(np.sum(neg_responsiveROIs)/len(neg_responsiveROIs)*100)
               
                                                   
                
                
                        
                #image_cond = (getattr(ep, 'Image-ID')==img_id)

        #
                #print("for session: %s" % filename)
                for cond, filter in zip(['all', 'aroused', 'still'],
                                        [aroused|~aroused, aroused, ~aroused]):
                        
                        if (np.sum(responsiveROIs)>=NMIN_ROIS) and \
                                (np.sum( filter)>= NMIN_EPISODES):
                                print("cond: %s-%s -> included %i ROIs and %i episodes" % (virus,cond, np.sum(responsiveROIs), np.sum( filter)))
                                print("cond: %s-%s -> %i ROIs out of %i ROIs are responsive" % (cond,virus, np.sum(responsiveROIs), len(responsiveROIs)))
                                
                                means['%s-%s' % (virus, cond)].append(
                                        ep.dFoF[filter, :, :][:, responsiveROIs, :])
                                
                                run_means['%s-%s' % (virus, cond)].append(
                                        ep.running_speed[filter, :])
                                pupil_means['%s-%s' % (virus, cond)].append(
                                        ep.pupil_diameter[filter, :])
                                
                                pos_means['%s-%s' % (virus, cond)].append(
                                        ep.dFoF[ filter, :, :][:, pos_responsiveROIs, :])
                                
                                neg_means['%s-%s' % (virus, cond)].append(
                                        ep.dFoF[filter, :, :][:, neg_responsiveROIs, :])
                                
                                
                                
                        
                                
                        else:
                                print("cond: %s -> [XX] response not included (%i ROIs, %i eps)" % (cond, np.sum(responsiveROIs), np.sum( filter)))
                

# now "means" is a list (over sessions) 
#    of responses of shape (episodes, responsiveROIs, time)

#%%
# 
if 'PN_cond-NDNF-CB1_WT-vs-KD' in folder:
       neuron= 'PN-cond-NDNF-CB1'
elif 'NDNF-cond-CB1_WT-vs-KD':
       neuron= 'NDNF-cond-CB1'

firgurename = 'dfof_beh_mod_all_natimg'+ neuron + '.svg'
figurepath = '/Users/macbookair/work/Figures/Natural_images/'

from scipy.stats import sem

fig, AX = pt.figure(axes=(3,1))

NMIN_SESSIONS = 0

for j, cond in enumerate(['all', 'aroused', 'still']):
    
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):


                session_responses = [np.mean(m,axis=(0,1))\
                        for m in means['%s-%s' % (virus, cond)]]

                if len(means['%s-%s' % (virus,cond)])>1:
                        #print(means['%s-%s' % (virus,cond)])
                        #print(virus,cond)
                        #print(session_responses)
                        pt.plot(ep.t, 
                                np.mean(session_responses,axis=0),
                                sy=sem(session_responses,axis=0),
                                color=color, ax=AX[j])
                        
                pt.annotate(AX[j],
                                'N=%i' % len(means['%s-%s'%(virus,cond)])+k*'\n',
                                (0.2,0.6), ha='right',
                                color=color, fontsize=6)
        
                pt.annotate(AX[j], cond, (0.5, 1))
     

        pt.set_plot(AX[j], 
                xlabel='time (s)',
                ylabel='$\\Delta$F/F' if j==0 else '')
pt.set_common_ylims(AX)   

#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')
#%%
baselineCond = (ep.t>-0.1) & (ep.t<0)

fig, AX = pt.figure(axes=(3,1))

NMIN_SESSIONS = 0

for j, cond in enumerate(['all', 'aroused', 'still']):
    
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):

                session_responses = []
                for m in means['%s-%s' % (virus, cond)]:
                        current_mean = np.mean(m,axis=(0,1))
                        current_mean_centered = current_mean - current_mean[baselineCond].mean()
                        session_responses.append(current_mean_centered)
                
                ##session_responses = [np.mean(m,axis=(0,1))\
                #        for m in means['%s-%s' % (virus, cond)]]
                
                

                if len(means['%s-%s' % (virus,cond)])>1:
                        #print(means['%s-%s' % (virus,cond)])
                        #print(virus,cond)
                        #print(session_responses)
                        pt.plot(ep.t, 
                                np.mean(session_responses,axis=0),
                                sy=sem(session_responses,axis=0),
                                color=color, ax=AX[j])
                        
                pt.annotate(AX[j],
                                'N=%i' % len(means['%s-%s'%(virus,cond)])+k*'\n',
                                (0.2,0.6), ha='right',
                                color=color, fontsize=6)
        
                pt.annotate(AX[j], cond, (0.5, 1))
     

        pt.set_plot(AX[j], 
                xlabel='time (s)',
                ylabel='$\\Delta$F/F' if j==0 else '')
pt.set_common_ylims(AX)   

#%%

baselineCond = (ep.t>-0.1) & (ep.t<0)

fig, AX = pt.figure(axes=(3,1))

NMIN_SESSIONS = 0

for j, cond in enumerate(['all', 'aroused', 'still']):
    
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):

                session_responses = [np.mean(m,axis=(0,1))\
                        for m in means['%s-%s' % (virus, cond)]]
                
                

                if len(means['%s-%s' % (virus,cond)])>1:
                        #print(means['%s-%s' % (virus,cond)])
                        #print(virus,cond)
                        #print(session_responses)
                        pt.plot(ep.t, 
                                np.mean(session_responses,axis=0) - np.mean(session_responses,axis=0)[baselineCond].mean(),
                                sy=sem(session_responses,axis=0),
                                color=color, ax=AX[j])
                        
                pt.annotate(AX[j],
                                'N=%i' % len(means['%s-%s'%(virus,cond)])+k*'\n',
                                (0.2,0.6), ha='right',
                                color=color, fontsize=6)
        
                pt.annotate(AX[j], cond, (0.5, 1))
     

        pt.set_plot(AX[j], 
                xlabel='time (s)',
                ylabel='$\\Delta$F/F' if j==0 else '')
pt.set_common_ylims(AX)   

#%%

fig, AX = pt.figure(axes=(3,1))

NMIN_SESSIONS = 0

for j, cond in enumerate(['all', 'aroused', 'still']):
    
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):


                session_responses = [np.mean(m,axis=(0))\
                        for m in run_means['%s-%s' % (virus, cond)]]

                if len(run_means['%s-%s' % (virus,cond)])>1:
                        #print(means['%s-%s' % (virus,cond)])
                        #print(virus,cond)
                        #print(session_responses)
                        pt.plot(ep.t, 
                                np.mean(session_responses,axis=0),
                                sy=sem(session_responses,axis=0),
                                color=color, ax=AX[j])
                        
                pt.annotate(AX[j],
                                'N=%i' % len(run_means['%s-%s'%(virus,cond)])+k*'\n',
                                (0.2,0.6), ha='right',
                                color=color, fontsize=6)
        
                pt.annotate(AX[j], cond, (0.5, 1))
     

        pt.set_plot(AX[j], 
                xlabel='time (s)',
                ylabel='speed(cm/s)' if j==0 else '')
pt.set_common_ylims(AX)   

#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')

#%%

fig, AX = pt.figure(axes=(3,1))

NMIN_SESSIONS = 0

for j, cond in enumerate(['all', 'aroused', 'still']):
    
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):


                session_responses = [np.mean(m,axis=(0))\
                        for m in pupil_means['%s-%s' % (virus, cond)]]

                if len(pupil_means['%s-%s' % (virus,cond)])>1:
                        #print(means['%s-%s' % (virus,cond)])
                        #print(virus,cond)
                        #print(session_responses)
                        pt.plot(ep.t, 
                                np.mean(session_responses,axis=0),
                                sy=sem(session_responses,axis=0),
                                color=color, ax=AX[j])
                        
                pt.annotate(AX[j],
                                'N=%i' % len(run_means['%s-%s'%(virus,cond)])+k*'\n',
                                (0.2,0.6), ha='right',
                                color=color, fontsize=6)
        
                pt.annotate(AX[j], cond, (0.5, 1))
     

        pt.set_plot(AX[j], 
                xlabel='time (s)',
                ylabel='pupil_diameter' if j==0 else '')
pt.set_common_ylims(AX)   

   
   

# %%
fig, AX = pt.figure(axes=(2,1))
firgurename = 'pie_all_natimg'+ neuron + '.svg'

NMIN_SESSIONS = 2

for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['blue','darkred']):
        
        
        perc_resp_ROI = np.mean(percentages['%s' % virus],axis=0)
        rest = 100-perc_resp_ROI
        print(perc_resp_ROI,'%s' % virus)
        pt.pie(data=[perc_resp_ROI,rest],
                COLORS=[color,'grey'],
                ext_labels = ['resp.','non-resp.'],
                
                pie_labels = ['%.1f'%perc_resp_ROI,'%.1f'%rest],
                title = '%s' % virus,
                ax=AX[k])
                #pt.annotate(AX[k][i],)
plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')
                                
# %%
# %% [markdown]
# ## Responsive neurons -  POSITIVE RESPONSES
#%%
firgurename = 'dfof_beh_mod_onlypos_natimg'+ neuron + '.svg'
figurepath = '/Users/macbookair/work/Figures/Natural_images/'


fig, AX = pt.figure(axes=(3,1))

NMIN_SESSIONS = 1
  
for j, cond in enumerate(['all', 'aroused', 'still']):
   
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):

                session_responses = [np.mean(m,axis=(0,1))\
                        for m in pos_means['%s-%s' % (virus, cond)]]
                
                if len(session_responses)>=NMIN_SESSIONS:
                        pt.plot(ep.t, 
                                np.mean(session_responses, axis=0),
                                #sy=sem(session_responses, axis=0),
                                color=color, ax=AX[j])
                        
                pt.annotate(AX[j],
                                'N=%i' % len(means['%s-%s'%(virus,cond)])+k*'\n',
                                (0.2,0.6), ha='right',
                                color=color, fontsize=6)
        
                pt.annotate(AX[j], cond, (0.5, 1))
     

        pt.set_plot(AX[j], 
                xlabel='time (s)',
                ylabel='$\\Delta$F/F' if j==0 else '',
                title= 'pos. resp.' if j==0 else '')
pt.set_common_ylims(AX)

plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')

# %% [markdown]
# ## Responsive neurons -  NEGATIVE RESPONSES
#%%

firgurename = 'dfof_beh_mod_onlyneg_natimg_'+ neuron + '.svg'
fig, AX = pt.figure(axes=(3,1))

NMIN_SESSIONS = 0
  
for j, cond in enumerate(['all', 'aroused', 'still']):
   
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):

                session_responses = [np.mean(m,axis=(0,1))\
                        for m in neg_means['%s-%s' % (virus, cond)]]
                
                if len(session_responses)>=NMIN_SESSIONS:
                        pt.plot(ep.t, 
                                np.mean(session_responses, axis=0),
                                sy=sem(session_responses, axis=0),
                                color=color, ax=AX[j])
                        
                pt.annotate(AX[j],
                                'N=%i' % len(means['%s-%s'%(virus,cond)])+k*'\n',
                                (0.2,0.6), ha='right',
                                color=color, fontsize=6)
        
                pt.annotate(AX[j], cond, (0.5, 1))
     

        pt.set_plot(AX[j], 
                xlabel='time (s)',
                ylabel='$\\Delta$F/F' if j==0 else '',
                title= 'neg. resp.' if j==0 else '')
pt.set_common_ylims(AX)

plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')

#%%
firgurename = 'pie_resp_neg_natimg_'+ neuron + '.svg'
fig, AX = pt.figure(axes=(2,1))

NMIN_SESSIONS = 1

for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['darkgrey','red']):
        
        
                perc_pos_resp_ROI = np.mean(pos_percentages['%s' % (virus)],axis=0)
                perc_neg_resp_ROI = np.mean(neg_percentages['%s' % (virus)],axis=0)
                perc_resp_ROI = perc_neg_resp_ROI + perc_pos_resp_ROI
                rest = 100-perc_resp_ROI
                print(perc_resp_ROI,'%s' % (virus))
                pt.pie(data=[perc_pos_resp_ROI,perc_neg_resp_ROI,rest],
                        COLORS=[color,'blue','lightgrey'],
                        #ext_labels = ['pos','neg','non-resp.'],
                        #ext_labels_distance=1.5,
                        pie_labels_distance=1,
                        
                        pie_labels = ['%.1f'%perc_pos_resp_ROI,'%.1f'%perc_neg_resp_ROI,'%.1f'%rest],
                        title = '%s' % (virus),
                        ax=AX[k])
                #pt.annotate(AX[k][i],)
plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')
# %%
