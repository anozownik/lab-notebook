# %%
import numpy as np
import os, sys
sys.path.append('../physion/src') # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('manuscript')
from scipy import stats
import matplotlib.pyplot as plt


    
# %%

"""" VISION SURVEY SUBPROTOCOLS 
static-patch' 'looming-stim' 'Natural-Images-4-repeats'
'drifting-grating' 'drifting-surround' 'moving-dots' 'grey-10min'
'black-2min' 'quick-spatial-mapping']"""

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
    with_computed_neuropil_fact=False)




# TO LOOP OVER NWB FILES WITH VISUAL STIMULUS --- DRIFITING GRATING ---  multisession

folder = os.path.join(os.path.expanduser('~'), 'DATA', 'Adrianna',
                        'NDNF_cond-CB1_WT-vs-KD', 'NWBs')

DATASET = physion.analysis.read_NWB.scan_folder_for_NWBfiles(folder,
                                        for_protocol='moving-dots')



# STATISTICS PROPERTIES --- DRIFTING GRATINGS ---
stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='ttest',
                       sign ='both')
response_significance_threshold =0.05

neg_stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='ttest',
                       sign ='negative')


pos_stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='ttest',
                       sign ='positive')
# PLOT PROPERTIES --- DRIFTING GRATINGS ---

plot_props = dict(column_key='angle',
                  with_annotation=True,
                  Ybar=0.5, Ybar_label="0.5$\Delta$F/F",
                  Xbar=0.5, Xbar_label="0.5s",
                  figsize=(9,1.8))


# RESPONSE ARGUMENTS --- DRIFTING GRATINGS ---



response_args = dict(quantity='dFoF')

summary_stats = []

RUNNING_SPEED_THRESHOLD = 0.5
NMIN_EPISODES = 1
NMIN_ROIS = 1


# %%

pos_means = {} # 
for virus in ['sgRosa', 'sgCnr1']:
      for cond in ['all', 'run', 'still']:
              for c, speed in enumerate([60, 90]):
                pos_means['%s-%s-c=%.1f' % (virus, cond, speed)] = []

neg_means = {} # 
for virus in ['sgRosa', 'sgCnr1']:
      for cond in ['all', 'run', 'still']:
              for c, speed in enumerate([60, 90]):
                neg_means['%s-%s-c=%.1f' % (virus, cond, speed)] = []

means = {} # 
for virus in ['sgRosa', 'sgCnr1']:
      for cond in ['all', 'run', 'still']:
              for c, speed in enumerate([60, 90]):
                means['%s-%s-c=%.1f' % (virus, cond, speed)] = []




percentages = {}
for virus in ['sgRosa', 'sgCnr1']:
      
        for c, speed in enumerate([60, 90]):
                percentages['%s-c=%.1f' % (virus, speed)] = []



pos_percentages = {}
for virus in ['sgRosa', 'sgCnr1']:
      
        for c, speed in enumerate([60, 90]):
                pos_percentages['%s-c=%.1f' % (virus, speed)] = []


neg_percentages = {}
for virus in ['sgRosa', 'sgCnr1']:
      
        for c, speed in enumerate([60, 90]):
                neg_percentages['%s-c=%.1f' % (virus, speed)] = []


run_means = {} # 
for virus in ['sgRosa', 'sgCnr1']:
      for cond in ['all', 'run', 'still']:
        for c, speed in enumerate([60, 90]):

                run_means['%s-%s-c=%.1f' % (virus, cond, speed)] =[]


# if you want to trouble shoot stuff on only a couple of files, use the following. 
# acces in DATASET only first 2 files

"""for filename in DATASET['files'][]:
    
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    
    print( '--', filename, '--', data.nROIs)
    # print(data.protocols)"""


for i,filename in enumerate(DATASET['files']):
    
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    
    print( i+1,'--', filename, '--', data.nROIs)
    # print(data.protocols)

    data.build_dFoF(**dFoF_options, verbose=True)
    data.build_pupil_diameter()
    data.build_facemotion()
    data.build_running_speed()
    
    if data.nROIs>0:

        ep= physion.analysis.episodes.build.EpisodeData(data, 
                                                        quantities=['dFoF', 'running_speed'],
                                                        protocol_name='moving-dots')
        
        # determine virus        
        if 'sgRosa' in data.nwbfile.virus:
                virus = 'sgRosa'
        elif 'sgCnr1':
                virus = 'sgCnr1'

                
                
                
                
        if 'speed' in ep.varied_parameters:
                


                        # 1) identify visually-responsive cells
                        evokedStats = ep.pre_post_statistics_over_cells(\
                                                                stat_test_props,
                                                                response_args=\
                                                                dict(quantity='dFoF'),
                                                                response_significance_threshold=response_significance_threshold,
                                                                )
                        

                        pos_evokedStats = ep.pre_post_statistics_over_cells(\
                                                                pos_stat_test_props,
                                                                response_args=\
                                                                dict(quantity='dFoF'),
                                                                response_significance_threshold=response_significance_threshold,
                                                                )

                        neg_evokedStats = ep.pre_post_statistics_over_cells(\
                                                                neg_stat_test_props,
                                                                response_args=\
                                                                dict(quantity='dFoF'),
                                                                response_significance_threshold=response_significance_threshold,
                                                                )

                        
                        # 2) split rest / run
                        withinEpisode = (ep.t>0) & (ep.t<ep.time_duration[0])
                        run = np.mean(ep.running_speed[:,withinEpisode], axis=1) > RUNNING_SPEED_THRESHOLD

                        for speed in ep.varied_parameters['speed']:
                        
                                # find responsive ROIs for this contrast (from summary stats)
                                speedCond = evokedStats['speed']==speed
                                responsiveROIs = evokedStats['significant'][:,speedCond].flatten()
                                pos_responsiveROIs = pos_evokedStats['significant'][:,speedCond].flatten()
                                neg_responsiveROIs = neg_evokedStats['significant'][:, speedCond].flatten()
                                percentages['%s-c=%.1f' % (virus, speed)].append(np.sum(responsiveROIs)/len(responsiveROIs)*100)
                                pos_percentages['%s-c=%.1f' % (virus, speed)].append(np.sum(pos_responsiveROIs)/len(pos_responsiveROIs)*100)
                                neg_percentages['%s-c=%.1f' % (virus, speed)].append(np.sum(neg_responsiveROIs)/len(neg_responsiveROIs)*100)
                
                                # build contrast condition on the single trial episodes
                                speed_cond = (ep.speed==speed)

                                #
                                print("for session: %s" % filename)
                                for cond, filter in zip(['all', 'run', 'still'],
                                                        [run|~run, run, ~run]):
                                        
                                        if (np.sum(responsiveROIs)>=NMIN_ROIS) and \
                                                (np.sum(speed_cond & filter)>= NMIN_EPISODES):
                                                print("cond: %s-%s-c=%.1f -> included %i ROIs and %i episodes" % (cond,virus,speed, np.sum(responsiveROIs), np.sum(speed_cond & filter)))
                                                print("cond: %s-%s-c=%.1f -> %i ROIs out of %i ROIs are responsive" % (cond,virus,speed, np.sum(responsiveROIs), len(responsiveROIs)))
                                                
                                                means['%s-%s-c=%.1f' % (virus, cond, speed)].append(
                                                        ep.dFoF[speed_cond & filter, :, :][:, responsiveROIs, :])
                                                
                                                pos_means['%s-%s-c=%.1f' % (virus, cond, speed)].append(
                                                        ep.dFoF[speed_cond & filter, :, :][:, pos_responsiveROIs, :])
                                                
                                                neg_means['%s-%s-c=%.1f' % (virus, cond, speed)].append(
                                                        ep.dFoF[speed_cond & filter, :, :][:, neg_responsiveROIs, :])
                                                
                                                run_means['%s-%s-c=%.1f' % (virus, cond, speed)].append(
                                                        ep.running_speed[speed_cond & filter,:])
                                                
                                        
                                                
                                        else:
                                                print("cond: %s -> [XX] response not included (%i ROIs, %i eps)" % (cond, np.sum(responsiveROIs), np.sum(speed_cond & filter)))
                        print()    

# now "means" is a list (over sessions) 
#    of responses of shape (episodes, responsiveROIs, time)



# %% [markdown]
# ## Responsive neurons - all responses
#%%

if 'PN_cond-NDNF-CB1_WT-vs-KD' in folder:
       neuron= 'PN-cond-NDNF-CB1'
elif 'NDNF-cond-CB1_WT-vs-KD':
       neuron= 'NDNF-cond-CB1'


figurepath = '/Users/macbookair/work/Figures/'
firgurename = 'dfof_beh_mod_movdots_'+ neuron + '.svg'

from scipy.stats import sem

fig, AX = pt.figure(axes=(3,2))

NMIN_SESSIONS = 1
  
for j, cond in enumerate(['all', 'run', 'still']):
    for i, speed in zip(range(2), ep.varied_parameters['speed']):
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):

                session_responses = [np.mean(m,axis=(0,1))\
                        for m in means['%s-%s-c=%.1f' % (virus, cond, speed)]]
                
                if len(session_responses)>=NMIN_SESSIONS:
                        pt.plot(ep.t, 
                                np.mean(session_responses, axis=0),
                                sy=sem(session_responses, axis=0),
                                color=color, ax=AX[i][j])
                        
                pt.annotate(AX[i][j],
                            'N=%i' % len(session_responses)+k*'\n',
                            (0,0.6), #ha='right',
                            color=color, fontsize=6)
        if i==0:
             pt.annotate(AX[i][j], cond, (0.5, 1))
        if j==0:
             pt.annotate(AX[i][j], 'speed=%.1f ' % speed,
                         (0,1), ha='right')

        pt.set_plot(AX[i][j], 
                    xlabel='time (s)' if i==2 else '',
                    ylabel='$\\Delta$F/F' if j==0 else '')
#pt.set_common_ylims(AX)

#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')

#%%
firgurename = 'speed_beh_mod_movdots_'+ neuron + '.svg'

from scipy.stats import sem

fig, AX = pt.figure(axes=(3,2))

NMIN_SESSIONS = 1
  
for j, cond in enumerate(['all', 'run', 'still']):
    for i, speed in zip(range(2), ep.varied_parameters['speed']):
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):

                session_responses = [np.mean(m,axis=(0))\
                        for m in run_means['%s-%s-c=%.1f' % (virus, cond, speed)]]
                
                if len(session_responses)>=NMIN_SESSIONS:
                        pt.plot(ep.t, 
                                np.mean(session_responses, axis=0),
                                #sy=sem(session_responses, axis=0),
                                color=color, ax=AX[i][j])
                        
                pt.annotate(AX[i][j],
                            'N=%i' % len(session_responses)+k*'\n',
                            (0,0.6), #ha='right',
                            color=color, fontsize=6)
        if i==0:
             pt.annotate(AX[i][j], cond, (0.5, 1))
        if j==0:
             pt.annotate(AX[i][j], 'speed=%.1f ' % speed,
                         (0,1), ha='right')

        pt.set_plot(AX[i][j], 
                    xlabel='time (s)' if i==2 else '',
                    ylabel='$\\Delta$F/F' if j==0 else '')
pt.set_common_ylims(AX)


# %% [markdown]
# ## Responsive neurons -  POSITIVE RESPONSES
#%%
firgurename = 'dfof_beh_mod_onlypos_'+ neuron + '.svg'
fig, AX = pt.figure(axes=(3,2))

NMIN_SESSIONS = 1
  
for j, cond in enumerate(['all', 'run', 'still']):
    for i, speed in zip(range(2), ep.varied_parameters['speed']):
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):

                session_responses = [np.mean(m,axis=(0,1))\
                        for m in pos_means['%s-%s-c=%.1f' % (virus, cond, speed)]]
                
                if len(session_responses)>=NMIN_SESSIONS:
                        pt.plot(ep.t, 
                                np.mean(session_responses, axis=0),
                                #sy=sem(session_responses, axis=0),
                                color=color, ax=AX[i][j])
                        
                pt.annotate(AX[i][j],
                            'N=%i' % len(session_responses)+k*'\n',
                            (0,0.6), #ha='right',
                            color=color, fontsize=6)
        if i==0:
             pt.annotate(AX[i][j], cond, (0.5, 1))
        if j==0:
             pt.annotate(AX[i][j], 'angle=%.1f ' % speed,
                         (0,1), ha='right')

        pt.set_plot(AX[i][j], 
                    xlabel='time (s)' if i==2 else '',
                    ylabel='$\\Delta$F/F' if j==0 else '',
                    title = "pos. responses" if j==0 else'')
pt.set_common_ylims(AX)

#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')

# %% [markdown]
# ## Responsive neurons -  NEGATIVE RESPONSES
#%%

firgurename = 'dfof_beh_mod_onlyneg_'+ neuron + '.svg'
fig, AX = pt.figure(axes=(3,2))

NMIN_SESSIONS = 1
  
for j, cond in enumerate(['all', 'run', 'still']):
    for i, speed in zip(range(3), ep.varied_parameters['speed']):
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):

                session_responses = [np.mean(m,axis=(0,1))\
                        for m in neg_means['%s-%s-c=%.1f' % (virus, cond, speed)]]
                
                if len(session_responses)>=NMIN_SESSIONS:
                        pt.plot(ep.t, 
                                np.mean(session_responses, axis=0),
                                #sy=sem(session_responses, axis=0),
                                color=color, ax=AX[i][j])
                        
                pt.annotate(AX[i][j],
                            'N=%i' % len(session_responses)+k*'\n',
                            (0,0.6), #ha='right',
                            color=color, fontsize=6)
        if i==0:
             pt.annotate(AX[i][j], cond, (0.5, 1))
        if j==0:
             pt.annotate(AX[i][j], 'speed=%.1f ' % speed,
                         (0,1), ha='right')

        pt.set_plot(AX[i][j], 
                    xlabel='time (s)' if i==2 else '',
                    ylabel='$\\Delta$F/F' if j==0 else '',
                    title = "neg. responses" if j==0 else'')
pt.set_common_ylims(AX)
#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')


# %% [markdown]
# ## pie charts -  ALL RESPONSES
#%%
firgurename = 'pie_resp_neurons_'+ neuron + '.svg'
fig, AX = pt.figure(axes=(2,2))

NMIN_SESSIONS = 1

for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['blue','darkred']):
        for i, contrast in zip(range(2), ep.varied_parameters['speed']):
        
                perc_resp_ROI = np.mean(percentages['%s-c=%.1f' % (virus,speed)],axis=0)
                rest = 100-perc_resp_ROI
                print(perc_resp_ROI,'%s-c=%.1f' % (virus,speed))
                pt.pie(data=[perc_resp_ROI,rest],
                        COLORS=[color,'grey'],
                        ext_labels = ['resp.','non-resp.'],
                        
                        pie_labels = ['%.1f'%perc_resp_ROI,'%.1f'%rest],
                        title = '%s-a=%.1f' % (virus,speed),
                        ax=AX[k][i])
                #pt.annotate(AX[k][i],)

#plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')

# %% [markdown]
# ## pie charts - SPLIT BY POS AND NEG
#%% 
firgurename = 'pie_resp_posneg_'+ neuron + '.svg'
fig, AX = pt.figure(axes=(3,2))

NMIN_SESSIONS = 1

for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['darkgrey','red']):
        for i, contrast in zip(range(2), ep.varied_parameters['speed']):
        
                perc_pos_resp_ROI = np.mean(pos_percentages['%s-c=%.1f' % (virus,contrast)],axis=0)
                perc_neg_resp_ROI = np.mean(neg_percentages['%s-c=%.1f' % (virus,contrast)],axis=0)
                perc_resp_ROI = perc_neg_resp_ROI + perc_pos_resp_ROI
                rest = 100-perc_resp_ROI
                print(perc_resp_ROI,'%s-c=%.1f' % (virus,contrast))
                pt.pie(data=[perc_pos_resp_ROI,perc_neg_resp_ROI,rest],
                        COLORS=[color,'blue','lightgrey'],
                        #ext_labels = ['pos','neg','non-resp.'],
                        #ext_labels_distance=1.5,
                        pie_labels_distance=1.5,
                        
                        pie_labels = ['%.1f'%perc_pos_resp_ROI,'%.1f'%perc_neg_resp_ROI,'%.1f'%rest],
                        title = '%s-c=%.1f' % (virus,contrast),
                        ax=AX[k][i])
                #pt.annotate(AX[k][i],)
                               
                        
plt.savefig(os.path.join(figurepath+firgurename),transparent=True, format='svg')

#%%


fig, AX = pt.figure(axes=(2,1))

baselineCond = (epGrating.t>-0.1) & (epGrating.t<0)
for cond, color in zip(['run', 'rest'], ['tab:orange', 'tab:blue']):
        session_responses = [np.mean(m,axis=(0,1))\
                        for m in means['%s-%s-c=%.1f' % (virus, cond, contrast)]]
    
    
    
    
    
    
    pt.plot(ep, np.mean(Responses[cond], axis=0), 
            sy = stats.sem(Responses[cond], axis=0), color=color, no_set=True, ax=AX[0])

    Responses['baselineSubstr_%s' % cond] = [\
                r-r[baselineCond].min() for r in Responses[cond]]
    
    pt.plot(Responses['t'], np.mean(Responses['baselineSubstr_%s' % cond], axis=0), 
            sy = stats.sem(Responses['baselineSubstr_%s' % cond], axis=0), color=color, no_set=True, ax=AX[1])
    
pt.set_plot(AX[0], ylabel='$\delta$ $\Delta$F/F', 
            xlabel='time from stim. (s)' , xlim=[-1, Responses['t'][-1]],
            title='N=%i sessions' % len(Responses['run']))


