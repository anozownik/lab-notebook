import numpy as np
import os, sys
sys.path.append('../physion/src') # add src code directory for physion
import physion
import physion.utils.plot_tools as pt
pt.set_style('ticks')
from scipy import stats

    
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

response_significance_threshold = 0.05

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
      for cond in ['all', 'run', 'still']:
              for c, img_id in enumerate([1., 2., 3., 4., 5.]):
                means['%s-%s-%s' % (virus, cond, img_id)] = []
percentages = {}
for virus in ['sgRosa', 'sgCnr1']:
      
        for c, img_id in enumerate([1., 2., 3., 4., 5.]):
                percentages['%s-%s' % (virus, img_id)] = []

for i, filename in enumerate(DATASET['files']):
    
    data = physion.analysis.read_NWB.Data(filename,
                                    verbose=False)
    
    print(i+1, '--', filename, '--', data.nROIs)
    # print(data.protocols)

    data.build_dFoF(**dFoF_options, verbose=True)
    data.build_pupil_diameter()
    data.build_facemotion()
    data.build_running_speed()
    
    if data.nROIs>0:

        ep = physion.analysis.episodes.build.EpisodeData(data, 
                                                        quantities=['dFoF', 'running_speed'],
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
                                                        )
                
                # 2) split rest / run
                withinEpisode = (ep.t>0) & (ep.t<ep.time_duration[0])
                run = np.mean(ep.running_speed[:,withinEpisode], axis=1) > RUNNING_SPEED_THRESHOLD

                for img_id in ep.varied_parameters['Image-ID']:
                
                        # find responsive ROIs for this ImageID (from summary stats)
                        imageCond = evokedStats['Image-ID']==img_id
                        responsiveROIs = evokedStats['significant'][:,imageCond].flatten()
                        percentages['%s-%s' % (virus, img_id)].append(np.sum(responsiveROIs)/len(responsiveROIs)*100)
                        
                        
                               
                        image_cond = (getattr(ep, 'Image-ID')==img_id)

                #
                        #print("for session: %s" % filename)
                        for cond, filter in zip(['all', 'run', 'still'],
                                                [run|~run, run, ~run]):
                                
                                if (np.sum(responsiveROIs)>=NMIN_ROIS) and \
                                        (np.sum(image_cond & filter)>= NMIN_EPISODES):
                                        print("cond: %s-%s-%s -> included %i ROIs and %i episodes" % (virus,cond,img_id, np.sum(responsiveROIs), np.sum(image_cond & filter)))
                                        #print("cond: %s-%s-%s -> %i ROIs out of %i ROIs are responsive" % (cond,virus,img_id, np.sum(responsiveROIs), len(responsiveROIs)))
                                        
                                        means['%s-%s-%s' % (virus, cond, img_id)].append(
                                                ep.dFoF[image_cond & filter, :, :][:, responsiveROIs, :])
                                        
                                        
                                        
                                
                                        
                                else:
                                        print("cond: %s -> [XX] response not included (%i ROIs, %i eps)" % (cond, np.sum(responsiveROIs), np.sum(image_cond & filter)))
                

# now "means" is a list (over sessions) 
#    of responses of shape (episodes, responsiveROIs, time)

#%%
# 

from scipy.stats import sem

fig, AX = pt.figure(axes=(3,5))

NMIN_SESSIONS = 2

for j, cond in enumerate(['all', 'run', 'still']):
    for i, img_id in enumerate([1., 2., 3., 4., 5.]):
        for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['grey','darkred']):

                session_responses = [np.mean(m,axis=(0,1))\
                        for m in means['%s-%s-%s' % (virus,cond,img_id)]]
                
                if len(session_responses)>=NMIN_SESSIONS:
                        pt.plot(ep.t, 
                                np.mean(session_responses,axis=0),
                                sy=sem(session_responses, axis=0),
                                color=color, ax=AX[i][j])
                        
                pt.annotate(AX[i][j],
                            'N=%i' % len(session_responses)+k*'\n',
                            (0,0), ha='right',
                            color=color, fontsize=6)
        if i==0:
             pt.annotate(AX[i][j], cond, (0.5, 1))
        if j==0:
             pt.annotate(AX[i][j], 'Image-ID= %s ' % img_id,
                         (0,1), ha='right')

        pt.set_plot(AX[i][j], 
                xlabel='time (s)' if i==2 else '',
                ylabel='$\\Delta$F/F' if j==0 else '')
#pt.set_common_ylims(AX)        

#%%

fig, AX = pt.figure(axes=(5,2))


for k, virus, color in zip(range(2), ['sgRosa', 'sgCnr1'], ['blue','darkred']):
        for i, img_id in zip(range(5), ep.varied_parameters['Image-ID']):
        
                perc_resp_ROI = np.mean(percentages['%s-%s' % (virus,img_id)],axis=0)
                rest = 100-perc_resp_ROI
                print(perc_resp_ROI,'%s-%s' % (virus,img_id))
                pt.pie(data=[perc_resp_ROI,rest],
                        COLORS=[color,'grey'],
                        ext_labels = ['resp.','non-resp.'],
                        
                        pie_labels = ['%.1f'%perc_resp_ROI,'%.1f'%rest],
                        title = '%s-%s' % (virus,img_id),
                        ax=AX[k][i])
                #pt.annotate(AX[k][i],)
                                
  
   
   

# %%
