import numpy as np
import os, sys
sys.path.append('../../physion/src') # add src code directory for physion
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
                evokedStats = ep.pre_post_statistics(\
                                                        stat_test_props,
                                                        response_args=\
                                                        dict(quantity='dFoF'),
                                                        response_significance_threshold=response_significance_threshold,
                                                        loop_over_cells=True,
                                                        repetition_keys=['Image-ID','repeat']
                                                        )
                pos_evokedStats = ep.pre_post_statistics(\
                                                        pos_stat_test_props,
                                                        response_args=\
                                                        dict(quantity='dFoF'),
                                                        response_significance_threshold=response_significance_threshold,
                                                        loop_over_cells=True,
                                                        repetition_keys=['Image-ID','repeat']
                                                        )

                neg_evokedStats = ep.pre_post_statistics(\
                                                        neg_stat_test_props,
                                                        response_args=\
                                                        dict(quantity='dFoF'),
                                                        response_significance_threshold=response_significance_threshold,
                                                        loop_over_cells=True,
                                                        repetition_keys=['Image-ID','repeat']
                                                        )


                
                # 2) split rest / run
                withinEpisode = (ep.t<0) & (ep.t<ep.time_duration[0]) # 2 conditions: episode must have values over zero and ???

                
                Ep_run_speed = ep.running_speed[:,withinEpisode].mean(axis=1)
                run = np.mean(ep.running_speed[:,withinEpisode], axis=1) > RUNNING_SPEED_THRESHOLD
                
                

    
       
                
                responsiveROIs = evokedStats['significant'].flatten()
                pos_responsiveROIs = pos_evokedStats['significant'].flatten()
                neg_responsiveROIs = neg_evokedStats['significant'].flatten()
                responsive = np.sum(responsiveROIs)/len(responsiveROIs)*100
                percentages['%s' % virus].append(responsive)
                pos_percentages['%s' % virus].append(np.sum(pos_responsiveROIs)/len(pos_responsiveROIs)*100)
                neg_percentages['%s' % virus].append(np.sum(neg_responsiveROIs)/len(neg_responsiveROIs)*100)
               
                                                   
                
                
                        
                #image_cond = (getattr(ep, 'Image-ID')==img_id)

        #
                #print("for session: %s" % filename)
                for cond, filter in zip(['all', 'aroused', 'still'],
                                        [run|~run, run, ~run]):
                        
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

if 'PN_cond-NDNF-CB1_WT-vs-KD' in folder:
       neuron= 'PN'
       lieu= 'PN-cond-NDNF-CB1_WT-vs-KD'
elif 'NDNF-cond-CB1_WT-vs-KD':
       lieu= 'NDNF-cond-CB1_WT-vs-KD'
       neuron= 'NDNF'


filename = 'Natural-Images-4-repeats_' + neuron
filepath = os.path.join('/Users/macbookair/DATA/Adrianna/'+lieu + filename)

file =np.save(os.path.join(filepath+filename), means, allow_pickle=True )




# %%


npyfile = np.load('/Users/macbookair/DATA/Adrianna/PN-cond-NDNF-CB1_WT-vs-KDNatural-Images-4-repeats_PNNatural-Images-4-repeats_PN.npz',allow_pickle=True)
npyfile.files
# %%
