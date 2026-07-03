# custom dFoF:
dFoF_options = dict(\
    method_for_F0='sliding_minmax',
    roi_to_neuropil_fluo_inclusion_factor=1.6,
    roi_to_neuropil_fluo_inclusion_factor_metric='std',
    sliding_window= 300,
    neuropil_correction_factor=0.7, 
    with_computed_neuropil_fact=True)

# STATISTICS PROPERTIES --- NATURAL IMAGES ---
stat_test_props = dict(interval_pre=[-1.,0],                                   
                       interval_post=[1.,2.],                                   
                       test='wilcoxon',
                       sign='both')

response_significance_threshold = 0.05

# RESPONSE ARGUMENTS --- NATURAL IMAGES ---

response_args = dict(quantity='dFoF')
dt_sampling = 1/30*1e3
summary_stats = []

RUNNING_SPEED_THRESHOLD = 0.5
NMIN_ROIS = 3
NMIN_EPISODES = 2
NMIN_SESSIONS = 1