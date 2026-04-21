import numpy as np
from scipy.optimize import minimize

def remove_empty_sessions(d:dict):

    d_clean = {key : []  for key in d.keys()}
    
    for key in d.keys():
        for i in range(len(d[key])):
            if d[key][i].shape[1] :
                d_clean[key].append(d[key][i])
    
    return d_clean

def define_trials_arousal_state(episode, cond='speed', 
                                speed_threshold=0.5):
     
    if cond == 'speed':
        run, _ = get_trials_with_running(episode, speed_threshold)
        return ['all', 'run', 'still'], [run|~run, run, ~run]
    
    elif cond == 'pupil':
        dilated, _ = get_trials_with_dilated_pupil(episode, speed_threshold)
        return ['all', 'dilated', 'constricted'], [dilated|~dilated, dilated, ~dilated]
    
    elif cond == 'speed & pupil':
        dilated, _ = get_trials_with_dilated_pupil(episode, speed_threshold)
        return ['all', 'dilated & run', 'dilated & rest', 'constr. & rest', 'constr. & run'], \
            [run|~run, dilated & run, dilated & ~run, ~dilated & ~run, ~dilated & run]
    
    else:
        raise Exception("condition not defined. Choose between 'speed', 'pupil' or 'speed & pupil'")

def get_trials_with_running(episode, speed_threshold):

    if hasattr(episode, 'running_speed'):
         
        # Defining stimuli interval, excluding pre-stimulus and post-stimulus intervals
        withinEpisode = (episode.t>0) & (episode.t<episode.time_duration[0]) 

        Ep_run_speed = episode.running_speed[:,withinEpisode].mean(axis=1)

        # Trials with running
        run = Ep_run_speed > speed_threshold
    
    else :
        raise Exception("episode does not have running_speed attribute")

    return run, Ep_run_speed

def get_trials_with_dilated_pupil(episode, speed_threshold=0.5):

    _, Ep_run_speed = get_trials_with_running(episode, speed_threshold)

    if hasattr(episode, 'pupil_diameter'):
         
        # Defining stimuli interval, excluding pre-stimulus and post-stimulus intervals
        withinEpisode = (episode.t>0) & (episode.t<episode.time_duration[0]) 
        Ep_pupil_size = episode.pupil_diameter[:,withinEpisode].mean(axis=1)

        #------ Determine pupil threshold for arousal state classification ------

        # binning the data according to pupil level for analysis:
        pupil_bins = np.linspace(Ep_pupil_size.min(), Ep_pupil_size.max(), 15)

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

        # Trials with dilated pupil (aroused state)
        dilated = Ep_pupil_size > pupil_threshold
    
    else :
        raise Exception("episode does not have pupil_diameter attribute")

    return dilated, Ep_pupil_size