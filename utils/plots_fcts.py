import os, sys, tempfile
sys.path.append('../../physion/src') # add src code directory for physion

import physion.utils.plot_tools as pt
pt.set_style('ticks')

import numpy as np
from scipy.stats import sem
import matplotlib.pyplot as plt
from itertools import product
import matplotlib.colors as mcolors

# plot variables
color_virus = {'sgRosa' : 'grey', 
               'sgCnr1': "darkred"}

def plot_average_response(t, responses, 
                          viruses, states, varied_parameter=[],
                          vparam_name='',
                          included_mice=None,
                          nmin_sessions=1,
                          baselineSubtraction=False,
                          baselineCond=None,
                          annotation_props=dict(xy=(0.05,1), ha='left', fontsize=4),
                          savepath=None):
    
    if len(varied_parameter) == 0:
        fig, AX = plot_average_response_no_vparam(t, responses, 
                                                  viruses, states, 
                                                  included_mice, nmin_sessions, 
                                                  baselineSubtraction, baselineCond,
                                                  annotation_props)
    else :
        fig, AX = plot_average_response_vparam(t, responses, 
                                               viruses, states, varied_parameter, vparam_name,
                                               included_mice, nmin_sessions, 
                                               baselineSubtraction, baselineCond,
                                               annotation_props)

    if savepath is not None:
        plt.savefig(os.path.join(savepath), transparent=True, format='svg')

    return fig, AX

def plot_average_response_no_vparam(t, responses, 
                                    viruses, states, 
                                    included_mice=None,
                                    nmin_sessions=1,
                                    baselineSubtraction=False,
                                    baselineCond=None,
                                    annotation_props=dict(xy=(0.05,1), ha='left', fontsize=4)):
    
    fig, AX = pt.figure(axes=(len(states), 1))

    if len(states) == 1:
        AX = np.reshape(AX, (1))

    for j, state in enumerate(states):
        for k, virus in enumerate(viruses):

            key = f'{virus}-{state}'

            session_responses = np.array([np.mean(r, axis=(0,1)) for r in responses[key]])
            
            if baselineSubtraction:
                if baselineCond is None:
                    baselineCond = (t<0)
                if len(session_responses.shape) == 2:
                    session_responses = session_responses - session_responses[:, baselineCond].mean(axis=1, keepdims=True)

            if np.shape(session_responses)[0] >= nmin_sessions :

                if np.shape(session_responses)[0] == 1 :
                    pt.plot(t, session_responses[0],
                            color=color_virus[virus], ax=AX[j])
                else :
                    pt.plot(t, np.mean(session_responses,axis=0),
                            sy=sem(session_responses,axis=0),
                            color=color_virus[virus], ax=AX[j])
            
            if len(responses[key]) > 0:

                nb_eps, nb_rois = np.column_stack([(r.shape[0], r.shape[1]) for r in responses[key]])

                if included_mice is not None:
                    nb_mice = np.unique(included_mice[key]).shape[0]
                    pt.annotate(AX[j],
                                'N=%i (%i mice, %i rois, %i eps)' % (len(responses[key]), 
                                                                        nb_mice,
                                                                        nb_rois.sum(), 
                                                                        nb_eps.sum())
                                                                +k*'\n',
                                                                
                                color=color_virus[virus], **annotation_props)
                else :
                    pt.annotate(AX[j],
                                'N=%i (%i rois, %i eps)' % (len(responses[key]), 
                                                            nb_rois.sum(), 
                                                            nb_eps.sum())
                                                            +k*'\n',
                                                                
                                color=color_virus[virus], **annotation_props)          
            else :
                pt.annotate(AX[j], 'N=0' +k*'\n', color=color_virus[virus], **annotation_props)
    
            pt.annotate(AX[j], state, (0.5, 1.3), ha='center') 
        
        pt.set_plot(AX[j], xlabel='time (s)', ylabel='$\\Delta$F/F' if j==0 else '')
    
    #pt.set_common_ylims(AX)

    return fig, AX

def plot_average_response_vparam(t, responses, 
                                 viruses, states, varied_parameter=[],
                                 vparam_name='',
                                 included_mice=None,
                                 nmin_sessions=1,
                                 baselineSubtraction=False,
                                 baselineCond=None,
                                 annotation_props=dict(xy=(0.05,1), ha='left', fontsize=4)):
    

    fig, AX = pt.figure(axes=(len(states), len(varied_parameter)))

    if len(states) == 1:
        AX = np.reshape(AX, (-1, 1))

    for j, state in enumerate(states):
        for i, vparam in enumerate(varied_parameter):
            for k, virus in enumerate(viruses):

                key = f'{virus}-{state}-{vparam}'

                session_responses = np.array([np.mean(r, axis=(0,1)) for r in responses[key]])
                
                if baselineSubtraction and len(session_responses) > 0:
                    if baselineCond is None:
                        baselineCond = (t<0)
                    if len(session_responses.shape) == 2:
                        session_responses= session_responses - session_responses[:, baselineCond].mean(axis=1, keepdims=True)

                if np.shape(session_responses)[0] >= nmin_sessions :

                    if np.shape(session_responses)[0] == 1 :
                        pt.plot(t, session_responses[0],
                                color=color_virus[virus], ax=AX[i][j])
                    else :
                        pt.plot(t, np.mean(session_responses,axis=0),
                                sy=sem(session_responses,axis=0),
                                color=color_virus[virus], ax=AX[i][j])
                
                if len(responses[key]) > 0:

                    nb_eps, nb_rois = np.column_stack([(r.shape[0], r.shape[1]) for r in responses[key]])

                    if included_mice is not None:
                        nb_mice = np.unique(included_mice[key]).shape[0]
                        pt.annotate(AX[i][j],
                                    'N=%i (%i mice, %i rois, %i eps)' % (len(responses[key]), 
                                                                            nb_mice,
                                                                            nb_rois.sum(), 
                                                                            nb_eps.sum())
                                                                    +k*'\n',
                                                                    
                                    color=color_virus[virus], **annotation_props)
                    else :
                        pt.annotate(AX[i][j],
                                    'N=%i (%i rois, %i eps)' % (len(responses[key]), 
                                                                nb_rois.sum(), 
                                                                nb_eps.sum())
                                                                +k*'\n',
                                                                    
                                    color=color_virus[virus], **annotation_props)
                    
                else :
                    pt.annotate(AX[i][j], 'N=0' +k*'\n', color=color_virus[virus], **annotation_props)


                if i==0:
                    pt.annotate(AX[i][j], state, (0.5, 1.3), ha='center')
                if j==0:
                    pt.annotate(AX[i][j], '%s=%.1f ' % (vparam_name, vparam),
                                (0,1.1), ha='right')
                    
            pt.set_plot(AX[i][j], 
                        xlabel='time (s)' if i==len(varied_parameter)-1 else '', 
                        ylabel='$\\Delta$F/F' if j==0 else '')
    
    #pt.set_common_ylims(AX)

    return fig, AX

def pie_chart_responsive_neurons(percentages, viruses, varied_parameter=[], savepath=None):

    if len(varied_parameter) == 0:
        fig, AX = pie_chart_responsive_neurons_no_vparam(percentages, viruses)
    else :
        fig, AX = pie_chart_responsive_neurons_vparam(percentages, viruses, varied_parameter)

    if savepath is not None:
        plt.savefig(os.path.join(savepath), transparent=True, format='svg')

    return fig, AX

def pie_chart_responsive_neurons_no_vparam(percentages, viruses):

    fig, AX = pt.figure(axes=(1,len(viruses)))

    if len(viruses) == 1:
        AX = np.array([AX])

    for k, virus in enumerate(viruses):

        perc_resp_ROI = np.mean(percentages[virus], axis=0)
        rest = 100 - perc_resp_ROI
        #print(perc_resp_ROI,'%s' % virus)

        pt.pie(data=[perc_resp_ROI, rest],
               COLORS=[color_virus[virus], 'lightgrey'],
               ext_labels = ['resp.', 'non-resp.'],
               pie_labels = ['%.1f'%perc_resp_ROI,'%.1f'%rest],
               title = virus,
               ax=AX[k])
        
    return fig, AX

def pie_chart_responsive_neurons_vparam(percentages, viruses, varied_parameter):

    fig, AX = pt.figure(axes=(len(viruses), len(varied_parameter)))

    if len(viruses) == 1:
        AX = np.reshape(AX, (-1, 1))

    for k, virus in enumerate(viruses):

        perc_resp_ROI = np.mean(np.array(percentages[virus]), axis=0)

        for i, vparam in enumerate(varied_parameter):

            rest = 100 - perc_resp_ROI[i]
            #print(perc_resp_ROI[i],'%s-a=%.1f' % (virus, vparam))

            pt.pie(data=[perc_resp_ROI[i], rest],
                COLORS=[color_virus[virus], 'lightgrey'],
                ext_labels = ['resp.', 'non-resp.'],
                pie_labels = ['%.1f'%perc_resp_ROI[i], '%.1f'%rest],
                title = '%s-%.1f' % (virus, vparam),
                ax=AX[i][k])
        
    return fig, AX

def pie_chart_responsive_neurons_pos_neg(pos_percentages, neg_percentages, viruses, varied_parameter=[], savepath=None):

    if len(varied_parameter) == 0:
        fig, AX = pie_chart_responsive_neurons_pos_neg_no_vparam(pos_percentages, neg_percentages, viruses)
    else :
        fig, AX = pie_chart_responsive_neurons_pos_neg_vparam(pos_percentages, neg_percentages, viruses, varied_parameter)

    if savepath is not None:
        plt.savefig(os.path.join(savepath), transparent=True, format='svg')

    return fig, AX

def pie_chart_responsive_neurons_pos_neg_no_vparam(pos_percentages, neg_percentages, viruses):

    fig, AX = pt.figure(axes=(len(viruses), 1))

    if len(viruses) == 1:
        AX = np.array([AX])

    for k, virus in enumerate(viruses):

        perc_pos_resp_ROI = np.mean(pos_percentages[virus],axis=0)
        perc_neg_resp_ROI = np.mean(neg_percentages[virus],axis=0)
        perc_resp_ROI = perc_neg_resp_ROI + perc_pos_resp_ROI
        rest = 100 - perc_resp_ROI
        #print(perc_resp_ROI,'%s' % virus)

        pt.pie(data=[perc_pos_resp_ROI,perc_neg_resp_ROI,rest],
               COLORS=['tomato', 'royalblue','lightgrey'],
               #ext_labels = ['pos','neg','non-resp.'],
               # #ext_labels_distance=1.5,
               pie_labels_distance=1,
               pie_labels = ['%.1f'%perc_pos_resp_ROI, '%.1f'%perc_neg_resp_ROI, '%.1f'%rest],
               title = virus,
               ax=AX[k])
        
    return fig, AX

def pie_chart_responsive_neurons_pos_neg_vparam(pos_percentages, neg_percentages, viruses, varied_parameter):

    fig, AX = pt.figure(axes=(len(viruses), len(varied_parameter)))

    if len(viruses) == 1:
        AX = np.reshape(AX, (-1, 1))

    for k, virus in enumerate(viruses):

        perc_pos_resp_ROI = np.mean(np.array(pos_percentages[virus]), axis=0)
        perc_neg_resp_ROI = np.mean(np.array(neg_percentages[virus]), axis=0)

        for i, vparam in enumerate(varied_parameter):

            perc_resp_ROI = perc_neg_resp_ROI[i] + perc_pos_resp_ROI[i] 
            rest = 100 - perc_resp_ROI
            #print(perc_resp_ROI[i],'%s-a=%.1f' % (virus, vparam))

            pt.pie(data=[perc_pos_resp_ROI[i], perc_neg_resp_ROI[i], rest],
                    COLORS=['tomato', 'royalblue', 'lightgrey'],
                    #ext_labels = ['pos','neg','non-resp.'],
                    #ext_labels_distance=1.5,
                    pie_labels_distance=1.2,
                    
                    pie_labels = ['%.1f'%perc_pos_resp_ROI[i], '%.1f'%perc_neg_resp_ROI[i], '%.1f'%rest],
                    title = '%s-%.1f' % (virus, vparam),
                    ax=AX[i][k])
        
    return fig, AX

def plot_average_behavior(t, behaviors, viruses, states, 
                          nmin_sessions=1,
                          annotation_props=dict(xy=(0.05,1), ha='left', fontsize=5),
                          ylabel='',
                          savepath=None):

    fig, AX = pt.figure(axes=(len(states), 1))

    for j, state in enumerate(states):
        for k, virus in enumerate(viruses):

            key = f'{virus}-{state}'

            session_responses = np.array([np.mean(b, axis=(0)) for b in behaviors[key]])
            
            nb_eps = np.array([b.shape[0] for b in behaviors[key]])

            if np.shape(session_responses)[0] >= nmin_sessions :

                if np.shape(session_responses)[0] == 1 :
                    pt.plot(t, session_responses[0],
                            color=color_virus[virus], ax=AX[j])
                else :
                    pt.plot(t, np.mean(session_responses,axis=0),
                            sy=sem(session_responses,axis=0),
                            color=color_virus[virus], ax=AX[j])
            
            if len(behaviors[key]) > 0:
                nb_eps = np.array([b.shape[0] for b in behaviors[key]])

                pt.annotate(AX[j],
                            'N=%i ( %i eps)' % (len(behaviors[key]), nb_eps.sum())
                                                            +k*'\n',
                                                            
                            color=color_virus[virus], **annotation_props)          
            else :
                pt.annotate(AX[j], 'N=0' +k*'\n', color=color_virus[virus], **annotation_props)
    
            pt.annotate(AX[j], state, (0.5, 1.3), ha='center') 
        
        pt.set_plot(AX[j], xlabel='time (s)', ylabel=ylabel if j==0 else '')
    
    pt.set_common_ylims(AX)

    if savepath is not None:
        plt.savefig(os.path.join(savepath), transparent=True, format='svg')

    return fig, AX

def rastermap_session(session_id, ep, responses, viruses,
                      baselineSubtraction=False, baselineCond=None,
                      state_cond='all', 
                      savepath=None):
    
    responses_over_rois = []

    for i, virus in enumerate(viruses):

        if session_id[virus] >= len(responses[f'{virus}-{state_cond}']):
            raise Exception(f"Session {session_id[virus]} does not exist for {virus}-{state_cond}. Must be below {len(responses[f'{virus}-{state_cond}'])} ")

        if baselineSubtraction :
            responses_over_rois.append(np.subtract(np.mean(responses[f'{virus}-{state_cond}'][session_id[virus]], axis=1).T, 
                                                   np.mean(responses[f'{virus}-{state_cond}'][session_id[virus]][:, :, baselineCond], axis=(1,2))).T)
        else :
            responses_over_rois.append(np.mean(responses[f'{virus}-{state_cond}'][session_id[virus]], axis=1))

    responses_over_rois = np.array(responses_over_rois)

    vmin = np.min(responses_over_rois)
    vmax = np.max(responses_over_rois)

    fig, AX = plt.subplots(1, len(viruses), figsize=(8, 5))

    for i, virus in enumerate(viruses):

        im = AX[i].pcolormesh(ep.t, np.arange(responses[f'{virus}-{state_cond}'][session_id[virus]].shape[0]), 
                              responses_over_rois[i], vmin=vmin, vmax=vmax, cmap='magma')
        AX[i].set_title(f'{virus}-{state_cond}')

    fig.colorbar(im)

    if savepath is not None:
        plt.savefig(os.path.join(savepath), transparent=True, format='svg')

    return fig, AX

def sort_rois_by_mean_response(responses, response_window=None):
    if response_window is None:
        response_window = np.ones(responses.shape[1], dtype=bool)

    idx_sorted = np.argsort(np.mean(responses[:, response_window], axis=1))
    return idx_sorted

def plot_rastermap(responses, ep, viruses, state_cond='all', varied_parameters=[],
                   baselineSubtraction=False, baselineCond=None,
                   sort_fct=sort_rois_by_mean_response, sort_fcts_options=dict(response_window=None),
                   savepath=None):
    
    if len(varied_parameters) == 0:
        fig, AX = plot_rastermap_no_vparam(responses, ep, viruses, state_cond,
                                           baselineSubtraction, baselineCond,
                                           sort_fct, sort_fcts_options,
                                           savepath)
    else :
        fig, AX = plot_rastermap_vparam(responses, ep, viruses, state_cond, varied_parameters,
                                        baselineSubtraction, baselineCond,
                                        sort_fct, sort_fcts_options,
                                        savepath)

    if savepath is not None:
        plt.savefig(os.path.join(savepath), transparent=True, format='svg')

    return fig, AX

def plot_rastermap_no_vparam(responses, ep, viruses, state_cond='all',
                   baselineSubtraction=False, baselineCond=None,
                   sort_fct=sort_rois_by_mean_response, sort_fcts_options=dict(response_window=None),
                   savepath=None):
    
    responses_over_rois = {}
    vmin, vmax = np.inf, -np.inf

    for i, virus in enumerate(viruses):

        if baselineSubtraction :
            responses_over_rois[f'{virus}-{state_cond}'] = np.subtract(responses[f'{virus}-{state_cond}'].T, 
                                                                       np.mean(responses[f'{virus}-{state_cond}'][:, baselineCond], axis=1)).T
        else :
            responses_over_rois[f'{virus}-{state_cond}'] = responses[f'{virus}-{state_cond}']

        if np.min(responses_over_rois[f'{virus}-{state_cond}']) < vmin:
            vmin = np.min(responses_over_rois[f'{virus}-{state_cond}'])
        if np.max(responses_over_rois[f'{virus}-{state_cond}']) > vmax:
            vmax = np.max(responses_over_rois[f'{virus}-{state_cond}'])

        norm = mcolors.TwoSlopeNorm(
            vmin=vmin, 
            vcenter=0.0,
            vmax=vmax  
        )

    fig, AX = plt.subplots(1, len(viruses), figsize=(10, 5), sharex=True)

    for i, virus in enumerate(viruses):

        idx_sorted = sort_fct(responses_over_rois[f'{virus}-{state_cond}'], **sort_fcts_options)

        im = AX[i].pcolormesh(ep.t, np.arange(responses[f'{virus}-{state_cond}'].shape[0]), 
                              responses_over_rois[f'{virus}-{state_cond}'][idx_sorted, :], cmap='PiYG', 
                              norm=norm)
        AX[i].set_title(f'{virus}-{state_cond}')

        AX[i].vlines(0, -0.5, responses[f'{virus}-{state_cond}'].shape[0]-0.5, color='k', linewidth=1, linestyle='--')

        AX[i].vlines(ep.time_duration[0], -0.5, responses[f'{virus}-{state_cond}'].shape[0]-0.5, color='k', linewidth=1, linestyle='--')

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.38, 0.007, 0.5])
    fig.colorbar(im, cax=cbar_ax)

    if savepath is not None:
        plt.savefig(os.path.join(savepath), transparent=True, format='svg')

    return fig, AX

def plot_rastermap_vparam(responses, ep, viruses, state_cond='all', varied_parameters=[],
                   baselineSubtraction=False, baselineCond=None,
                   sort_fct=sort_rois_by_mean_response, sort_fcts_options=dict(response_window=None),
                   savepath=None):
    
    responses_over_rois = {}
    norms = []

    for i, vparam in enumerate(varied_parameters):

        vmin, vmax = np.inf, -np.inf

        for j, virus in enumerate(viruses): 

            key = f'{virus}-{state_cond}-{vparam}'

            if baselineSubtraction :
                responses_over_rois[key] = np.subtract(responses[key].T, 
                                                       np.mean(responses[key][:, baselineCond], axis=1)).T
            else :
                responses_over_rois[key] = responses[key]

            if responses_over_rois[key].shape[0] != 0:
                if np.min(responses_over_rois[key]) < vmin:
                    vmin = np.min(responses_over_rois[key])
                if np.max(responses_over_rois[key]) > vmax:
                    vmax = np.max(responses_over_rois[key])

        if responses_over_rois[key].shape[0] != 0:
            norm = mcolors.TwoSlopeNorm(
                vmin=vmin, 
                vcenter=0.0,
                vmax=vmax
            )
        else :
            norm = None

        norms.append(norm)

    fig, AX = plt.subplots(len(varied_parameters), len(viruses), figsize=(10, 5*len(varied_parameters)), sharex=True)

    for i, vparam in enumerate(varied_parameters):
        for j, virus in enumerate(viruses):

            key = f'{virus}-{state_cond}-{vparam}'

            if responses_over_rois[key].shape[0] != 0:
                idx_sorted = sort_fct(responses_over_rois[key], **sort_fcts_options)

                im = AX[i][j].pcolormesh(ep.t, np.arange(responses[key].shape[0]), 
                                    responses_over_rois[key][idx_sorted, :], cmap='PiYG', 
                                    norm=norms[i])
                

                AX[i][j].vlines(0, -0.5, responses[key].shape[0]-0.5, color='k', linewidth=1, linestyle='--')

                AX[i][j].vlines(ep.time_duration[0], -0.5, responses[key].shape[0]-0.5, color='k', linewidth=1, linestyle='--')

            AX[i][j].set_title(key)

        if responses_over_rois[key].shape[0] != 0: 
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.82, AX[i][0].get_position().y0, 0.007, AX[i][0].get_position().height])
            fig.colorbar(im, cax=cbar_ax)


    if savepath is not None:
        plt.savefig(os.path.join(savepath), transparent=True, format='svg')

    return fig, AX

def plot_dist_reliability(reliability, viruses, varied_parameter=[], vparam_name='', plot_type='violin', 
                          only_significant=True, significant_threshold=0.01, savepath=None):

    if len(varied_parameter) == 0:
        fig, AX = plot_dist_reliability_no_vparam(reliability, viruses, plot_type, only_significant, significant_threshold)
    else :
        fig, AX = plot_dist_reliability_vparam(reliability, viruses, varied_parameter, vparam_name, plot_type, only_significant, significant_threshold)

    if savepath is not None:
        plt.savefig(os.path.join(savepath), transparent=True, format='svg')

    return fig, AX

def plot_dist_reliability_no_vparam(reliability, viruses, plot_type='violin',
                          only_significant=True, significant_threshold=0.01):

    r_values = {virus: [] for virus in viruses}
    for v in viruses:

        for i in range(len(reliability[v])):

            if only_significant:
                significant = reliability[v][i]['pval'] <= significant_threshold
                r_values[v].append(reliability[v][i]['r'][significant])
            else :
                r_values[v].append(reliability[v][i]['r'])
        
        r_values[v] = np.concatenate(r_values[v])

    if plot_type == 'violin':

        fig, AX = plt.subplots(1, 1, figsize=(1.5, 1.5))

        for k, virus in enumerate(viruses):
            pt.violin(r_values[virus], x=k*1, color=color_virus[virus], ax=AX)

        pt.set_plot(AX, ['left'])

    elif plot_type == 'hist' or plot_type == 'histogram':

        fig, AX = pt.figure(axes=(len(viruses), 1))

        for k, virus in enumerate(viruses):
            ymax = np.histogram(r_values[virus], bins=20)[0].max()

            AX[k].hist(r_values[virus], bins=20, color=color_virus[virus])
            AX[k].vlines(np.mean(r_values[virus]), 0, ymax, color='black', linewidth=0.5, linestyle='dashed')
            AX[k].annotate('mean=%.2f' % np.mean(r_values[virus]), xy=(np.mean(r_values[virus])-0.01, ymax), 
                        ha='right', fontsize=4)
            AX[k].set_title(virus+ ' (n=%d)' % len(r_values[virus]), fontsize=6)
            AX[k].set_xlim(round(np.min(r_values[virus])*10)*0.1 - 0.1, 1)
    
    else :
        raise ValueError('plot_type value not recognized, should be "violin" or "hist"')

    return fig, AX

def plot_dist_reliability_vparam(reliability, viruses, varied_parameter, vparam_name, 
                                 plot_type='violin', only_significant=True, significant_threshold=0.01):

    r_values = {f"{virus}-{vparam}": [] for virus, vparam in product(viruses, varied_parameter)}

    for v in viruses:

        for k, vparam in enumerate(varied_parameter):

            for i in range(len(reliability[v])):
        
                if only_significant:
                    significant = reliability[v][i]['pval'][:, k] <= significant_threshold
                    r_values[f"{v}-{vparam}"].append(reliability[v][i]['r'][:, k][significant])
                else :
                    r_values[f"{v}-{vparam}"].append(reliability[v][i]['r'][:, k])
        
            r_values[f"{v}-{vparam}"] = np.concatenate(r_values[f"{v}-{vparam}"])

    if plot_type == 'violin':

        fig, AX = plt.subplots(len(varied_parameter), 1, figsize=(1.5, 1.5*len(varied_parameter)))

        for i, vparam in enumerate(varied_parameter):

            for k, virus in enumerate(viruses):
                pt.violin(r_values[f"{virus}-{vparam}"], x=k*1, color=color_virus[virus], ax=AX[i])

            pt.set_plot(AX[i], ['left'])
            pt.annotate(AX[i], '%s=%.0f ' % (vparam_name, vparam),
                                (0,1.05), ha='right', fontsize=7)

    elif plot_type == 'hist' or plot_type == 'histogram':

        fig, AX = pt.figure(axes=(len(viruses), len(varied_parameter)))

        for i, vparam in enumerate(varied_parameter):

            for k, virus in enumerate(viruses):
                ymax = np.histogram(r_values[f"{virus}-{vparam}"], bins=20)[0].max()

                AX[i][k].hist(r_values[f"{virus}-{vparam}"], bins=20, color=color_virus[virus])
                AX[i][k].vlines(np.mean(r_values[f"{virus}-{vparam}"]), 0, ymax, color='black', linewidth=0.5, linestyle='dashed')
                AX[i][k].annotate('mean=%.2f' % np.mean(r_values[f"{virus}-{vparam}"]), xy=(np.mean(r_values[f"{virus}-{vparam}"])-0.01, ymax), 
                            ha='right', fontsize=4)
                AX[i][k].set_title(virus+ ' (n=%d)' % len(r_values[f"{virus}-{vparam}"]), fontsize=6)
                AX[i][k].set_xlim(round(np.min(r_values[f"{virus}-{vparam}"])*10)*0.1 - 0.1, 1)

                if k == 0 :
                    pt.annotate(AX[i][k], '%s=%.0f ' % (vparam_name, vparam),
                                    (-0.25,1.), ha='right', fontsize=7)
    
    else :
        raise ValueError('plot_type value not recognized, should be "violin" or "hist"')

    return fig, AX
