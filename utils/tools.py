def remove_empty_sessions(d:dict):

    d_clean = {key : []  for key in d.keys()}
    
    for key in d.keys():
        for i in range(len(d[key])):
            if d[key][i].shape[1] :
                d_clean[key].append(d[key][i])
    
    return d_clean