#!/usr/bin/env python3
# identify select matches to deal with multimatching svs

import pandas as pd

# *** selectMatches
# ** deal with multi-matching **
# create a dictioary of matches for a given id
# explore matching graph in a greedy manner focusing on events that are most similar and by caller_order and expanding fully connected match graphs
# outliers that are not fully connected will be considered non-matching with a given clique

def select(matches1,out_dir,sample,caller_order,verbose=True):

    outFile1 = out_dir + sample + '-matches2.txt'    

    if verbose:
        print('selecting best sv matches...', flush=True)            

    # ** identify select matches that represent selection for multi-matching svs **
    if verbose:
        print('building sv match graph...', end='',flush=True)            
    matches1.loc[:,'score'] = matches1[['diff1','diff2']].apply(lambda x: max(x),axis=1)  # used to weight matchIDs when multi-matching
    matches1 = matches1.sort_values(by=['score','matchID'])
    # ** greedy matching based on perfect matches and then by matchID **

    # ** build a graph of matches **
    graph = dict()
    for idx,row in matches1.iterrows():
        rowDict = row.to_dict()
        id1,id2,score = rowDict['id'],rowDict['id_b'],rowDict['score']

        if id1 not in graph:
            graph[id1] = {id2:'N'}
        else:
            graph[id1][id2] = 'N'

        if id2 not in graph:
            graph[id2] = {id1:'N'}
        else:
            graph[id2][id1] = 'N'
    # **
    
    # ** search for fully connected cliques **
    if verbose:
        print('identifying cliques...', end='',flush=True)
    matches1['isSelect'] = 'N'
    matches1['cliqueID'] = 'NA'
    matches1['other_variant_ids'] = 'NA'    
    for idx,row in matches1.iterrows():
        rowDict = row.to_dict()
        matchID,id1,id2,score = rowDict['matchID'],rowDict['id'],rowDict['id_b'],rowDict['score']
        partners1 = sorted(sorted(graph[id1].keys()) + [id1])
        partners2 = sorted(sorted(graph[id2].keys()) + [id2])
        common = sorted([id for id in partners1 if id in partners2])
        union = sorted(list(set(partners1 + partners2)))

        # identify representatively cliqueID by caller order
        cliqueID = None
        for caller in caller_order:
            cliqueIDs = sorted([id for id in common if caller in id])
            if len(cliqueIDs) >= 1:
                cliqueID = cliqueIDs[0]
                break
        assert cliqueID != None
        other_variant_ids = ','.join([id for id in common if id != cliqueID])

        # look to see if other connected matches are already part of a clique
        commonSelect = [id for id in common if (id in graph[id1] and graph[id1][id] == 'Y') or (id in graph[id2] and graph[id2][id] == 'Y')]
        unionSelect = [id for id in union if (id in graph[id1] and graph[id1][id] == 'Y') or (id in graph[id2] and graph[id2][id] == 'Y')]

        # scenario 1: isolated clique with all the same partners 
        # ** if clique then all partners will be the same **
        if partners1 == partners2:
            matches1.at[idx,'isSelect'] = 'Y'
            matches1.at[idx,'cliqueID'] = cliqueID
            matches1.at[idx,'other_variant_ids'] = other_variant_ids
            graph[id1][id2] = 'Y'
            graph[id2][id1] = 'Y'
            # print('1')

        # else partners are not the same
        else:
            # scenario 2: triplicate where not fully connected
            # if common partners are only id1 and id2 and neither is part of a clique then by order we start a clique
            if len(common) == 2 and len(unionSelect) == 0:
                matches1.at[idx,'isSelect'] = 'Y'
                matches1.at[idx,'cliqueID'] = cliqueID
                matches1.at[idx,'other_variant_ids'] = other_variant_ids                
                graph[id1][id2] = 'Y'
                graph[id2][id1] = 'Y'
                # print('2')
            # scenario 3: new clique and one partner is connected to an extension
            # common partners are id1,id2, and others then 
            elif len(common) > 2 and len(unionSelect) == 0:
                matches1.at[idx,'isSelect'] = 'Y'
                matches1.at[idx,'cliqueID'] = cliqueID
                matches1.at[idx,'other_variant_ids'] = other_variant_ids                                
                graph[id1][id2] = 'Y'
                graph[id2][id1] = 'Y'
                # print('3')
            # scenario 4: previously identified clique and one partner is connected to an extension
            # start a clique
            elif len(common) > 2 and len(commonSelect) > 0:
                matches1.at[idx,'isSelect'] = 'Y'
                matches1.at[idx,'cliqueID'] = cliqueID
                matches1.at[idx,'other_variant_ids'] = other_variant_ids
                graph[id1][id2] = 'Y'
                graph[id2][id1] = 'Y'
                # print('4')
            # scenario 5: extension of a clique and this edge is not directly in the clique (since common select is not)
            # prevent false extension
            elif len(common) >= 2 and len(unionSelect) > 0:
                continue
            else:
                print('i no understand')
                raise
            
    matches1.to_csv(outFile1,sep="\t",index = False)
    
    if verbose:
        print('done')
    return(matches1)
