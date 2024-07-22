import json
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import random
import sys

def connected_components(adj_mat):
    nV = len(adj_mat)
    visited = [False for _ in range(nV)]
    cc = []
    for start_vertex in range(nV):
        if not visited[start_vertex]:
            stack = [start_vertex]
            comp = []
            while stack:
                v = stack.pop()
                if not visited[v]:
                    visited[v] = True
                    comp.append(v)
                    for i in range(nV): 
                        if adj_mat[v, i] == 1 and not visited[i]:
                            stack.append(i)
            cc.append(comp)
    return cc

def network_ids(network):
    m_size = len(network['Node_list'])
    key_list = np.array(list(network['Node_list'].keys()))
    key_map = {key_list[i]:i for i in range(len(key_list))}
    map_key = {j:i for i,j in key_map.items()}
    return m_size,key_list,key_map,map_key

def alter_patterns(patterns, map_key,key_map,network,printing_flag = False,fname=None):
    removed_ids = []
    for i in range(len(patterns)):
        if len(patterns[i]) > 1:
            edges = []
            for v in [map_key[x] for x in patterns[i]]:
                for key in network['Vessel_list'].keys():
                    if v in network['Vessel_list'][key]:
                        edges.append(key)
            edges = list(np.unique(edges))
            edge_data = pd.DataFrame({'source':[network['Vessel_list'][x][0] for x in edges],'target':[network['Vessel_list'][x][1] for x in edges],'edge_key':edges})
            G = nx.MultiGraph()
            for _, row in edge_data.iterrows():
                G.add_edge(row['source'], row['target'], key=row['edge_key'])
            if printing_flag:
                fig,axs = plt.subplots(1,2)
                pos = nx.spring_layout(G)
                nx.draw_networkx_nodes(G, pos, node_color = ['red' if key_map[x] in patterns[i] else 'blue' for x in G.nodes()], ax=axs[0])
                nx.draw_networkx_labels(G, pos, font_size=8, ax=axs[0])
                nx.draw_networkx_edges(G, pos, ax=axs[0])
                labels = {tuple(sorted(network['Vessel_list'][i])):f"{i}={len(network['centerline'][i])}" for i in edges}
                nx.draw_networkx_edge_labels(G,pos,labels, ax=axs[0])
            while(len(list(nx.bridges(G))) < len(G.edges())):
                removable_edges = [x for x in G.edges(keys=True) if (x[0],x[1]) not in nx.bridges(G)]
                nodes_by_degree = sorted(nx.degree(G,np.unique([x for y in removable_edges for x in (y[0],y[1])])), key = lambda x: x[1], reverse = True)
                biggest_degree = nodes_by_degree[0][1]
                removable_nodes = [x[0] for x in nodes_by_degree if x[1] == biggest_degree]
                removable_edges = [x for x in removable_edges if x[0] in removable_nodes or x[1] in removable_nodes]
                removed_edge = random.choice(removable_edges)
                removed_ids.append(removed_edge[2])
                edges.remove(removed_edge[2])
                G.remove_edge(*removed_edge)
#                print(f"removed edge {removed_edge[2]}")
            if printing_flag:
                pos = nx.spring_layout(G)
                nx.draw_networkx_nodes(G, pos, node_color = ['red' if key_map[x] in patterns[i] else 'blue' for x in G.nodes()], ax=axs[1])
                nx.draw_networkx_labels(G, pos, font_size=8, ax=axs[1])
                nx.draw_networkx_edges(G, pos, ax=axs[1])
                labels = {tuple(sorted(network['Vessel_list'][i])):f"{i}={len(network['centerline'][i])}" for i in edges}
                nx.draw_networkx_edge_labels(G,pos,labels, ax=axs[1])
                if printing_flag.lower() == "save":
                    plt.savefig(f'{fname}/{fname}_pattern_{i}.png')
                elif printing_flag.lower() == "plot":
                    plt.show()
            plt.close()
    return removed_ids

def make_adjacency(network, mode = "Multi"):
    m_size = len(network['Node_list'])
    key_list = np.array(list(network['Node_list'].keys()))
    key_map = {key_list[i]:i for i in range(len(key_list))}
    adj_mat = np.zeros((m_size,m_size))
    for key in network['Vessel_list'].keys():
        vi,vf = network['Vessel_list'][key]
        adj_mat[key_map[vf],key_map[vi]] += 1
        adj_mat[key_map[vi],key_map[vf]] += 1
    if mode != "Multi":
        adj_mat[adj_mat != 0] = 1
    return adj_mat

def degree_counts(adj_mat):
    degrees = np.sum(adj_mat,axis=0)
    count_degree = {i:j for i,j in zip(*np.unique(degrees,return_counts=True))}
    print(count_degree)
    return degrees

def merge_2_degree(network,degrees,map_key):
    for remnode in [map_key[i] for i in range(len(degrees)) if degrees[i] == 2]:
        E = []
        for vid in network['Vessel_list'].keys():
            if remnode in network['Vessel_list'][vid]:
                E.append(vid)
        if len(E) == 2:
            inid = E[0]
            outid = E[1]
            if network['Vessel_list'][inid][0] == remnode:
                if network['Vessel_list'][outid][0] == remnode:
                    network['centerline'][outid] = network['centerline'][outid][::-1]
                    network['Vessel_list'][outid][0],network['Vessel_list'][outid][1] = network['Vessel_list'][outid][1],network['Vessel_list'][outid][0]
                inid,outid=outid,inid
            elif network['Vessel_list'][outid][1] == remnode:
                network['centerline'][outid] = network['centerline'][outid][::-1]
                network['Vessel_list'][outid][0],network['Vessel_list'][outid][1] = network['Vessel_list'][outid][1],network['Vessel_list'][outid][0]
#            print(f"For Node {remnode} - edges {inid},{outid}\nBefore:\n{inid}:{network['Vessel_list'][inid]}\n{outid}:{network['Vessel_list'][outid]}")
            network['centerline'][inid].append(network['Node_list'][remnode])
            network['centerline'][inid].extend(network['centerline'][outid])
            network['Vessel_list'][inid][1] = network['Vessel_list'][outid][1]
            network['Vessel_list'][outid] = network['Vessel_list'][inid]
            network['centerline'][outid] = network['centerline'][inid]
#            print(f"After:{outid}:{network['Vessel_list'][outid]}\n")
            del network['Vessel_list'][inid]
            del network['centerline'][inid]
            del network['Node_list'][remnode]

def remove_small_vessels(network, map_key, degrees, vessel_cutoff = 0):
    removed_ids = {}
    for node in [map_key[index] for index in range(len(degrees)) if degrees[index] == 1]:
        for vid in network['Vessel_list'].keys():
            if node in network['Vessel_list'][vid]:
                if len(network['centerline'][vid]) <= vessel_cutoff:
                    removed_ids[node] = vid
    for node,vessel in removed_ids.items():
        del network['Vessel_list'][vessel]
        del network['centerline'][vessel]
        del network['Node_list'][node]

def offending_nodes(network,degrees,map_key,cutoff = 4):
    high_nodes = [map_key[x] for x in range(len(degrees)) if degrees[x] == cutoff]
    edge_list = {}
    for i in high_nodes:
        edge_list[i] = []
        for g in network['Vessel_list'].keys():
                if i in network['Vessel_list'][g]:
                        edge_list[i].append(g)
                        print(g)
                        print(network['Vessel_list'][g])
                        print(len(network['centerline'][g]))
    return high_nodes,edge_list

def make_clusters(network,high_nodes,edge_list):
    T1 = [[False,False,False],
      [True,False,False],
      [False,True,False]]
    T2 = [[False,False,False],
        [False,False,False],
        [True,True,False]]
    T3 = [[False,False,False],
        [True,False,False],
        [True,False,False]]
    Q1 = [[False,False,False,False],
        [True,False,False,False],
        [False,True,False,False],
        [True,False,True,False]]
    Q2 = [[False,False,False,False],
        [False,False,False,False],
        [True,True,False,False],
        [True,True,False,False]]
    Q3 = [[False,False,False,False],
        [True,False,False,False],
        [True,False,False,False],
        [False,True,True,False]]
    current_node = int(sorted(list(network['Node_list'].keys()),key=lambda x: int(x.lstrip('v')))[-1].lstrip('v'))
    current_edge = int(sorted(list(network['Vessel_list'].keys()),key=lambda x: int(x.lstrip('e')))[-1].lstrip('e'))
    for problem_node in high_nodes:
        new_nodes = []
        for edge in edge_list[problem_node]:
            if len(network['centerline'][edge]):
                vi,vf = network['Vessel_list'][edge]
                index = 0 if vi == problem_node else -1
                vn_coord = network['centerline'][edge][index]
                # Add new nodes
                current_node+=1
                vn = f'v{current_node}'
                new_nodes.append(vn)
                network['Node_list'][vn] = vn_coord
                # Alter current edge
                network['Vessel_list'][edge] = [vi,vn]
                network['centerline'][edge] = network['centerline'][edge][1+index:len(network['centerline'][edge])+index]
                # Add new edge
                current_edge += 1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [vn,vf]
                network['centerline'][en] = []
        numnodes = len(new_nodes)
        dist_m = np.zeros((numnodes,numnodes))
        for i in range(numnodes):
            for j in range(i+1):
                dist_m[i,j] = np.sqrt(np.sum([(x-y)**2 for x,y in zip(network['Node_list'][new_nodes[i]],network['Node_list'][new_nodes[j]])]))
        if numnodes == 4:
            connectivity = [np.sum(dist_m[x]) for x in [Q1,Q2,Q3]]
            case = connectivity.index(min(connectivity))
            if case == 0:
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[0],new_nodes[1]]
                network['centerline'][en] = []
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[1],new_nodes[2]]
                network['centerline'][en] = []
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[2],new_nodes[3]]
                network['centerline'][en] = []
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[3],new_nodes[0]]
                network['centerline'][en] = []
            elif case == 1:
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[0],new_nodes[2]]
                network['centerline'][en] = []
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[2],new_nodes[1]]
                network['centerline'][en] = []
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[1],new_nodes[3]]
                network['centerline'][en] = []
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[3],new_nodes[0]]
                network['centerline'][en] = []
            else:
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[0],new_nodes[2]]
                network['centerline'][en] = []
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[2],new_nodes[3]]
                network['centerline'][en] = []
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[3],new_nodes[1]]
                network['centerline'][en] = []
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[1],new_nodes[0]]
                network['centerline'][en] = []
        elif numnodes == 3:
            connectivity = [np.sum(dist_m[x]) for x in [T1,T2,T3]]
            case = connectivity.index(min(connectivity))
            if case == 0:
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[0],new_nodes[1]]
                network['centerline'][en] = []
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[1],new_nodes[2]]
                network['centerline'][en] = []
            elif case == 1:
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[0],new_nodes[2]]
                network['centerline'][en] = []
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[2],new_nodes[1]]
                network['centerline'][en] = []
            else:
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[0],new_nodes[2]]
                network['centerline'][en] = []
                current_edge+=1
                en = f'e{current_edge}'
                network['Vessel_list'][en] = [new_nodes[0],new_nodes[1]]
                network['centerline'][en] = []
        elif numnodes == 2:
            current_edge+=1
            en = f'e{current_edge}'
            network['Vessel_list'][en] = [new_nodes[0],new_nodes[1]]
            network['centerline'][en] = []

def network_status(net,patterns,fname = 'print'):
    m_size = len(net['Node_list'])
    adj_mat = np.zeros((m_size,m_size))
    key_list = np.array(list(net['Node_list'].keys()))
    key_map = {key_list[i]:i for i in range(len(key_list))}
    for key in net['Vessel_list'].keys():
        vi,vf = net['Vessel_list'][key]
        adj_mat[key_map[vf],key_map[vi]] += 1
        adj_mat[key_map[vi],key_map[vf]] += 1
    degs = np.sum(adj_mat,axis=0)
    count_deg = {i:0 for i in np.unique(degs)}
    for i in np.sum(adj_mat,axis=0):
        count_deg[i]+=1
    ##### Report #####
    out_buffer = ""
    out_buffer+=f'Number of nodes = {len(net["Node_list"])}\n'
    out_buffer+=f'Number of vessels = {len(net["Vessel_list"])}\n\n'
    out_buffer+=f'Vessel Lengths\n'
    sizes = np.unique([len(x) for x in net['centerline'].values()],return_counts = True)
    out_buffer+="\t".join(sizes[0].astype('str'))+'\n'
    out_buffer+="\t".join(sizes[1].astype('str'))+'\n'
    out_buffer+=f'\nNetwork node degree count:\n'
    out_buffer+='\n'.join([f'{x}:{y}' for x,y in count_deg.items()])
    out_buffer+=f'\nNumber of conflicting patterns: {len(patterns)}'
    if fname == 'print':
        print(out_buffer)
    else:
        with open(fname,'w') as f:
            f.write(out_buffer)

def main():
#    fname = "P8_Segmented"
    fname = sys.argv[1]
    with open(f'{fname}_primary_network.json','r') as f:
        network = json.load(f)

    m_size,key_list,key_map,map_key = network_ids(network)

    perc_tot = m_size**2/2
    perc_count = 0
    dist_m = np.zeros((m_size,m_size))
    for i in range(m_size):
        for j in range(i,m_size):
            dist_m[i,j] = np.sqrt(np.sum([(x-y)**2 for x,y in zip(network['Node_list'][key_list[i]],network['Node_list'][key_list[j]])]))
            dist_m[j,i] = dist_m[i,j]
            perc_count+=1
            print(f'{100*perc_count/perc_tot:.2f}%',end="\r")

    print("Distances Calculated")

    adj_mat = dist_m.copy()
    adj_mat[adj_mat<2] = 1
    for i in range(m_size):
        adj_mat[i,i] = 0

    adj_mat[adj_mat>=2] = 0

    patterns = connected_components(adj_mat)
    removed_ids = alter_patterns(patterns,map_key,key_map,network)#,printing_flag="save",fname=fname)
    removed_ids = np.unique(removed_ids)

    print("Patterns Saved\nBefore Removal:")

    adj_mat = make_adjacency(network)
    degrees = degree_counts(adj_mat)

    for i in removed_ids:
        del network['Vessel_list'][i]
        del network['centerline'][i]

    print("After altering patterns (remove edges)")
    adj_mat = make_adjacency(network)
    degrees = degree_counts(adj_mat)
    m_size,key_list,key_map,map_key = network_ids(network)

    merge_2_degree(network,degrees,map_key)

    print("merge 2-degree nodes")
    adj_mat = make_adjacency(network)
    degrees = degree_counts(adj_mat)
    m_size,key_list,key_map,map_key = network_ids(network)

    remove_small_vessels(network,map_key,degrees,vessel_cutoff=3)

    print("After removal of nodes with degree 1 and centerline <= 3")
    adj_mat = make_adjacency(network)
    degrees = degree_counts(adj_mat)
    m_size,key_list,key_map,map_key = network_ids(network)

    merge_2_degree(network,degrees,map_key)

    print("Merged new 2 degree")
    adj_mat = make_adjacency(network)
    degrees = degree_counts(adj_mat)
    m_size,key_list,key_map,map_key = network_ids(network)

    high_nodes,edge_list = offending_nodes(network,degrees,map_key)
    make_clusters(network,high_nodes,edge_list)

    adj_mat = make_adjacency(network)
    degrees = degree_counts(adj_mat)
    m_size,key_list,key_map,map_key = network_ids(network)

    for i in range(len(degrees)):
        if degrees[i] < 4:
            adj_mat[i] = np.zeros(m_size)

    patterns = connected_components(adj_mat)

    print('Alter Patterns')
    removed_ids = alter_patterns(patterns,map_key,key_map,network)#,printing_flag="plot")
    adj_mat = make_adjacency(network)
    degrees = degree_counts(adj_mat)
    for i in removed_ids:
        del network['Vessel_list'][i]
        del network['centerline'][i]

    adj_mat = make_adjacency(network)
    degrees = degree_counts(adj_mat)
    m_size,key_list,key_map,map_key = network_ids(network)
    print('Merge 2 Degree')
    merge_2_degree(network,degrees,map_key)
    adj_mat = make_adjacency(network)
    degrees = degree_counts(adj_mat)
    m_size,key_list,key_map,map_key = network_ids(network)

    print('Remove Small Vessels')
    remove_small_vessels(network,map_key,degrees,vessel_cutoff=3)
    adj_mat = make_adjacency(network)
    degrees = degree_counts(adj_mat)
    m_size,key_list,key_map,map_key = network_ids(network)
    print('Merge 2 Degree')
    merge_2_degree(network,degrees,map_key)
    adj_mat = make_adjacency(network)
    degrees = degree_counts(adj_mat)
    m_size,key_list,key_map,map_key = network_ids(network)

#    network_status(network,patterns,f"{fname}/{fname}_network_status.txt")

    with open(f"{fname}_final_network.json",'w') as f:
        json.dump(network,f)

if __name__=="__main__":
    main()