import numpy as np
import networkit as nk
from networkit import *
from mip import *
import math
import sys
import pandas as pd
try:
    import regex as re
except ImportError:
    import re
import copy
import time
import multiprocessing
import networkx as nx
import math
from subprocess import Popen, PIPE, STDOUT
import subprocess
import json

def get_git_revision_hash() -> str:
    return subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()


def merge_dicts(first, second):
    for key in second:
        if key in first:
            first[key].append(second[key])
        else:
            first[key] = [second[key]]

def networkit_graph_to_string(g: nk.Graph):
    g = nk.components.ConnectedComponents.extractLargestConnectedComponent(g,True)
    result_string = ''
    result_string += "{}\n".format(g.numberOfNodes())

    for (u,v) in g.iterEdges():
        result_string += "{} {}\n".format(u,v)

    return result_string


def dijkstra_algorithm(graph, start_node):
    unvisited_nodes = [n for n in graph.iterNodes()]
 
    shortest_path = {}

    paths = {}
    paths[start_node] = [start_node]
 
    max_value = sys.maxsize
    for node in unvisited_nodes:
        shortest_path[node] = max_value

    shortest_path[start_node] = 0
    
    while unvisited_nodes:
        current_min_node = None
        for node in unvisited_nodes:
            if current_min_node == None:
                current_min_node = node
            elif shortest_path[node] < shortest_path[current_min_node]:
                current_min_node = node
                
        neighbors = [nb for nb in graph.iterNeighbors(current_min_node)]
        for neighbor in neighbors:
            tentative_value = shortest_path[current_min_node] + graph.weight(current_min_node, neighbor)
            if tentative_value < shortest_path[neighbor]:
                shortest_path[neighbor] = tentative_value
                paths[neighbor] = paths[current_min_node]+[neighbor]
            elif tentative_value == shortest_path[neighbor]:
                new_path=paths[current_min_node]+[neighbor]

                if len(new_path) < len(paths[neighbor]):
                    paths[neighbor] = new_path
        unvisited_nodes.remove(current_min_node)
    
    return paths, shortest_path

def deepcopy_graph(g: nk.Graph):
    ng = nk.Graph(g.numberOfNodes(), weighted=True, directed=True)
    if g.isDirected():
        for (u,v,w) in g.iterEdgesWeights():
            ng.addEdge(u,v,w=w)
    else:
        for (u,v,w) in g.iterEdgesWeights():
            ng.addEdge(u,v,w=w)
            ng.addEdge(v,u,w=w)
    return ng

class ShortCut:
    def __init__(self, var, start, end):
        self.var = var
        self.start = start 
        self.end = end


def build_edge_var_name(start,end, u, v):
    return '{}-{}@{}-{}'.format(start,end,u,v)

def build_shortcut_rep_var_name(s,start,end):
    return build_shortcut_rep_var_name_from_string(s.name,start,end)

def build_shortcut_rep_var_name_from_string(string, start, end):
    return "{}-{}@".format(start,end)+string


def build_shortcut_name(u,v):
    return "shortcut-{}-{}".format(u,v)

def build_k_rho_finder(start: int, end: int, graph: nk.Graph, k: int, model: Model, splen: float, shortcuts: list[ShortCut], shortcuts_weight, apsp_exec):
    allvars = []
    weights = {}
    for u,v,w in graph.iterEdgesWeights():
        if w > splen:
            continue
        edge_var_name = build_edge_var_name(start,end,u,v)
        weights[edge_var_name] = w
        var = model.add_var(name=edge_var_name,var_type=BINARY)
        allvars.append(var)

    for s in shortcuts:
        if shortcuts_weight[s.var.name] > splen:
            continue
        if apsp_exec.getDistance(start,s.start)+apsp_exec.getDistance(s.start,s.end)+apsp_exec.getDistance(s.end,end) > splen:
            continue
        svn = build_shortcut_rep_var_name(s.var,start,end)
        var = model.add_var(name=svn,var_type=BINARY)
        model+= var <= s.var
        allvars.append(var)
        if svn in weights:
            raise Exception("Shortcut already in weights")
        weights[svn] = shortcuts_weight[s.var.name]


    flow_in = [[] for _ in range(graph.numberOfNodes())]
    flow_out = [[] for _ in range(graph.numberOfNodes())]

    for u in graph.iterNodes():
        for v in graph.iterNodes():
            if u == v:
                continue
            if graph.hasEdge(u,v):
                flow_out_name = build_edge_var_name(start,end,u,v)
                fout = model.var_by_name(flow_out_name)
                if fout is not None:
                    flow_out[u].append(fout)
            else:
                svn = build_shortcut_rep_var_name_from_string(build_shortcut_name(u,v), start, end)
                fout = model.var_by_name(svn)
                if fout is not None:
                    flow_out[u].append(fout)

            if graph.hasEdge(v,u):
                flow_in_name = build_edge_var_name(start,end,v,u)
                fin = model.var_by_name(flow_in_name)
                if fin is not None:
                    flow_in[u].append(fin)
            else:
                svn = build_shortcut_rep_var_name_from_string(build_shortcut_name(v,u), start, end)
                fin = model.var_by_name(svn)
                if fin is not None:
                    flow_in[u].append(fin)

    for s in range(graph.numberOfNodes()):
        if s == end:
            model+= -xsum(flow_out[s]) - xsum(flow_in[s])== -1
        elif s == start:
            model+= xsum(flow_out[s]) - xsum(flow_in[s]) == 1
        elif len(flow_out[s]) == 0 and len(flow_in[s]) == 0:
            continue
        else:
            model+= xsum(flow_out[s]) == xsum(flow_in[s])

    model+=xsum(v for v in allvars) <= k
    model+=xsum(v*weights[v.name] for v in allvars) == splen


def k_rho_ilp_solver(graph: nk.Graph, k: int, rho: int):
    if graph.isWeighted() == False:
        raise Exception("Graph is not weighted")
    elif graph.isDirected() == False:
        raise Exception("Graph is not directed")

    m = Model()

    APSP = nk.distance.APSP(graph)
    APSP.run()
    shortcuts = []
    shortcuts_weight = {}
    for x in range(graph.numberOfNodes()):
        for y in range(graph.numberOfNodes()):
            if x == y:
                continue
            if not graph.hasEdge(x,y):
                dist = APSP.getDistance(x,y)
                if dist == sys.float_info.max:
                    raise Exception("Distance is bigger than max path length, graph is not connected!")
                first_name = build_shortcut_name(x,y)
                shortcuts_weight[first_name] = dist
                var = m.add_var(name=first_name,var_type=BINARY)
                shortcuts.append(ShortCut(var,x,y))


    listoflist = APSP.getDistances()
    withid = []
    for lk in listoflist:
        withid.append([(i,k) for i,k in enumerate(lk)])
    for i in range(graph.numberOfNodes()):
        withid[i].pop(i)

    withid_sorted = [sorted(sublist, key=lambda val: val[1]) for sublist in withid]

    
    for i in range(graph.numberOfNodes()):
        dj = nk.distance.Dijkstra(graph,i,storePaths=True)
        dj.run()
        ids = withid_sorted[i][:rho]
        for id,_ in ids:
            if len(dj.getPath(id)) > 1:
                build_k_rho_finder(i,id,graph,k,m,APSP.getDistance(i,id),shortcuts,shortcuts_weight,APSP)

    # Minimize the number of shortcuts
    min_shorts = xsum(v.var for v in shortcuts)

    m.objective = minimize(min_shorts)
    m.threads = multiprocessing.cpu_count()

    status = m.optimize(max_seconds=1800)
    print(status)
    
    transformed_graph = deepcopy_graph(graph)
    prog = re.compile("shortcut-(\d+)-(\d+)")
    if status == OptimizationStatus.OPTIMAL:
        count = 0
        for v in m.vars:
            if v.x > 0:
                if v.name.startswith("shortcut"):
                    result = prog.match(v.name)
                    transformed_graph.addEdge(int(result.group(1)),int(result.group(2)), w=APSP.getDistance(int(result.group(1)),int(result.group(2))))
                    count += 1
        if is_k_rho_graph(transformed_graph,k,rho) == False:
            print("Algorithm produced wrong solution; OPT")
            sys.exit(-1)
        return (count,count)

    elif status == OptimizationStatus.FEASIBLE:
        count = 0
        for v in m.vars:
            if v.x > 0:
                if v.name.startswith("shortcut"):
                    count += 1
        return (-1,count)
    elif status == OptimizationStatus.NO_SOLUTION_FOUND:
        print("Error: There should be a solution for any graph")
        sys.exit(-1)
    else:
        return (-1,-1)



def is_k_rho_graph(graph: nk.Graph, k: int, rho: int):
    apsp = nk.distance.APSP(graph)
    apsp.run()

    for x in range(graph.numberOfNodes()):
        paths,dist = dijkstra_algorithm(graph,x)
        sorted_dist = sorted([(n,dist) for (n,dist) in dist.items()], key=lambda x: x[1])
        rho_nearest = sorted_dist[:(rho+1)]
        for (n,_) in rho_nearest:
            if len(paths[n]) == 1:
                continue
            if len(paths[n]) > k+1:
                return False

            total_weight = 0.0
            last = None
            for inner in paths[n]:
                if last == None:
                    last = inner
                else:
                    total_weight += graph.weight(last,inner)
                    last = inner
            if total_weight > apsp.getDistance(x,n):
                return False

    return True


def read_network_repository_graph(path, n):
    file = open("/data/networks/network-repository.com/"+path,'r')
    lines = file.readlines()

    graph = nk.Graph(n,directed=False,weighted=False) 

    for (id,line) in enumerate(lines):
        if id == 0:
            continue
        else:
            split = line.split(",")
            f = int(split[0])
            t = int(split[1])
            graph.addEdge(f,t)
    return graph






sys.setrecursionlimit(100000)

if len(sys.argv) < 4:
    if len(sys.argv) != 3 or sys.argv[1] != 'real_world':
        print("Usage: main.py [setup] [output-file-name] [generator]\n setup    --- small_ilp,large,real_world\n generator --- gilbert,hyperbolic,powerlaw,euclidean")
        sys.exit(-1)

n_vals = []
use_ilp = True
exact = True
if sys.argv[1] == 'small_ilp':
    n_vals = [20,25,30,35,40,45,50,55,60,65,70]
    k_vals = [2]
    rho_vals = [lambda x: x-1]
    exact = True
elif sys.argv[1] == 'large':
    use_ilp = False
    n_vals = [2000,4000,6000,8000]
    k_vals = [2,3,4,5]
    exact = False
    rho_vals = [lambda x: x-1, lambda x: int(x/2), lambda x: int(x/3), lambda x: int(x/4), lambda x: int(x/10), lambda x: int(math.sqrt(x))]
elif sys.argv[1] == 'real_world':
    data_g = pd.read_csv('/data/networks/network-repository.com/networks.csv')
    print(data_g[data_g.nodes < 10000]) 

    data = {}
    k=3
    for _,row in data_g[data_g.nodes < 10000].iterrows():
        g = read_network_repository_graph(row['file'],row['nodes'])


        largest_connected_component = nk.components.ConnectedComponents.extractLargestConnectedComponent(g,True)
        print(f"Processing {row['file']}")
        rho = largest_connected_component.numberOfNodes()-1
        p = Popen(['src/rs/target/release/rs', str(k), str(rho)], stdout=PIPE, stdin=PIPE, stderr=PIPE, text=True)
        (stdout_data,stderr_data) = p.communicate(input=networkit_graph_to_string(largest_connected_component))

        if p.returncode != 0:
            print(stderr_data)
            sys.exit(-1)

        ret_data = json.loads(stdout_data)
        if ret_data["build_hash"] != get_git_revision_hash():
            print("Error rust binary git hash does not match python hash, please recompile!",file=sys.stderr)
            sys.exit(-1)
        merge_dicts(data,ret_data)
        merge_dicts(data, {"n": largest_connected_component.numberOfNodes(), "m": largest_connected_component.numberOfEdges()*2, "rho": rho, "k": k, "file": row["file"]})

        df = pd.DataFrame(data)
        df.to_csv(sys.argv[2], index=False)
    sys.exit(0)
else:
    raise Exception("Unimplemented")



data = {}
for k in k_vals:
    for rho_fn in rho_vals:
        for n in n_vals:
            count = 0
            while count < 50:
                    p = (3/n)
                    if sys.argv[3] == 'gilbert':
                        g = nk.generators.ErdosRenyiGenerator(n, p, directed=False).generate()
                    elif sys.argv[3] == 'hyperbolic':
                        g = nk.generators.HyperbolicGenerator(n, k=3).generate()
                    elif sys.argv[3] == 'powerlaw':
                        g = nk.generators.HyperbolicGenerator(n, k=3).generate()
                        switcher = nk.randomization.EdgeSwitching(g,numberOfSwapsPerEdge=100.0)
                        switcher.run()
                        g = switcher.getGraph()
                    else:
                        print("Unknown generator!")
                        sys.exit(-1)


                    largest_connected_component = nk.components.ConnectedComponents.extractLargestConnectedComponent(g,True)
                    if exact and largest_connected_component.numberOfNodes() < n:
                        continue

                    rho = rho_fn(largest_connected_component.numberOfNodes())


                    if use_ilp:
                        weighted = deepcopy_graph(nk.graphtools.toWeighted(largest_connected_component))

                        time_ilp = time.monotonic_ns()
                        (ilp,feasable) = k_rho_ilp_solver(weighted,k,rho)
                        time_ilp_ms = (time.monotonic_ns()-time_ilp)/(1000*1000)

                        merge_dicts(data,{"size_ilp": ilp, "speed_ilp_ms": time_ilp_ms, "size_ilp_feasable": feasable})

                    p = Popen(['src/rs/target/release/rs', str(k), str(rho)], stdout=PIPE, stdin=PIPE, stderr=PIPE, text=True)
                    (stdout_data,stderr_data) = p.communicate(input=networkit_graph_to_string(largest_connected_component))

                    if p.returncode != 0:
                        print(stderr_data)
                        sys.exit(-1)

                    ret_data = json.loads(stdout_data)
                    if ret_data["build_hash"] != get_git_revision_hash():
                        print("Error rust binary git hash does not match python hash, please recompile!",file=sys.stderr)
                        sys.exit(-1)
                    merge_dicts(data,ret_data)
                    merge_dicts(data, {"n": largest_connected_component.numberOfNodes(), "m": largest_connected_component.numberOfEdges()*2, "rho": rho, "k": k, "generator": sys.argv[3]})

                    df = pd.DataFrame(data)
                    df.to_csv(sys.argv[2], index=False)
                    count+=1
