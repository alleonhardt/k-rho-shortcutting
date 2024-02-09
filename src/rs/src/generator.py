import networkit as nk
import sys

def networkit_graph_to_string(g: nk.Graph):
    g = nk.components.ConnectedComponents.extractLargestConnectedComponent(g,True)
    print("{}".format(g.numberOfNodes()))

    for (u,v) in g.iterEdges():
        print("{} {}".format(u,v))

if len(sys.argv) < 4:
    print("Missing arguments")
    sys.exit(-1)

generator = sys.argv[1]
n = int(sys.argv[2])
average_number_of_edges = int(sys.argv[3])

if generator == 'gilbert':
    p = average_number_of_edges/n
    g = nk.generators.ErdosRenyiGenerator(n, p, directed=False).generate()
elif generator == 'hyperbolic':
    g = nk.generators.HyperbolicGenerator(n, k=average_number_of_edges).generate()
elif generator == 'powerlaw':
    g = nk.generators.HyperbolicGenerator(n, k=average_number_of_edges).generate()
    switcher = nk.randomization.EdgeSwitching(g,numberOfSwapsPerEdge=100.0)
    switcher.run()
    g = switcher.getGraph()
else:
    print("Unkown graph generator!")
    sys.exit(-1)

g = nk.components.ConnectedComponents.extractLargestConnectedComponent(g,True)
print("{}".format(g.numberOfNodes()))

for (u,v) in g.iterEdges():
    print("{} {}".format(u,v))
