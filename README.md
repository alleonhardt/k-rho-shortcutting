# $(k,\rho)$-shortcutting algorithms
This repository contains the code and accompanied scripts for the analysis and experimental part of the paper "Insights into $(k,\rho)$-shortcutting algorithms".

The code is split into several parts, `src/main.py` which contains the optimal ILP solver and `src/rs/*` which contains the implementation of all heuristics.
Normally the python code is used to interact with the rust code as the rust code provides functionality to read a graph from the standard input.
For testing purposes the rust binary also implements a routine to generate a single graph with either an hyperbolic random graph generator, the edge switchting markov chain or the gilbert random graph generator.
To use this functionality python should be installed alongside some other modules.

```
# Construct a new virtual environment to install our packages in
python -m venv analyse

# install the required libraries
analyse/bin/pip install -r requirements.txt

# Build the rust binary
cd src/rs && cargo build --release

# Return to root directory
cd ../..

# Generate large random graphs and solve them with all heuristics where the final argument is either hyperbolic, gilbert or powerlaw for the respective random graph model investigated
analyse/bin/python src/main.py large output.csv hyperbolic 

# This invokes the ILP solver in addition to all heuristics on very small graphs. The final argument is again the random graph model used, hence it is has to be one of the following hyperbolic, gilbert or powerlaw
analyse/bin/python src/main.py small_ilp output.csv gilbert

# This invokes all heuristics on graphs from the network repository with less than 10000 nodes, observe that for that to work there has to be a copy of the data at /data/networks/network-repository.com
# Furthermore there needs to be a CSV file /data/networks/network-repository.com/networks.csv tracking at least the number of nodes and relative file location within the data
analyse/bin/python src/main.py real_world out.csv


# Generates a random graph of said model [1] with n nodes [2] with average degree a [3] and running all heuristics with parameters k [4] and rho [5]
# Observe, if rho exceeds n it is reset to the highest value possible n-1
cd src/rs && cargo run --release -- hyperbolic 3000 4 5 3000

# Runs all heuristics with parameters k [1] and rho [2]. Expects that a connected undirected unweighted graph is supplied via stdin in the following form:
# n
# 2 3
# 4 5
# 1 2
# The aforementioned input results in a graph with n nodes and edges {{2,3},{4,5},{1,2}}.

cargo run --release -- 4 3000
```
