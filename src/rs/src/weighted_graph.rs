use crate::dijkstra::DijkstraFewestHops;


#[derive(Clone)]
pub struct Edge {
    pub target: u32,
    pub weight: u32,
}

pub struct WeightedDirectedGraph {
    nodes: u32,
    edges: Vec<Vec<Edge>>,
    edge_index: std::collections::HashSet<u64>,
}

impl WeightedDirectedGraph {
    pub fn add_edge(&mut self, source: u32, edge: Edge) {
        debug_assert!(!self.has_edge(source,edge.target), "Edge exists already, inserted multi-edges!");
        self.edge_index.insert(self.edge_uuid(source,edge.target));
        self.edges[source as usize].push(edge);
    }

    pub fn clone(&self) -> WeightedDirectedGraph {
        WeightedDirectedGraph {
            nodes: self.nodes,
            edges: self.edges.clone(),
            edge_index: self.edge_index.clone()
        }
    }

    pub fn is_k_rho_graph(&self, k: u32, rho: u32, num_threads: usize) -> bool {
        let work_size = self.num_nodes() as usize/num_threads as usize;
        let is_k_rho_graph_flag = std::sync::atomic::AtomicBool::new(true);
        std::thread::scope(|s| {
            for i in 0..num_threads {
                let cloned = self;
                let atomic_borrow = &is_k_rho_graph_flag;
                s.spawn(move || {
                    if i != num_threads-1 {
                        let mut dijk = DijkstraFewestHops::new(cloned);
                    	for node in (i*work_size) as u32..((i+1)*work_size) as u32 {
                            dijk.run(node as u32);
                            if !dijk.is_k_rho_ball(k,rho) {
                                atomic_borrow.store(false, std::sync::atomic::Ordering::Relaxed);
                            }
                    	}
                    }
                    else {
                        let mut dijk = DijkstraFewestHops::new(cloned);
                    	for node in (i*work_size) as u32..cloned.num_nodes() {
                            dijk.run(node as u32);
                            if !dijk.is_k_rho_ball(k,rho) {
                                atomic_borrow.store(false, std::sync::atomic::Ordering::Relaxed);
                            }
                    	}
                    }
                });
            }
        });
        is_k_rho_graph_flag.load(std::sync::atomic::Ordering::Relaxed)
    }

    pub fn has_edge(&self, source: u32, target: u32) -> bool {
	self.edge_index.contains(&(self.edge_uuid(source,target)))
    }


    pub fn edge_uuid(&self, source: u32, target: u32) -> u64 {
        return source as u64*self.nodes as u64 + target as u64
    }

    pub fn uuid_to_edge(&self, uuid: u64) -> (u32,u32) {
        let target = uuid  % self.num_nodes() as u64;
        let source = uuid / self.num_nodes() as u64;
        (source as u32,target as u32)
    }

    pub fn num_nodes(&self) -> u32 {
        self.nodes
    }

    pub fn degree_of(&self, node: u32) -> usize {
        self.edges[node as usize].len()
    }

    pub fn num_edges(&self) -> u32 {
        self.edge_index.len() as u32
    }

    pub fn neighbors(&self, node: u32) -> &[Edge] {
        &self.edges[node as usize]
    }

    pub fn from_generator(input: &str) -> Result<WeightedDirectedGraph,()> {
        let mut lines = input.lines();
        let nodes = lines.next().unwrap().parse::<u32>().unwrap();

        let mut edges = Vec::with_capacity(nodes as usize);
        for _ in 0..nodes {
            edges.push(Vec::new());
        }

        let mut graph = WeightedDirectedGraph {
            nodes,
            edges,
            edge_index: std::collections::HashSet::new(),
        };

        for ln in lines {
            let mut iterspaces = ln.split(" ");
            let first = iterspaces.next().unwrap().parse::<usize>().unwrap();
            let second = iterspaces.next().unwrap().parse::<usize>().unwrap();
            graph.add_edge(first as u32, Edge { target: second as u32, weight: 1 });
            graph.add_edge(second as u32, Edge { target: first as u32, weight: 1 });
        }
        Ok(graph)
   }
}
