use std::cmp::Ordering;
use crate::weighted_graph::{WeightedDirectedGraph};

#[derive(Copy, Clone, Eq, PartialEq)]
pub struct DijkstraDistState {
    tenative_distance: u32,
    vertex: u32,
}

impl Ord for DijkstraDistState {
    fn cmp(&self, other: &Self) -> Ordering {
        other.tenative_distance.cmp(&self.tenative_distance)
            .then_with(|| self.vertex.cmp(&other.vertex))
    }
}

impl PartialOrd for DijkstraDistState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Clone)]
struct DijkstraParentEntry {
    dist_from_root: u32,
    parent: u32,
}

pub struct VertexDistanceDescription {
    pub vertex: u32,
    pub distance: u32
}


pub struct DijkstraFewestHops<'a> {
    source: Option<u32>,
    graph: &'a WeightedDirectedGraph,
    distances: Vec<u32>,
    parent_array: Vec<DijkstraParentEntry>,
}

impl<'a> DijkstraFewestHops<'a> {
    pub fn new(graph: &'a WeightedDirectedGraph) -> Self {
        let distances: Vec<u32> = vec![u32::MAX; graph.num_nodes() as usize];

        DijkstraFewestHops {
            source: None,
            graph,
            distances,
            parent_array: vec![DijkstraParentEntry{ dist_from_root: u32::MAX, parent: 0}; graph.num_nodes() as usize]
        }
    }

    pub fn get_nodes_sorted_by_distance(&self) -> Vec<VertexDistanceDescription> {
        let mut sorted: Vec<VertexDistanceDescription> = self.distances.iter().enumerate().map(|(a,b)| VertexDistanceDescription{vertex: a as u32, distance: *b}).collect();
        sorted.sort_by(|a,b| a.distance.cmp(&b.distance));
        sorted
    }

    pub fn is_k_rho_ball(&self, k: u32, rho: u32) -> bool {
        self.get_nodes_sorted_by_distance()[1..rho as usize+1].iter().all(|VertexDistanceDescription{vertex,..}| self.parent_array[*vertex as usize].dist_from_root <= k)
    }

    pub fn get_path(&self, mut other_node: u32) -> Vec<u32> {
        let mut entry = &self.parent_array[other_node as usize];
        let mut path = vec![0;entry.dist_from_root as usize+1];
        while entry.dist_from_root != 0 {
            path[entry.dist_from_root as usize] = other_node;
            other_node = entry.parent;
            entry = &self.parent_array[other_node as usize];
        }
        path[0] = entry.parent;
        path
    }

    #[allow(dead_code)]
    pub fn is_ancestor_of(&self, mut node: u32, potential_ancestor: u32) -> bool {
        let mut entry = &self.parent_array[node as usize];
        while entry.dist_from_root != 0 {
            node = entry.parent;
            if node == potential_ancestor {
                return true;
            }
            entry = &self.parent_array[node as usize];
        }
        false
    }

    pub fn parent(&self, other_node: u32) -> u32 {
        self.parent_array[other_node as usize].parent
    }

    pub fn distance(&self, other_node: u32) -> u32 {
        return self.distances[other_node as usize];
    }

    pub fn hop_distance_from_root(&self, other_node: u32) -> u32 {
        self.parent_array[other_node as usize].dist_from_root
    }

    pub fn run(&mut self, source: u32) {
        if self.source.is_some() {
            self.distances.fill(u32::MAX);
            self.parent_array.fill(DijkstraParentEntry{ dist_from_root: u32::MAX, parent: 0});
        }
        self.source = Some(source);
        let mut pq = std::collections::BinaryHeap::with_capacity(self.graph.num_edges() as usize);

        pq.push(DijkstraDistState { tenative_distance: 0, vertex: source });
        self.distances[source as usize] = 0;
        self.parent_array[source as usize] = DijkstraParentEntry { dist_from_root: 0, parent: source };

        while let Some(DijkstraDistState {tenative_distance, vertex}) = pq.pop() {
            if tenative_distance > self.distances[vertex as usize] {
                continue;
            }

            for edge in self.graph.neighbors(vertex) {
            	let next = DijkstraDistState { tenative_distance: tenative_distance + edge.weight, vertex: edge.target };

                if next.tenative_distance < self.distances[next.vertex as usize] {
                    pq.push(next);
                    self.distances[next.vertex as usize] = next.tenative_distance;
                    self.parent_array[next.vertex as usize] = DijkstraParentEntry{ dist_from_root: self.parent_array[vertex as usize].dist_from_root+1, parent: vertex};
                }
                else if next.tenative_distance == self.distances[next.vertex as usize] && self.parent_array[vertex as usize].dist_from_root + 1 < self.parent_array[next.vertex as usize].dist_from_root {
                    pq.push(next);
                    self.parent_array[next.vertex as usize] = DijkstraParentEntry{ dist_from_root: self.parent_array[vertex as usize].dist_from_root+1, parent: vertex};
                }
            }
        }
    }
}
