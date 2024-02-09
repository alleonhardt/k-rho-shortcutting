use crate::weighted_graph::WeightedDirectedGraph;
use crate::dijkstra::{DijkstraFewestHops,VertexDistanceDescription};

pub fn k_rho_generic_worker<'a, I, S, F: Fn(u32,S, &Vec<Vec<u32>>, &DijkstraFewestHops)>(graph: &WeightedDirectedGraph, node_iterator: I, rho: u32, shortcut_set: S, apply: F)
where
    I: Iterator<Item=u32>, S: Copy {
    let mut dijk = DijkstraFewestHops::new(graph);

    for source in node_iterator {
        dijk.run(source);
        let sorted = dijk.get_nodes_sorted_by_distance();
        let mut children_tree: Vec<Vec<u32>> = vec![Vec::new(); graph.num_nodes() as usize];
        for VertexDistanceDescription{vertex,..} in &sorted[1..rho as usize+1] {
            let parent = dijk.parent(*vertex as u32);
            children_tree[parent as usize].push(*vertex as u32);
        }

        if children_tree[source as usize].len() == 0 {
            panic!("Graph is not connected!");
        }

	apply(source, shortcut_set, &children_tree, &dijk);
    }
}

#[derive(Copy,Clone,Debug)]
pub struct EdgeFrequencyDescription {
    pub num_times: u32,
    pub distinct_vertices: u32,
    pub fractional_distinct_vertices: f64,
    pub distance: u32
}
