use dashmap::DashMap;
use rand::thread_rng;
use rand::seq::SliceRandom;

use crate::weighted_graph::WeightedDirectedGraph;
use crate::util::k_rho_generic_worker;
use crate::Edge;

pub fn k_rho_blelloch_single_node(source: u32, k: u32, children: &Vec<Vec<u32>>) -> Vec<u32> {
    	    let mut levels: Vec<Vec<u32>> = vec![vec![source]];
    	    let mut current_level = 0;
    	    loop {
        	    let mut next_level: Vec<u32> = Vec::new();
        	    for child in &levels[current_level] {
            	    	for cc in &children[*child as usize] {
                	    	next_level.push(*cc);
            	    	}
        	    }
        	    if next_level.is_empty() {
            	    	break;
        	    }
        	    current_level+=1;
        	    levels.push(next_level);
    	    }
    	    let mut score_map: std::collections::HashMap<u32,Vec<u32>> = std::collections::HashMap::new();
    	    for x in (1..current_level+1).rev() {
        	    for vert in &levels[x] {
            	    	let mut default_score = Vec::new();
            	    	for t in 0..k {
                	    	let mut score_vec = Vec::new();
                	    	if t == 0 {
                    	    		score_vec.push(*vert);
                        	    	for child in &children[*vert as usize] {
                            	    		let id = *child*(k+1)+1;
                            	    		score_vec.extend_from_slice(&score_map.get(&id).unwrap()[..]);
                        	    	}
                    	    		let id = *vert*(k+1)+k;
                        	    	default_score = score_vec.clone();
                        	    	score_map.insert(id,score_vec);
                	    	}
                	    	else {
                        	    	for child in &children[*vert as usize] {
                            	    		// (child,t+1)
                            	    		let id = *child*(k+1)+(t+1);
                                	    	// F(w,
                            	    		score_vec.extend(score_map.get(&id).unwrap().iter());
                        	    	}
                    	    		let id = *vert*(k+1)+t;
                    	    		if score_vec.len() > default_score.len() {
                                	    	score_map.insert(id,default_score.clone());
                    	    		}
                    	    		else {
                                	    	score_map.insert(id,score_vec);
                    	    		}
                	    	}
            	    	}
        	    }
    	    }

   	    let mut final_vec = Vec::new();
    	    for child in &children[source as usize] {
        	    let id = child*(k+1)+k;
        	    final_vec.extend_from_slice(&score_map.get(&id).unwrap()[1..]);
    	    }
    	    final_vec
}

pub fn k_rho_blelloch_worker<'a, I>(graph: &WeightedDirectedGraph, node_iterator: I, k: u32, rho: u32, shortcut_set: &DashMap<u64,u32>)
where
    I: Iterator<Item=u32>
{
	k_rho_generic_worker(graph,node_iterator, rho, shortcut_set, |source,set,children,dijk| {
                let list = k_rho_blelloch_single_node(source,k,children);
                for y in list {
                    set.insert(graph.edge_uuid(source, y), dijk.distance(y));
                }
	});
}

pub fn k_rho_blelloch_worker_bottom_up<'a, I>(graph: &WeightedDirectedGraph, node_iterator: I, k: u32, rho: u32, shortcut_set: &DashMap<u64,u32>)
where
    I: Iterator<Item=u32>
{
	k_rho_generic_worker(graph,node_iterator, rho, shortcut_set, |source,set,children,dijk| {
            let list = k_rho_blelloch_single_node(source,k,children);
            for y in list {
                set.insert(graph.edge_uuid(source, y), dijk.distance(y));
            }
	});
}


pub fn k_rho_blelloch_star_worker<'a, I>(graph: &WeightedDirectedGraph, node_iterator: I, k: u32, rho: u32, shortcut_set: &DashMap<u64,u32>)
where
    I: Iterator<Item=u32>
{

    k_rho_generic_worker(graph,node_iterator, rho, shortcut_set, |source,set,children, dijkstra| {
            let list = k_rho_blelloch_single_node(source,k,children);
            for y in list {
                set.insert(graph.edge_uuid(source, y), dijkstra.distance(y));
            }
    });
}

pub fn k_rho_blelloch_bot_par(graph: &WeightedDirectedGraph, k: u32, rho: u32, num_threads: usize) -> DashMap<u64,u32> {
    let final_set = DashMap::<u64,u32>::new();
    let work_size = graph.num_nodes() as usize/num_threads as usize;


    std::thread::scope(|s| {
        for i in 0..num_threads {
            let borrowed_set = &final_set;
            s.spawn(move || {
                if i != num_threads-1 {
                	k_rho_blelloch_worker_bottom_up(graph,(i*work_size) as u32..((i+1)*work_size) as u32,k,rho,borrowed_set);
                }
                else {
                	k_rho_blelloch_worker_bottom_up(graph,(i*work_size) as u32..graph.num_nodes(),k,rho,borrowed_set);
                }
            });
        }
    });

    final_set
}

pub fn k_rho_blelloch_par(graph: &WeightedDirectedGraph, k: u32, rho: u32, num_threads: usize) -> DashMap<u64,u32> {
    let final_set = DashMap::<u64,u32>::new();
    let work_size = graph.num_nodes() as usize/num_threads as usize;


    std::thread::scope(|s| {
        for i in 0..num_threads {
            let borrowed_set = &final_set;
            s.spawn(move || {
                if i != num_threads-1 {
                	k_rho_blelloch_worker(graph,(i*work_size) as u32..((i+1)*work_size) as u32,k,rho,borrowed_set);
                }
                else {
                	k_rho_blelloch_worker(graph,(i*work_size) as u32..graph.num_nodes(),k,rho,borrowed_set);
                }
            });
        }
    });

    final_set
}



pub fn k_rho_blelloch_star_par(graph: &WeightedDirectedGraph, k: u32, rho: u32, phases: u32, num_threads: usize) -> DashMap<u64,u32> {
    let final_set = DashMap::<u64,u32>::new();
    let mut cloned_graph = graph.clone();
    let mut all_nodes: Vec<u32> = (0..graph.num_nodes()).collect();

    all_nodes.shuffle(&mut thread_rng());
    let phase_work_size = graph.num_nodes()/phases;


    for ph in 0..phases {
        let current_nodes = if ph != phases-1 {
            &all_nodes[(phase_work_size*ph) as usize..(phase_work_size*(ph+1)) as usize]
        }
        else {
            &all_nodes[(phase_work_size*ph) as usize..]
        };
        let borrowed = &cloned_graph;
        std::thread::scope(|s| {
            let work_size = current_nodes.len()/num_threads;
            for i in 0..num_threads {
                let borrowed_set = &final_set;
                s.spawn(move || {
                    if i != num_threads-1 {
                    	k_rho_blelloch_star_worker(borrowed,current_nodes[(i*work_size)..((i+1)*work_size)].into_iter().copied(),k,rho,borrowed_set);
                    }
                    else {
                    	k_rho_blelloch_star_worker(borrowed,current_nodes[(i*work_size)..].iter().copied(),k,rho,borrowed_set);
                    }
                });
            }
        });

        for x in final_set.iter() {
            let (source,target) = cloned_graph.uuid_to_edge(*x.key());

            if !cloned_graph.has_edge(source as u32,target as u32) {
                cloned_graph.add_edge(source as u32, Edge{ weight: *x.value(), target: target as u32});
            }
        }
    }

    final_set
}




pub fn k_rho_blelloch_greedy_single_node(source: u32, t: u32, k: u32, children_tree: &Vec<Vec<u32>>) -> Vec<u32> {
    if t == k+1 {
        let mut others = k_rho_blelloch_greedy_single_node(source,1,k,children_tree);
        others.push(source);
        return others;
    }
    else {
        if children_tree[source as usize].len() == 0 {
            return Vec::new();
        }
        let mut second_list = Vec::new();
        for x in &children_tree[source as usize] {
            let mut result = k_rho_blelloch_greedy_single_node(*x,t+1,k,children_tree);
            second_list.append(&mut result);
        }
        return second_list;
    }
}

pub fn k_rho_blelloch_greedy_worker<'a, I>(graph: &WeightedDirectedGraph, node_iterator: I, k: u32, rho: u32, shortcut_set: &DashMap<u64,u32>)
where
    I: Iterator<Item=u32>
{
    k_rho_generic_worker(graph,node_iterator, rho, shortcut_set, |source,set,children, dijk| {
            for x in &children[source as usize] {
                let list = k_rho_blelloch_greedy_single_node(*x,1,k,children);
                for y in list {
                    set.insert(graph.edge_uuid(source, y), dijk.distance(y));
                }
            }
    });
}

pub fn k_rho_blelloch_greedy_par(graph: &WeightedDirectedGraph, k: u32, rho: u32, num_threads: usize) -> DashMap<u64,u32> {
    let final_set = DashMap::<u64,u32>::new();
    let work_size = graph.num_nodes() as usize/num_threads as usize;


    std::thread::scope(|s| {
        for i in 0..num_threads {
            let borrowed_set = &final_set;
            s.spawn(move || {
                if i != num_threads-1 {
                	k_rho_blelloch_greedy_worker(graph,(i*work_size) as u32..((i+1)*work_size) as u32,k,rho,borrowed_set);
                }
                else {
                	k_rho_blelloch_greedy_worker(graph,(i*work_size) as u32..graph.num_nodes(),k,rho, borrowed_set);
                }
            });
        }
    });

    final_set
}
