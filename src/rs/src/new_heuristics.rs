use dashmap::{DashMap,DashSet};
use dashmap::mapref::entry::Entry;
use std::sync::Barrier;


use crate::dijkstra::{DijkstraFewestHops};
use crate::weighted_graph::WeightedDirectedGraph;
use crate::Edge;
use crate::{MinHash,MinHashSignature, EdgeFrequencyDescription, VertexDistanceDescription};
use crate::blelloch::{k_rho_blelloch_par,k_rho_blelloch_star_par,k_rho_blelloch_single_node};
use crate::util::k_rho_generic_worker;

pub fn k_rho_leonhardt_worker_dp<'a, I>(graph: &WeightedDirectedGraph, node_iterator: I, k: u32, rho: u32, optimal_scuts: &DashMap<u64,(u32,u32)>)
where
    I: Iterator<Item=u32>
{
	k_rho_generic_worker(graph,node_iterator, rho, 0, |source,_,children,dijk| {
    	    let mut extended_list = k_rho_blelloch_single_node(source,k,children);
                if extended_list.len() > 0 {
                    let mut shortcut_ancestor_succ_depth = std::collections::HashMap::<u32,u32>::new();
                    extended_list.insert(0,source);
                    extended_list.iter().for_each(|x| { shortcut_ancestor_succ_depth.insert(*x, 0);});

		    // Vertex and closest ancestor shortcut
                    let mut current_frontier: Vec<(u32,u32)> = vec![(source,source)];
                    let mut next_frontier:Vec<(u32,u32)> = Vec::new();
                    let mut depth = 1;

                    while current_frontier.len() > 0 {
                        let mut kill: bool = false;
                        for (ancestor,node) in current_frontier.iter() {
                            	for child in &children[*node as usize] {
                                	if shortcut_ancestor_succ_depth.contains_key(child) {
                                    		// No depth control as has direct connection from source.
                                    		next_frontier.push((*child,*child));
                                	}
                                	else {
                                    		let entry = shortcut_ancestor_succ_depth.get_mut(ancestor).unwrap();
                                    		*entry = (*entry).max(depth-dijk.hop_distance_from_root(*ancestor));
                                    		if *entry > k {
                                        		println!("NA: depth {} and ancestor {} dist {}",depth, *ancestor, dijk.distance(*child));
                                        		kill = true;
                                    		}
                                    		next_frontier.push((*ancestor,*child));
                                	}

                            	}
                        }
                        depth += 1;
                        if kill {
                            panic!("Should never happen!");
                        }
                        current_frontier = next_frontier;
                        next_frontier = Vec::new();
                    }

		    let mut double_counting_protector = std::collections::HashSet::<u64>::new();
                    for x in extended_list.iter() {
        	    	let key = graph.edge_uuid(source,*x);
        	    	optimal_scuts.entry(key).and_modify(|e| {e.0 += 1;}).or_insert((1,dijk.distance(*x)));
                    }

                    for optimal_short in shortcut_ancestor_succ_depth.iter() {
                        if *optimal_short.1 < k && *optimal_short.0 != source {
                            let path  = dijk.get_path(*optimal_short.0);
                            let filtered_path: Vec<(usize,u32)> = path.iter().enumerate().map(|(u,x)| (u,*x)).filter(|(_,vert)| shortcut_ancestor_succ_depth.contains_key(vert)).collect();

                            let logrho = rho.ilog2() as usize;
                            let final_logrho = if filtered_path.len() > logrho {&filtered_path[filtered_path.len()-logrho..]} else {&filtered_path};

			    for dist_from_source in 1..(k-optimal_short.1) {
    			    	for (index,vert) in final_logrho.iter() {
        			    		if *vert == *optimal_short.0 {
                			    		continue;
        			    		}
        			    		else if *vert == source {
                			    		let start = path[index+dist_from_source as usize];
                			    		if start == dijk.parent(*optimal_short.0) {
                    			    			continue;
                			    		}
                        			    	let key = graph.edge_uuid(start,*optimal_short.0);

                            			    	let inc = if double_counting_protector.contains(&key) {
                                			    	0
                            			    	}
                            			    	else {
                                			    	double_counting_protector.insert(key);
                                			    	1
                            			    	};
                        			    	optimal_scuts.entry(key).and_modify(|e| {e.0 += inc;}).or_insert((1,dijk.distance(*optimal_short.0)-dijk.distance(start)));
        			    		}
        			    		else {
                			    		let start = path[index+(dist_from_source as usize-1)];
                			    		if start == dijk.parent(*optimal_short.0) {
                    			    			continue;
                			    		}
                        			    	let key = graph.edge_uuid(start,*optimal_short.0);

                            			    	let inc = if double_counting_protector.contains(&key) {
                                			    	0
                            			    	}
                            			    	else {
                                			    	double_counting_protector.insert(key);
                                			    	1
                            			    	};
                        			    	optimal_scuts.entry(key).and_modify(|e| {e.0 += inc;}).or_insert((1,dijk.distance(*optimal_short.0)-dijk.distance(start)));
    			    		}
    			    	}
        		    	}
                        }
                    }
                }
	});
}


pub fn k_rho_leonhardt_worker_dp_minhash<'a, I>(graph: &WeightedDirectedGraph, node_iterator: I, k: u32, rho: u32, optimal_scuts: &DashMap<u64,(u32,u32)>,fingerprint: &std::sync::Arc<DashMap<u64,MinHashSignature>>, minhash: &MinHash)
where
    I: Iterator<Item=u32>
{
	k_rho_generic_worker(graph,node_iterator, rho, 0, |source,_,children,dijk| {
    	    let mut extended_list = k_rho_blelloch_single_node(source,k,children);
                if extended_list.len() > 0 {
                    let mut shortcut_ancestor_succ_depth = std::collections::HashMap::<u32,u32>::new();
                    extended_list.insert(0,source);
                    extended_list.iter().for_each(|x| { shortcut_ancestor_succ_depth.insert(*x, 0);});

		    // Vertex and closest ancestor shortcut
                    let mut current_frontier: Vec<(u32,u32)> = vec![(source,source)];
                    let mut next_frontier:Vec<(u32,u32)> = Vec::new();
                    let mut depth = 1;

                    while current_frontier.len() > 0 {
                        let mut kill: bool = false;
                        for (ancestor,node) in current_frontier.iter() {
                            	for child in &children[*node as usize] {
                                	if shortcut_ancestor_succ_depth.contains_key(child) {
                                    		// No depth control as has direct connection from source.
                                    		next_frontier.push((*child,*child));
                                	}
                                	else {
                                    		let entry = shortcut_ancestor_succ_depth.get_mut(ancestor).unwrap();
                                    		*entry = (*entry).max(depth-dijk.hop_distance_from_root(*ancestor));
                                    		if *entry > k {
                                        		println!("NA: depth {} and ancestor {} dist {}",depth, *ancestor, dijk.distance(*child));
                                        		kill = true;
                                    		}
                                    		next_frontier.push((*ancestor,*child));
                                	}

                            	}
                        }
                        depth += 1;
                        if kill {
                            panic!("Should never happen!");
                        }
                        current_frontier = next_frontier;
                        next_frontier = Vec::new();
                    }

		    let mut double_counting_protector = std::collections::HashSet::<u64>::new();
                    for x in extended_list.iter() {
        	    	let key = graph.edge_uuid(source,*x);
        	    	optimal_scuts.entry(key).and_modify(|e| {e.0 += 1;}).or_insert((1,dijk.distance(*x)));
        		fingerprint.entry(key).and_modify(|e| minhash.add(key,e)).or_insert(minhash.generate_signature(key));
                    }

                    for optimal_short in shortcut_ancestor_succ_depth.iter() {
                        if *optimal_short.1 < k && *optimal_short.0 != source {
        	            let outer_key = graph.edge_uuid(source, *optimal_short.0);
                            let path  = dijk.get_path(*optimal_short.0);
                            let filtered_path: Vec<(usize,u32)> = path.iter().enumerate().map(|(u,x)| (u,*x)).filter(|(_,vert)| shortcut_ancestor_succ_depth.contains_key(vert)).collect();

                            let logrho = rho.ilog2() as usize;
                            let final_logrho = if filtered_path.len() > logrho {&filtered_path[filtered_path.len()-logrho..]} else {&filtered_path};

			    for dist_from_source in 1..(k-optimal_short.1) {
    			    	for (index,vert) in final_logrho.iter() {
        			    		if *vert == *optimal_short.0 {
                			    		continue;
        			    		}
        			    		else if *vert == source {
                			    		let start = path[index+dist_from_source as usize];
                			    		if start == dijk.parent(*optimal_short.0) {
                    			    			continue;
                			    		}
                        			    	let key = graph.edge_uuid(start,*optimal_short.0);

                            			    	let inc = if double_counting_protector.contains(&key) {
                                			    	0
                            			    	}
                            			    	else {
                                			    	double_counting_protector.insert(key);
                                			    	1
                            			    	};
                        			    	optimal_scuts.entry(key).and_modify(|e| {e.0 += inc;}).or_insert((1,dijk.distance(*optimal_short.0)-dijk.distance(start)));
                                        		fingerprint.entry(key).and_modify(|e| minhash.add(outer_key,e)).or_insert(minhash.generate_signature(outer_key));
        			    		}
        			    		else {
                			    		let start = path[index+(dist_from_source as usize-1)];
                			    		if start == dijk.parent(*optimal_short.0) {
                    			    			continue;
                			    		}
                        			    	let key = graph.edge_uuid(start,*optimal_short.0);

                            			    	let inc = if double_counting_protector.contains(&key) {
                                			    	0
                            			    	}
                            			    	else {
                                			    	double_counting_protector.insert(key);
                                			    	1
                            			    	};
                        			    	optimal_scuts.entry(key).and_modify(|e| {e.0 += inc;}).or_insert((1,dijk.distance(*optimal_short.0)-dijk.distance(start)));
                                        		fingerprint.entry(key).and_modify(|e| minhash.add(outer_key,e)).or_insert(minhash.generate_signature(outer_key));
    			    		}
    			    	}
        		    	}
                        }
                    }
                }
	});
}


pub struct SyncStats {
    fraction_distinct_vertices_sum: f64,
    fraction_distinct_vertices_squared_sum: f64
}

impl SyncStats {
    pub fn new() -> SyncStats {
        SyncStats {
            fraction_distinct_vertices_sum: 0.0,
            fraction_distinct_vertices_squared_sum: 0.0
        }
    }

    pub fn update(&mut self, added: f64, before: f64) {
        self.fraction_distinct_vertices_sum += added;
        self.fraction_distinct_vertices_squared_sum += -(before*before) + ((before+added)*(before+added));
    }

    pub fn calculate_mean_stddev(&self, n_elem: usize) -> (f64,f64) {
        let avg = self.fraction_distinct_vertices_sum/n_elem as f64;
        let standard_deviation = ((self.fraction_distinct_vertices_squared_sum-2.0*self.fraction_distinct_vertices_sum*avg+n_elem as f64*avg*avg)/n_elem as f64).sqrt();
        (avg,standard_deviation)
    }
}



pub fn k_rho_leonhardt_worker_base_new_ud<'a, I>(graph: &WeightedDirectedGraph, node_iterator: I, k: u32, rho: u32, shortcut_set: &DashMap<u64,EdgeFrequencyDescription>, counter: &std::sync::atomic::AtomicU32, sync: &std::sync::RwLock<SyncStats>)
where
    I: Iterator<Item=u32>
{
    let mut dijk = DijkstraFewestHops::new(graph);

    let mut local_counter: u32 = 0;
    let mut local_sync: SyncStats = SyncStats::new();

    for source in node_iterator {
        let mut covered_pairs = std::collections::HashSet::<u64>::new();
        dijk.run(source);

        let sorted = dijk.get_nodes_sorted_by_distance();

        let mut unique_children = std::collections::HashMap::<u32,(u32,u32)>::new();
        let mut children_tree: Vec<Vec<u32>> = vec![Vec::new(); graph.num_nodes() as usize];
        for VertexDistanceDescription{vertex,..} in &sorted[1..rho as usize+1] {
            let parent = dijk.parent(*vertex as u32);
            children_tree[parent as usize].push(*vertex as u32);
        }

        let mut current_frontier: Vec<(u32,u32)> = vec![(source,source)];
        let mut next_frontier:Vec<(u32,u32)> = Vec::new();

        while current_frontier.len() > 0 {
            for (parent,node) in current_frontier.iter() {
                if children_tree[*node as usize].len() == 0 {
                    if !unique_children.contains_key(parent) {
                        unique_children.insert(*parent, (*node,*node));
                    }
                    else {
                        let res = unique_children.get_mut(parent).unwrap();
                        if res.0 == res.1 {
                            res.1 = *node;
                        }
                    }
                    continue;
                }
        	for child in &children_tree[*node as usize] {
        		next_frontier.push((*node,*child));
            	}
            }
            current_frontier = next_frontier;
            next_frontier = Vec::new();
        }


        let total: usize = sorted[1..rho as usize+1].iter().filter(|x| {
            let path = dijk.get_path(x.vertex);

            if children_tree[x.vertex as usize].len() == 0 {

                let ent = unique_children.get(&path[path.len()-2]).unwrap();
                if ent.0 != x.vertex && ent.1 != x.vertex {
                    return false;
                }
            }
            return true;
            
                }).map(|x| if dijk.get_path(x.vertex).len() <= k as usize+1 { 0 } else { 1 }).sum();
        for VertexDistanceDescription{vertex,..} in &sorted[1..rho as usize+1] {
            let path = dijk.get_path(*vertex as u32);

            if path.len()-1 <= k as usize {
                continue;
            }

            // Skip not unique children
            if children_tree[*vertex as usize].len() == 0 {
                let ent = unique_children.get(&path[path.len()-2]).unwrap();
                if (ent.0 != *vertex) && (ent.1 != *vertex) {
                    continue;
                }
            }

            local_counter+=1;

            for (index,node) in (&path[0..k as usize]).iter().enumerate() {
                let start_index = path.len()-(k as usize-index);
                if index+1 == start_index {
                    panic!("Should never happen!");
                }

                for pair in &path[start_index as usize..] {
                    let key = graph.edge_uuid(*node,*pair);
                    let distinct = if !covered_pairs.contains(&key) {
                        covered_pairs.insert(key);
			1
                    }
                    else {
                        0
                    };

                    let item = shortcut_set.entry(key);
                    match item {
                        Entry::Occupied(mut entry_lck) => {
                            let entry = entry_lck.get_mut();
                            entry.num_times+=1;
                            local_sync.update(1.0/total as f64, entry.fractional_distinct_vertices);
                            entry.fractional_distinct_vertices += 1.0/total as f64;
                            entry.distinct_vertices+=distinct;
                        }
                        Entry::Vacant(entry) => {
                            entry.insert(EdgeFrequencyDescription{num_times: 1, fractional_distinct_vertices: 1.0/total as f64, distinct_vertices: 1, distance: dijk.distance(*pair)-dijk.distance(*node)});
                            local_sync.update(1.0/total as f64, 0.0);
                        }
                    }
                }
            }
        }
    }
    counter.fetch_add(local_counter,std::sync::atomic::Ordering::Relaxed);
    let mut lck = sync.write().unwrap();
    lck.fraction_distinct_vertices_sum += local_sync.fraction_distinct_vertices_sum;
    lck.fraction_distinct_vertices_squared_sum += local_sync.fraction_distinct_vertices_squared_sum;
}


pub fn k_rho_leonhardt_base_ud_par(graph: &WeightedDirectedGraph, k: u32, rho: u32, num_threads: usize) -> DashMap<u64,u32> {
    let work_size = graph.num_nodes() as usize/num_threads as usize;
    let cumulative_shortcut_hits = DashMap::<u64,EdgeFrequencyDescription>::new();

    let atomic_counter = std::sync::atomic::AtomicU32::new(0);
    let sync = std::sync::RwLock::new(SyncStats::new());
    let barrier = Barrier::new(num_threads);

    let mut final_set = DashMap::<u64,u32>::new();


    std::thread::scope(|s| {
        for i in 0..num_threads {
            let borrowed_set = &cumulative_shortcut_hits;
            let borrowed_counter = &atomic_counter;
            let borrowed_sync = &sync;
            let borrowed_barrier = &barrier;
            let borrowed_final_set = &final_set;
            s.spawn(move || {
                if i != num_threads-1 {
                	k_rho_leonhardt_worker_base_new_ud(graph,(i*work_size) as u32..((i+1)*work_size) as u32,k,rho,borrowed_set, borrowed_counter, &borrowed_sync);
                }
                else {
                	k_rho_leonhardt_worker_base_new_ud(graph,(i*work_size) as u32..graph.num_nodes(),k,rho,borrowed_set, borrowed_counter, &borrowed_sync);
                }
                borrowed_barrier.wait();
                // First thread calculates average and standard deviation
                let sync_lock = borrowed_sync.read().unwrap();

                let (avg,mut standard_deviation) = sync_lock.calculate_mean_stddev(borrowed_set.len());
                standard_deviation*=3.0;

                let shards = borrowed_set.shards();
                let shards_len = shards.len();

                let total_work = shards_len/num_threads;
                let start = total_work*i;
                let end =  if i != num_threads-1 {total_work*(i+1)} else {shards_len};

                for sh in &shards[start..end] {
                    let lck = sh.read();

                    for (key, value_s) in lck.iter() {
                        let value = value_s.get();
                        if value.fractional_distinct_vertices > avg + standard_deviation && value.distinct_vertices > k && value.num_times > value.distinct_vertices*(k-1)/2 {
                            borrowed_final_set.insert(*key, value.distance);
                        }
                    }
                }
            });
        }
    });

    let mut cloned_graph = graph.clone();

    for entry in final_set.iter() {
        let target = *entry.key()  % graph.num_nodes() as u64;
        let source = *entry.key() / graph.num_nodes() as u64;

        if !cloned_graph.has_edge(source as u32,target as u32) {
            cloned_graph.add_edge(source as u32, Edge{ weight: *entry.value(), target: target as u32});
        }
        else {
            panic!("The edge should not exist!");
        }
    }

    let sol = k_rho_blelloch_par(&cloned_graph,k,rho,num_threads);
    final_set.extend(sol.iter().map(|entry| (*entry.key(), *entry.value())));
    final_set
}


pub fn k_rho_leonhardt_worker<'a, I>(graph: &WeightedDirectedGraph, node_iterator: I, k: u32, rho: u32, shortcut_set: &DashMap<u64,EdgeFrequencyDescription>, fingerprint: &DashMap<u64,MinHashSignature>, minhash: &MinHash, counter: &std::sync::atomic::AtomicU32, sync: &std::sync::RwLock<SyncStats>)
where
    I: Iterator<Item=u32>
{
    let mut dijk = DijkstraFewestHops::new(graph);
    let mut local_sync = SyncStats::new();
    let mut local_counter:u32 = 0;

    for source in node_iterator {
        let mut covered_pairs = std::collections::HashSet::<u64>::new();
        dijk.run(source);
        let sorted = dijk.get_nodes_sorted_by_distance();

        let mut unique_children = std::collections::HashMap::<u32,(u32,u32)>::new();
        let mut children_tree: Vec<Vec<u32>> = vec![Vec::new(); graph.num_nodes() as usize];
        for VertexDistanceDescription{vertex,..} in &sorted[1..rho as usize+1] {
            let parent = dijk.parent(*vertex as u32);
            children_tree[parent as usize].push(*vertex as u32);
        }

        let mut current_frontier: Vec<(u32,u32)> = vec![(source,source)];
        let mut next_frontier:Vec<(u32,u32)> = Vec::new();

        while current_frontier.len() > 0 {
            for (parent,node) in current_frontier.iter() {
                if children_tree[*node as usize].len() == 0 {
                    if !unique_children.contains_key(parent) {
                        unique_children.insert(*parent, (*node,*node));
                    }
                    else {
                        let res = unique_children.get_mut(parent).unwrap();
                        if res.0 == res.1 {
                            res.1 = *node;
                        }
                    }
                    continue;
                }
        	for child in &children_tree[*node as usize] {
        		next_frontier.push((*node,*child));
            	}
            }
            current_frontier = next_frontier;
            next_frontier = Vec::new();
        }


        let total: usize = sorted[1..rho as usize+1].iter().filter(|x| {
            let path = dijk.get_path(x.vertex);

            if children_tree[x.vertex as usize].len() == 0 {

                let ent = unique_children.get(&path[path.len()-2]).unwrap();
                if ent.0 != x.vertex && ent.1 != x.vertex {
                    return false;
                }
            }
            return true;
            
                }).map(|x| if dijk.get_path(x.vertex).len() <= k as usize+1 { 0 } else { 1 }).sum();
        for VertexDistanceDescription{vertex,..} in &sorted[1..rho as usize+1] {
            let path = dijk.get_path(*vertex as u32);

            if path.len()-1 <= k as usize {
                continue;
            }

            // Skip not unique children
            if children_tree[*vertex as usize].len() == 0 {
                let ent = unique_children.get(&path[path.len()-2]).unwrap();
                if (ent.0 != *vertex) && (ent.1 != *vertex) {
                    continue;
                }
            }

            local_counter+=1;
            let outer_key = graph.edge_uuid(source, *vertex);

            for (index,node) in (&path[0..k as usize]).iter().enumerate() {
                let start_index = path.len()-(k as usize-index);
                if index+1 == start_index {
                    panic!("Should never happen!");
                }

                for pair in &path[start_index..] {
                    let key = graph.edge_uuid(*node,*pair);

		    fingerprint.entry(key).and_modify(|e| minhash.add(outer_key,e)).or_insert(minhash.generate_signature(outer_key));

                    let distinct = if !covered_pairs.contains(&key) {
                        covered_pairs.insert(key);
			1
                    }
                    else {
                        0
                    };

                    let item = shortcut_set.entry(key);
                    match item {
                        Entry::Occupied(mut entry_lck) => {
                            let entry = entry_lck.get_mut();
                            entry.num_times+=1;
                            local_sync.update(1.0/total as f64, entry.fractional_distinct_vertices);
                            entry.fractional_distinct_vertices += 1.0/total as f64;
                            entry.distinct_vertices+=distinct;
                        }
                        Entry::Vacant(entry) => {
                            entry.insert(EdgeFrequencyDescription{num_times: 1, fractional_distinct_vertices: 1.0/total as f64, distinct_vertices: 1, distance: dijk.distance(*pair)-dijk.distance(*node)});
                            local_sync.update(1.0/total as f64, 0.0);
                        }
                    }
                }
            }
        }
    }

    counter.fetch_add(local_counter,std::sync::atomic::Ordering::Relaxed);
    let mut lck = sync.write().unwrap();
    lck.fraction_distinct_vertices_sum += local_sync.fraction_distinct_vertices_sum;
    lck.fraction_distinct_vertices_squared_sum += local_sync.fraction_distinct_vertices_squared_sum;
}

pub fn k_rho_leonhardt_perturbation_par(graph: &WeightedDirectedGraph, k: u32, rho: u32, num_threads: usize) -> DashMap<u64,u32> {
    let work_size = graph.num_nodes() as usize/num_threads as usize;

    let mut cloned_graph = graph.clone();

    let optimal_cuts = DashMap::<u64,(u32,u32)>::new();

    std::thread::scope(|s| {
        for i in 0..num_threads {
            let borrowed_cuts = &optimal_cuts;
            s.spawn(move || {
                if i != num_threads-1 {
                	k_rho_leonhardt_worker_dp(graph,(i*work_size) as u32..((i+1)*work_size) as u32,k,rho,borrowed_cuts);
                }
                else {
                	k_rho_leonhardt_worker_dp(graph,(i*work_size) as u32..graph.num_nodes(),k,rho,borrowed_cuts);
                }
            });
        }
    });

    let mut final_set = DashMap::<u64,u32>::new();

    let res : Vec<(u64,u32,u32)> = optimal_cuts.iter().map(|x| (*x.key(),x.value().0, x.value().1)).collect();

    for last in res.iter() {
        if last.1 < 2 {
            continue;
        }
        let (source,target) = cloned_graph.uuid_to_edge(last.0);

        if !cloned_graph.has_edge(source,target) {
            cloned_graph.add_edge(source, Edge{ weight: last.2, target: target});
        }
        else {
            panic!("The edge should not exist!");
        }
        final_set.insert(last.0,last.2);
    }

    let sol = k_rho_blelloch_par(&cloned_graph,k,rho,num_threads);
    final_set.extend(sol.iter().map(|x| (*x.key(),*x.value())));
    final_set
}


pub fn k_rho_leonhardt_perturbation_minhash_par(graph: &WeightedDirectedGraph, k: u32, rho: u32, num_threads: usize) -> DashMap<u64,u32> {
    let work_size = graph.num_nodes() as usize/num_threads as usize;

    let mut cloned_graph = graph.clone();

    let optimal_cuts = DashMap::<u64,(u32,u32)>::new();

    const MINHASH_NUM_BANDS:usize = 10;

    let mut rng = rand::thread_rng();
    let min_hash = MinHash::new(MINHASH_NUM_BANDS, &mut rng);
    let finger_print = std::sync::Arc::new(DashMap::<u64,MinHashSignature>::new());

    std::thread::scope(|s| {
        for i in 0..num_threads {
            let borrowed_cuts = &optimal_cuts;
            let cloned_fingerprint = finger_print.clone();
            let borrowed_minhash = &min_hash;
            s.spawn(move || {
                if i != num_threads-1 {
                	k_rho_leonhardt_worker_dp_minhash(graph,(i*work_size) as u32..((i+1)*work_size) as u32,k,rho, borrowed_cuts, &cloned_fingerprint, borrowed_minhash);
                }
                else {
                	k_rho_leonhardt_worker_dp_minhash(graph,(i*work_size) as u32..graph.num_nodes(),k,rho, borrowed_cuts, &cloned_fingerprint, borrowed_minhash);
                }
            });
        }
    });

    let mut final_set = DashMap::<u64,u32>::new();

    let res : Vec<(u64,u32,u32)> = optimal_cuts.iter().map(|x| (*x.key(),x.value().0, x.value().1)).collect();
    let mut jumped = std::collections::HashSet::<(u32,u32)>::new();

    for last in res.iter() {
        if last.1 < 2 {
            continue;
        }

        if let Some(sig) = finger_print.get(&last.0) {
            let hit:usize = sig.bands.iter().enumerate().map(|(idx,x)| if jumped.contains(&(idx as u32,x.signature)) {1} else {0}).sum();
            if (last.1*(MINHASH_NUM_BANDS as u32 - hit as u32)/MINHASH_NUM_BANDS as u32) < 2 {
                continue;
            }

            sig.bands.iter().enumerate().for_each(|(idx,x)| {jumped.insert((idx as u32,x.signature));});
        }
        let (source, target) = graph.uuid_to_edge(last.0);

        if !cloned_graph.has_edge(source as u32,target as u32) {
            cloned_graph.add_edge(source as u32, Edge{ weight: last.2, target: target as u32});
        }
        else {
            panic!("The edge should not exist!");
        }
        final_set.insert(last.0,last.2);
    }

    let sol = k_rho_blelloch_par(&cloned_graph,k,rho,num_threads);
    final_set.extend(sol.iter().map(|x| (*x.key(), *x.value())));
    final_set
}

pub fn k_rho_leonhardt_perturbation_minhash_multi_stage_par(graph: &WeightedDirectedGraph, k: u32, rho: u32, num_threads: usize) -> DashMap<u64,u32> {
    let work_size = graph.num_nodes() as usize/num_threads as usize;

    let mut cloned_graph = graph.clone();

    let optimal_cuts = DashMap::<u64,(u32,u32)>::new();

    const MINHASH_NUM_BANDS:usize = 10;

    let mut rng = rand::thread_rng();
    let min_hash = MinHash::new(MINHASH_NUM_BANDS, &mut rng);
    let finger_print = std::sync::Arc::new(DashMap::<u64,MinHashSignature>::new());

    std::thread::scope(|s| {
        for i in 0..num_threads {
            let borrowed_cuts = &optimal_cuts;
            let cloned_fingerprint = finger_print.clone();
            let borrowed_minhash = &min_hash;
            s.spawn(move || {
                if i != num_threads-1 {
                	k_rho_leonhardt_worker_dp_minhash(graph,(i*work_size) as u32..((i+1)*work_size) as u32,k,rho, borrowed_cuts, &cloned_fingerprint, borrowed_minhash);
                }
                else {
                	k_rho_leonhardt_worker_dp_minhash(graph,(i*work_size) as u32..graph.num_nodes(),k,rho, borrowed_cuts, &cloned_fingerprint, borrowed_minhash);
                }
            });
        }
    });

    let mut final_set = DashMap::<u64,u32>::new();

    let res : Vec<(u64,u32,u32)> = optimal_cuts.iter().map(|x| (*x.key(),x.value().0, x.value().1)).collect();
    let mut jumped = std::collections::HashSet::<(u32,u32)>::new();

    for last in res.iter() {
        if last.1 < 2 {
            continue;
        }

        if let Some(sig) = finger_print.get(&last.0) {
            let hit:usize = sig.bands.iter().enumerate().map(|(idx,x)| if jumped.contains(&(idx as u32,x.signature)) {1} else {0}).sum();
            if (last.1*(MINHASH_NUM_BANDS as u32 - hit as u32)/MINHASH_NUM_BANDS as u32) < 2 {
                continue;
            }

            sig.bands.iter().enumerate().for_each(|(idx,x)| {jumped.insert((idx as u32,x.signature));});
        }
        let (source, target) = graph.uuid_to_edge(last.0);

        if !cloned_graph.has_edge(source as u32,target as u32) {
            cloned_graph.add_edge(source as u32, Edge{ weight: last.2, target: target as u32});
        }
        else {
            panic!("The edge should not exist!");
        }
        final_set.insert(last.0,last.2);
    }

    let sol = k_rho_blelloch_star_par(&cloned_graph,k,rho,10,num_threads);
    final_set.extend(sol.iter().map(|x| (*x.key(), *x.value())));
    final_set
}

pub fn k_rho_leonhardt_fingerprint_par(graph: &WeightedDirectedGraph, k: u32, rho: u32, num_threads: usize) -> DashMap<u64,u32> {
    let (cloned_graph,mut final_set) = __internal_k_rho_leonhardt_pc_minhash(graph,k,rho,num_threads);

    let sol = k_rho_blelloch_par(&cloned_graph,k,rho,num_threads);
    final_set.extend(sol.iter().map(|entry| (*entry.key(), *entry.value())));
    final_set
}


fn __internal_k_rho_leonhardt_pc_minhash(graph: &WeightedDirectedGraph, k: u32, rho: u32, num_threads: usize) -> (WeightedDirectedGraph,DashMap<u64,u32>) {

    let shortcut_incidence_map = DashMap::<u64,EdgeFrequencyDescription>::new();
    let finger_print = DashMap::<u64,MinHashSignature>::new();
    let work_size = graph.num_nodes() as usize/num_threads as usize;


    const MINHASH_NUM_BANDS:usize = 5;

    let mut rng = rand::thread_rng();
    let min_hash = MinHash::new(MINHASH_NUM_BANDS,&mut rng);

    let atomic_counter = std::sync::atomic::AtomicU32::new(0);
    let sync = std::sync::RwLock::new(SyncStats::new());
    let barrier = Barrier::new(num_threads);
    let jumped = DashSet::<(u32,u32)>::new();


    let final_set = DashMap::<u64,u32>::new();
    std::thread::scope(|s| {
        for i in 0..num_threads {
            let borrowed_incidence_map = &shortcut_incidence_map;
            let borrowed_fingerprint_map = &finger_print;
            let borrowed_counter = &atomic_counter;
            let borrowed_minhash = &min_hash;
            let borrowed_sync = &sync;
            let borrowed_barrier = &barrier;
            let borrowed_final_set = &final_set;
            let borrowed_jumped = &jumped;
            s.spawn(move || {
                if i != num_threads-1 {
                	k_rho_leonhardt_worker(graph,(i*work_size) as u32..((i+1)*work_size) as u32,k,rho,borrowed_incidence_map, borrowed_fingerprint_map, borrowed_minhash, borrowed_counter, &borrowed_sync);
                }
                else {
                	k_rho_leonhardt_worker(graph,(i*work_size) as u32..graph.num_nodes(),k,rho, borrowed_incidence_map, borrowed_fingerprint_map, borrowed_minhash, borrowed_counter, &borrowed_sync);
                }

                borrowed_barrier.wait();
                // First thread calculates average and standard deviation
                let sync_lock = borrowed_sync.read().unwrap();

                let (avg,mut standard_deviation) = sync_lock.calculate_mean_stddev(borrowed_incidence_map.len());
                standard_deviation*=3.0;

                let shards = borrowed_incidence_map.shards();
                let shards_len = shards.len();

                let total_work = shards_len/num_threads;
                let start = total_work*i;
                let end =  if i != num_threads-1 {total_work*(i+1)} else {shards_len};

                for sh in &shards[start..end] {
                    let lck = sh.read();

                    for (key, value_s) in lck.iter() {
                        let value = value_s.get();
                        if value.fractional_distinct_vertices > avg + standard_deviation && value.distinct_vertices > k && value.num_times > value.distinct_vertices*(k-1)/2 {
                            if let Some(sig) = borrowed_fingerprint_map.get(key) {
                                let hit:usize = sig.bands.iter().enumerate().map(|(idx,x)| if borrowed_jumped.contains(&(idx as u32,x.signature)) {1} else {0}).sum();
                                if hit > 3 {
                                    continue;
                                }

                                sig.bands.iter().enumerate().for_each(|(idx,x)| {borrowed_jumped.insert((idx as u32,x.signature));});
                            }
                            borrowed_final_set.insert(*key, value.distance);
                        }
                    }
                }
            });
        }
    });

    let mut cloned_graph = graph.clone();

    for entry in final_set.iter() {
        let target = *entry.key()  % graph.num_nodes() as u64;
        let source = *entry.key() / graph.num_nodes() as u64;

        if !cloned_graph.has_edge(source as u32,target as u32) {
            cloned_graph.add_edge(source as u32, Edge{ weight: *entry.value(), target: target as u32});
        }
        else {
            panic!("The edge should not exist!");
        }
    }

    (cloned_graph,final_set)
}


pub fn k_rho_leonhardt_perturbation_minhash_pair_cutting_par(graph: &WeightedDirectedGraph, k: u32, rho: u32, num_threads: usize) -> DashMap<u64,u32> {
    let (cloned_graph,mut final_set) = __internal_k_rho_leonhardt_pc_minhash(graph,k,rho,num_threads);

    let sol = k_rho_leonhardt_perturbation_minhash_par(&cloned_graph,k,rho,num_threads);
    final_set.extend(sol.iter().map(|entry| (*entry.key(), *entry.value())));
    final_set
}

pub fn k_rho_leonhardt_perturbation_minhash_pair_cutting_multi_stage_par(graph: &WeightedDirectedGraph, k: u32, rho: u32, num_threads: usize) -> DashMap<u64,u32> {
    let (cloned_graph,mut final_set) = __internal_k_rho_leonhardt_pc_minhash(graph,k,rho,num_threads);

    let sol = k_rho_leonhardt_perturbation_minhash_multi_stage_par(&cloned_graph,k,rho,num_threads);
    final_set.extend(sol.iter().map(|entry| (*entry.key(), *entry.value())));
    final_set
}
