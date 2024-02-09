#![feature(isqrt)]
use std::env;
use std::process::Command;
use dashmap::DashMap;
use std::io::Read;
use serde::Serialize;
use serde_json;

mod min_hash;
mod weighted_graph;
mod dijkstra;
mod util;
mod blelloch;
mod new_heuristics;



use util::{EdgeFrequencyDescription};
use min_hash::{MinHash, MinHashSignature};
use weighted_graph::{Edge,WeightedDirectedGraph};
use dijkstra::{VertexDistanceDescription};
use blelloch::{k_rho_blelloch_greedy_par, k_rho_blelloch_star_par, k_rho_blelloch_par, k_rho_blelloch_bot_par};
use new_heuristics::{k_rho_leonhardt_fingerprint_par, k_rho_leonhardt_perturbation_par, k_rho_leonhardt_perturbation_minhash_par,k_rho_leonhardt_perturbation_minhash_pair_cutting_par,k_rho_leonhardt_base_ud_par,k_rho_leonhardt_perturbation_minhash_pair_cutting_multi_stage_par};


pub fn assert_valid_shortcut_set_k_rho_graph<I>(graph: &WeightedDirectedGraph, shortcut_set: I, k: u32, rho: u32, threads: usize)
where
    I: Iterator<Item=(u64,u32)>
{
    let mut cloned = graph.clone();
    for (edge_uuid, weight) in shortcut_set {
        let (source,target) = cloned.uuid_to_edge(edge_uuid);
        if !cloned.has_edge(source,target) {
            cloned.add_edge(source, Edge{ weight: weight, target: target});
        }
        else {
            panic!("The edge should not exist!");
        }
    }

    if !cloned.is_k_rho_graph(k,rho,threads) {
        panic!("The graph is not a ({},{})-graph!",k,rho);
    }
}

pub fn bench_and_check_algorithm<F: Fn(&WeightedDirectedGraph,u32,u32,usize) -> DashMap<u64,u32>>(name: &str, graph: &WeightedDirectedGraph, k: u32, rho: u32, func: F, threads: usize, surpress_output: bool) -> (usize, usize) {
    if !surpress_output {
        println!("Running now {}",name);
    }
    let now = std::time::Instant::now();
    let result = func(graph,k,rho,threads);
    let elapsed_ms = now.elapsed().as_millis();
    if !surpress_output {
        println!("Solution size {} in {}ms\n",result.len(),elapsed_ms);
    }

    let result_len = result.len();
    assert_valid_shortcut_set_k_rho_graph(graph, result.into_iter(), k, rho, threads);
    (result_len, elapsed_ms as usize)
}

#[derive(Serialize)]
pub struct ExperimentalResults {
    size_bl_dp: usize,
    speed_bl_dp_ms: usize,
    size_bl_dp_s: usize,
    speed_bl_dp_s_ms: usize,
    size_bl_dp_g: usize,
    speed_bl_dp_g_ms: usize,
    size_bl_dp_pc: usize,
    speed_bl_dp_pc_ms: usize,
    size_bl_dp_pc_min: usize,
    speed_bl_dp_pc_min_ms: usize,
    size_bl_dp_pert: usize,
    speed_bl_dp_pert_ms: usize,
    size_bl_dp_pert_min: usize,
    speed_bl_dp_pert_min_ms: usize,
    size_bl_dp_pert_min_pc: usize,
    speed_bl_dp_pert_min_pc_ms: usize,
    size_bl_dp_pert_min_pc_ms: usize,
    speed_bl_dp_pert_min_pc_ms_ms: usize,
    build_hash: String
}


fn main() {
    // Arguments program generator n avg_degree k rho output-file-name
    let args: Vec<String> = env::args().collect();

    // Read graph from stdin
    if args.len() == 3 {
        let mut reader = std::io::stdin();
        let mut data = String::new();
        match reader.read_to_string(&mut data) {
            Ok(num_bytes) => if num_bytes == 0 {panic!("Need to read more than 0 bytes");},
            Err(k) => panic!("Error {}",k),
        }


        let k:u32 = args[1].parse::<u32>().expect("Could not read k from the input");
        let rho:u32 = args[2].parse::<u32>().expect("Could not read rho from the input");
        let g = WeightedDirectedGraph::from_generator(&data).expect("Could not read graph from stdin");
        let threads = std::thread::available_parallelism().unwrap_or(unsafe{std::num::NonZeroUsize::new_unchecked(1)}).get();


        let (size_bl_dp, speed_bl_dp_ms) = bench_and_check_algorithm("Blelloch DP", &g, k, rho, |graph,kp,rp,th| k_rho_blelloch_par(graph,kp,rp,th),threads,true);
        let (size_bl_dp_s, speed_bl_dp_s_ms) = bench_and_check_algorithm("Blelloch DP*", &g, k, rho, |graph,kp,rp,th| k_rho_blelloch_star_par(graph,kp,rp,10,th),threads,true);
        let (size_bl_dp_g, speed_bl_dp_g_ms) = bench_and_check_algorithm("Blelloch Greedy", &g, k, rho, |graph,kp,rp,th| k_rho_blelloch_greedy_par(graph,kp,rp,th),threads,true);
        let (size_bl_dp_pc, speed_bl_dp_pc_ms) = bench_and_check_algorithm("Leonhardt DP-PC", &g, k, rho, |graph,kp,rp,th| k_rho_leonhardt_base_ud_par(graph,kp,rp,th),threads,true);
        let (size_bl_dp_pc_min, speed_bl_dp_pc_min_ms) = bench_and_check_algorithm("Leonhardt DP-PC + MinHash", &g, k, rho, |graph,kp,rp,th| k_rho_leonhardt_fingerprint_par(graph,kp,rp,th),threads,true);
        let (size_bl_dp_pert, speed_bl_dp_pert_ms) = bench_and_check_algorithm("Leonhardt Pert", &g, k, rho, |graph,kp,rp,th| k_rho_leonhardt_perturbation_par(graph,kp,rp,th),threads,true);
        let (size_bl_dp_pert_min, speed_bl_dp_pert_min_ms) = bench_and_check_algorithm("Leonhardt Pert + MinHash", &g, k, rho, |graph,kp,rp,th| k_rho_leonhardt_perturbation_minhash_par(graph,kp,rp,th),threads,true);
        let (size_bl_dp_pert_min_pc, speed_bl_dp_pert_min_pc_ms) = bench_and_check_algorithm("Leonhardt DP-PC + Pert + MinHash", &g, k, rho, |graph,kp,rp,th| k_rho_leonhardt_perturbation_minhash_pair_cutting_par(graph,kp,rp,th),threads,true);
        let (size_bl_dp_pert_min_pc_ms, speed_bl_dp_pert_min_pc_ms_ms) = bench_and_check_algorithm("Leonhardt DP-PC + Pert + MinHash + MS", &g, k, rho, |graph,kp,rp,th| k_rho_leonhardt_perturbation_minhash_pair_cutting_multi_stage_par(graph,kp,rp,th),threads,true);
        let result: ExperimentalResults = ExperimentalResults {
            size_bl_dp,
            speed_bl_dp_ms,
            size_bl_dp_s,
            speed_bl_dp_s_ms,
            size_bl_dp_g,
            speed_bl_dp_g_ms,
            size_bl_dp_pc,
            speed_bl_dp_pc_ms,
            size_bl_dp_pc_min,
            speed_bl_dp_pc_min_ms,
            size_bl_dp_pert,
            speed_bl_dp_pert_ms,
            size_bl_dp_pert_min,
            speed_bl_dp_pert_min_ms,
            size_bl_dp_pert_min_pc,
            speed_bl_dp_pert_min_pc_ms,
            size_bl_dp_pert_min_pc_ms,
            speed_bl_dp_pert_min_pc_ms_ms,
            build_hash: env!("GIT_HASH").to_string()
        };

        println!("{}",serde_json::to_string(&result).unwrap());
    }
    // Generate graph
    else if args.len() < 6 {
        panic!("Usage: program generator-name n avg_degree k rho\n generator-name --- gilbert,hyperbolic,powerlaw!");
    }
    else {
        let result = Command::new("../../analyse/bin/python")
            .args(["src/generator.py", &args[1],&args[2], &args[3]])
            .output()
            .expect("Error from generator");
        let k:u32 = args[4].parse::<u32>().expect("Could not read k from the input");
        let mut rho:u32 = args[5].parse::<u32>().expect("Could not read rho from the input");
        let result_str = std::str::from_utf8(&result.stdout).expect("Could not transform generator output to utf-8");
        if result.status.success() {
            // Generates undirected graphs which are converted to weighted directed graphs here
            let g = WeightedDirectedGraph::from_generator(result_str).expect("Could not read generated graph from stdout");
            println!("Largest connected component of generated graph has size {}",g.num_nodes());
            if rho >= g.num_nodes() {
                rho = g.num_nodes()-1;
                println!("Adjusting rho to {}\n",rho);
            }

            let threads = std::thread::available_parallelism().unwrap_or(unsafe{std::num::NonZeroUsize::new_unchecked(1)}).get();

            bench_and_check_algorithm("Blelloch DP", &g, k, rho, |graph,kp,rp,th| k_rho_blelloch_par(graph,kp,rp,th),threads,false);
            bench_and_check_algorithm("Blelloch DP Bot", &g, k, rho, |graph,kp,rp,th| k_rho_blelloch_bot_par(graph,kp,rp,th),threads,false);
            bench_and_check_algorithm("Blelloch DP*", &g, k, rho, |graph,kp,rp,th| k_rho_blelloch_star_par(graph,kp,rp,10,th),threads,false);
            bench_and_check_algorithm("Blelloch Greedy", &g, k, rho, |graph,kp,rp,th| k_rho_blelloch_greedy_par(graph,kp,rp,th),threads,false);
            bench_and_check_algorithm("Leonhardt DP-PC", &g, k, rho, |graph,kp,rp,th| k_rho_leonhardt_base_ud_par(graph,kp,rp,th),threads,false);
            bench_and_check_algorithm("Leonhardt DP-PC + MinHash", &g, k, rho, |graph,kp,rp,th| k_rho_leonhardt_fingerprint_par(graph,kp,rp,th),threads,false);
            bench_and_check_algorithm("Leonhardt Pert", &g, k, rho, |graph,kp,rp,th| k_rho_leonhardt_perturbation_par(graph,kp,rp,th),threads,false);
            bench_and_check_algorithm("Leonhardt Pert + MinHash", &g, k, rho, |graph,kp,rp,th| k_rho_leonhardt_perturbation_minhash_par(graph,kp,rp,th),threads,false);
            bench_and_check_algorithm("Leonhardt DP-PC + Pert + MinHash", &g, k, rho, |graph,kp,rp,th| k_rho_leonhardt_perturbation_minhash_pair_cutting_par(graph,kp,rp,th),threads,false);
            bench_and_check_algorithm("Leonhardt DP-PC + Pert + MinHash + MS", &g, k, rho, |graph,kp,rp,th| k_rho_leonhardt_perturbation_minhash_pair_cutting_multi_stage_par(graph,kp,rp,th),threads,false);
        }
        else {
            let result_err_str = std::str::from_utf8(&result.stderr).expect("Could not transform generator output to utf-8");
            println!("{} status: {:?}",result_err_str, result.status);
            std::process::exit(-1);
        }
    }
}
