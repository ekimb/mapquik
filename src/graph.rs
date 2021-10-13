use petgraph::prelude::*;
use petgraph::Graph;
use std::sync::atomic::{AtomicUsize, Ordering};
use dashmap::{DashMap, DashSet};
pub type Kmer = KmerVec;
pub type DbgIndex = u32;// heavily optimized assuming we won't get more than 2B kminmers of abundance <= 65535
pub type DbgAbundance = u16;
use crate::kmer_vec::KmerVec;
#[derive(Debug, Clone)]
pub struct DbgEntry {pub index: DbgIndex, pub mers: Vec<Kmer>, pub abundance: DbgAbundance, pub origin: Vec<String>} 

pub fn add_kminmers(seq_id: &str, ref_mers: &Vec<Kmer>, dbg_nodes: &DashMap<Vec<u64>, DbgEntry>, NODE_INDEX: &AtomicUsize, dbg_edges: &DashMap<DbgIndex, Vec<Kmer>>) {
    for i in 0..ref_mers.len() {
        let mut cur_node_index: DbgIndex = 0 as DbgIndex;
        let mut contains_key;
        contains_key = dbg_nodes.contains_key(&ref_mers[i].normalize().0.hashes);
        if contains_key {
            let mut entry_mut = dbg_nodes.get_mut(&ref_mers[i].normalize().0.hashes).unwrap();
            cur_node_index = entry_mut.index;
            entry_mut.abundance += 1;
            if !entry_mut.origin.contains(&seq_id.to_string()) {
                entry_mut.origin.push(seq_id.to_string());
            }
            entry_mut.mers.push(ref_mers[i].clone());
        }
        else {
            cur_node_index = NODE_INDEX.fetch_add(1, Ordering::Relaxed) as DbgIndex;
            dbg_nodes.insert(ref_mers[i].normalize().0.hashes.clone(), DbgEntry{index: cur_node_index, abundance: 1, mers: vec![ref_mers[i].clone()], origin: vec![seq_id.to_string()]}); 
            contains_key = true;
        }
        if i < ref_mers.len() - 1 {
            let mut contains_edge;
            contains_edge = dbg_edges.contains_key(&cur_node_index);
            if contains_edge {
                let mut entry_mut = dbg_edges.get_mut(&cur_node_index).unwrap();
                entry_mut.push(ref_mers[i+1].clone());
            }
            else {
                dbg_edges.insert(cur_node_index, vec![ref_mers[i+1].clone()]); 
            }
        }
    }

}

