use std::collections::HashSet;
use std::collections::VecDeque;
use std::collections::HashMap;

pub struct DAG {
    graph: Option<HashMap<usize, Vec<usize>>>,
    weights: HashMap<(usize, usize), usize>,
    count: usize,
}

impl DAG {
    pub fn new(graph_info: Vec<(usize, usize, usize)>, node_count: usize) -> Self {
        // DirectedGraph { graph: None }
        let mut adjacency_list: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut weights: HashMap<(usize, usize), usize> = HashMap::new();

        let graph = graph_info.get(0..);
        for value in graph.unwrap() {

            let source_vertex = &mut adjacency_list.entry(value.0).or_insert(vec![]);
            source_vertex.push(value.1);
            *weights.entry((value.0, value.1)).or_insert(0) = value.2;
        }
        let the_graph = DAG {
            graph: Some(adjacency_list),
            weights: weights,
            count: node_count,
        };
        the_graph
    }

    pub fn get_topological_order(&self, node: &usize) -> Vec<usize> {
        let mut stack: Vec<usize> = vec![];
        let mut visited: Vec<bool> = vec![false; self.count];
        if !visited[*node] {
            self.get_order(&node, &mut stack, &mut visited);
        }
        return stack;
    }

    pub fn get_order(&self, node: &usize, stack: &mut Vec<usize>, visited: &mut Vec<bool>) {
        visited[*node] = true;
        let receiving_nodes = self.graph.as_ref().unwrap().get(node);
        // if visited.get(node) == None {
        if receiving_nodes != None {
            for value in receiving_nodes.unwrap() {
                if !visited[*value] {
                    self.get_order(value, stack, visited);
                }
            }
        }
        if !stack.contains(node) {
            stack.push(*node);
        }
        // }
    }

    pub fn longest_paths(&self) -> Vec<(usize, Vec<usize>, usize)> {
        let mut paths_per_node = Vec::<(usize, Vec<usize>, usize)>::new();
        let source_nodes = self.graph.as_ref().unwrap().keys();
        for node in source_nodes {
            let V = self.count;
            let mut stack = self.get_topological_order(node);
            let mut scores: Vec<i32> = vec![i32::MIN; V];
            let mut parents: Vec<i32> = vec![-1; V];
            scores[*node] = 0 as i32;
            while !stack.is_empty() {
                let mut u = stack.pop().unwrap();
                if scores[u] != i32::MIN {
                    let receiving_nodes = self.graph.as_ref().unwrap().get(&u);
                    if receiving_nodes != None {
                        for v in receiving_nodes.unwrap() {
                            if scores[*v] < scores[u] + self.weights[&(u, *v)] as i32 {
                                scores[*v] = scores[u] + self.weights[&(u, *v)] as i32;
                                parents[*v] = u as i32; 
                            }
                        }  
                    }
                }
            }
            let (mut k, s) = scores.iter().enumerate().max_by(|a, b| a.1.cmp(&b.1)).unwrap();
            let mut final_path = Vec::<usize>::new();
            while parents[k] != -1 {
                final_path.push(k);
                k = parents[k] as usize;
            }
            final_path.reverse();
            paths_per_node.push((*node, final_path, *s as usize));
        }
        paths_per_node
    }
}