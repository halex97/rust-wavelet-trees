use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;
use petgraph::EdgeType;
use petgraph::graph::{Graph, Neighbors, IndexType};
use serde::{Serialize, Deserialize};
use crate::WaveletTree;
use crate::pointerless::PointerlessWaveletTree;

#[derive(Serialize, Deserialize)]
pub struct GraphWaveletTree {
    /// The wavelet tree represents the concatenated adjacency lists of the graph. The node indices
    /// are stored as usize integers.
    tree: PointerlessWaveletTree<usize>,
    /// The bitmap marks the beginnings of the adjacency lists
    bitmap: RankSelect
}

impl GraphWaveletTree {

    /// Creates a new GraphWaveletTree from a graph. The resulting GraphWaveletTree's indices are
    /// guaranteed to be in the same order as the graph's indices, but they are NOT guaranteed to
    /// be THE SAME indices.
    pub fn from_graph<N,E,Ty:EdgeType,Ix:IndexType> (graph: Graph<N,E,Ty,Ix>) -> Self {
        let nodes = graph.node_indices();

        let mut lists = Vec::with_capacity(nodes.len());

        for node in nodes {
            lists.push(graph.neighbors(node));
        }

        Self::from_adjacency_lists(lists)
    }

    pub fn from_adjacency_lists<E,Ix:IndexType> (lists: Vec<Neighbors<E, Ix>>) -> Self {
        // Build the sequence (which is the concatenation of the adjacency lists) and the
        // corresponding bit vector marking the beginnings of the individual lists in the sequence.
        let mut sequence : Vec<usize> = Vec::new();
        let mut lists_bv : BitVec<u8> = BitVec::new();

        // Iterate over all adjacency lists
        for neighbors in lists {
            // A '1' in the bit vector separates two adjacency lists.
            lists_bv.push(true);

            // Create a sub sequence to collect the current adjacency list.
            // This needs to be done because neighbors are listed in REVERSE order.
            let mut subseq = Vec::new();

            // Iterate over all neighbors (i.e. adjacent nodes that appear in the list)
            for n in neighbors {
                // Add the neighbor's index to the end of the sub sequence 
                subseq.push(n.index());

                // '0's in the bit vector represent node indices in the adjacency lists
                lists_bv.push(false);
            }

            // Append the reversed sub sequence.
            subseq.reverse();
            sequence.append(&mut subseq);
        }

        // Return a pointerless wavelet tree representing the sequence together with the bitmap.
        GraphWaveletTree {
            tree: PointerlessWaveletTree::from_slice(sequence.as_slice()),
            bitmap: RankSelect::new(lists_bv, 1)
        }
    }

    /// Access the i-th neighbour (i.e. successor) of the node given by index v
    pub fn access_neighbor(&self, v: usize, i: usize) -> Option<usize> {
        // First we need to find the start of the adjacency list l.
        // Then we need to take into account that the bitmap is longer than the actual sequence
        // because of the additional '1's.
        let l = self.bitmap.select_1((v+1) as u64)
            .map(|x| x - v as u64);
        
        // Finally we can add the index i to find the correct neighbor's index in the sequence of
        // concatenated adjacency lists.
        let n = l.map(|l| l + (i-1) as u64);

        // Now we look for the end of adjacency list. This is either the start of the next 
        // adjacency list or the end of the list of all concatenated adjacency lists.
        let m = self.bitmap.select_1((v+2) as u64)
            .map(|l| l - (v+1) as u64)
            .unwrap_or(self.sequence_length());

        //println!("Adjacency list: [{:?}, {:?}]. n: {:?}", l, m, n);

        // If n is no actual index inside the correct adjacency list, no neighbor can be returned.
        // Otherwise, we can look for the neighbor's node index in the wavelet tree.
        if n.is_none() || n.unwrap() >= m {
            None
        } else {
            n.and_then(|x| self.tree.access(x)).map(|x| x.clone())
        }        
    }

    /// Access the i-th reverse neighbor (i.e. predecessor) of node v
    pub fn access_reverse_neighbor(&self, v: usize, i: usize) -> Option<usize> {
        // First we need to find the i-th occurrence of v in the sequence of concatenated adjacency
        // lists, which is represented by the wavelet tree.
        let p = self.tree.select(&v, i as u64);

        // Then we need to find the (p+1)-th occurence of a '0' in the bitmap. This corresponds to the
        // position of v in the bitmap representing the concatenated adjacency lists.
        let p0 = p.and_then(|x| self.bitmap.select_0(x+1));

        // Now we can count the number of '1's up until position p0 in order to find out which
        // adjacency list the i-th occurence of v is in.
        // Remember that we need to subtract 1 because we are returning an index.
        p0.and_then(|x| self.bitmap.rank_1(x)).map(|x| (x-1) as usize)
    }

    /// Returns a vector containing all indices of the neighbors (successors) of the node given by
    /// index v. If v has no neighbors, an empty vector is returned.
    pub fn neighbors(&self, v: usize) -> Vec<usize> {
        // Find the start and end of the adjacency list in the sequence (see access_neighbor for details)
        let l = self.bitmap.select_1((v+1) as u64)
            .map(|x| x - v as u64);

        let m = self.bitmap.select_1((v+2) as u64)
            .map(|l| l - (v+1) as u64)
            .unwrap_or(self.sequence_length());

        // Compute how many neighbors v has.
        let num_neighbors = l.map(|l| (m-l) as usize).unwrap_or(0);

        // Fill a vector with all indices of v's neighbors.
        let mut neighbors : Vec<usize> = Vec::with_capacity(num_neighbors);

        for i in 1..num_neighbors+1 {
            neighbors.push(self.access_neighbor(v, i).unwrap());
        }

        neighbors
    }

    /// Returns a vector containing all indices of the reverse neighbors (predecessors) of the node given by
    /// index v. If v has no reverse neighbors, an empty vector is returned.
    pub fn reverse_neighbors(&self, v: usize) -> Vec<usize> {
        // Get the number of occurences of v in the sequence. This is equal to the amount of v's reverse neighbors.
        let num_v = self.tree.rank(&v, self.sequence_length()-1).unwrap_or(0) as usize;

        // Fill a vector with all indices of v's reverse neighbors.
        let mut reverse_neighbors : Vec<usize> = Vec::with_capacity(num_v);

        for i in 1..num_v+1 {
            reverse_neighbors.push(self.access_reverse_neighbor(v, i).unwrap());
        }

        reverse_neighbors
    }

    /// Returns whether an edge exists between the nodes given by the indices 'from' and 'to'.
    pub fn edge_exists(&self, from: usize, to: usize) -> bool {
        self.neighbors(from).contains(&to)
    }

    fn sequence_length(&self) -> u64 {
        let bitmap_length = self.bitmap.bits().len();
        bitmap_length - self.bitmap.rank_1(bitmap_length-1).unwrap_or(0)
    }

}


#[cfg(test)]
mod tests {
    use super::*;
    use super::super::*;
    use bv::bit_vec;

    fn example_graph() -> Graph<&'static str, &'static str> {
        let mut graph = Graph::<&str, &str>::new();

        let v1 = graph.add_node("v1");
        let v2 = graph.add_node("v2");
        let v3 = graph.add_node("v3");
        let v4 = graph.add_node("v4");
        let v5 = graph.add_node("v5");
        let v6 = graph.add_node("v6");

        graph.extend_with_edges(&[
            (v1,v2), (v1, v4),
            (v2,v1), (v2, v4), (v2, v3),
            (v4,v3),
            (v5,v1), (v5,v4)
        ]);

        graph
    }

    #[test]
    fn test_from_graph() {
        let gwt = GraphWaveletTree::from_graph(example_graph());

        let expected_bits = bit_vec![
            true, false, false,
            true, false, false, false,
            true,
            true, false,
            true, false, false,
            true
        ];

        assert_eq!(expected_bits, gwt.bitmap.bits());

        assert_eq!(Some(&1), gwt.tree.access(0));
        assert_eq!(Some(&3), gwt.tree.access(1));
        assert_eq!(Some(&0), gwt.tree.access(2));
        assert_eq!(Some(&3), gwt.tree.access(3));
        assert_eq!(Some(&2), gwt.tree.access(4));
        assert_eq!(Some(&2), gwt.tree.access(5));
        assert_eq!(Some(&0), gwt.tree.access(6));
        assert_eq!(Some(&3), gwt.tree.access(7));
    }

    #[test]
    fn test_access_neighbor() {
        let gwt = GraphWaveletTree::from_graph(example_graph());

        assert_eq!(Some(1), gwt.access_neighbor(0, 1));
        assert_eq!(Some(3), gwt.access_neighbor(0, 2));
        assert_eq!(None, gwt.access_neighbor(0, 3));

        assert_eq!(Some(0), gwt.access_neighbor(1, 1));
        assert_eq!(Some(3), gwt.access_neighbor(1, 2));
        assert_eq!(Some(2), gwt.access_neighbor(1, 3));
        assert_eq!(None, gwt.access_neighbor(1, 4));

        assert_eq!(None, gwt.access_neighbor(2, 1));

        assert_eq!(Some(2), gwt.access_neighbor(3, 1));
        assert_eq!(None, gwt.access_neighbor(3, 2));

        assert_eq!(Some(0), gwt.access_neighbor(4, 1));
        assert_eq!(Some(3), gwt.access_neighbor(4, 2));
        assert_eq!(None, gwt.access_neighbor(4, 3));

        assert_eq!(None, gwt.access_neighbor(5, 1));

    }

    #[test]
    fn test_access_reverse_neighbor() {
        let gwt = GraphWaveletTree::from_graph(example_graph());

        assert_eq!(Some(1), gwt.access_reverse_neighbor(0, 1));
        assert_eq!(Some(4), gwt.access_reverse_neighbor(0, 2));
        assert_eq!(None, gwt.access_reverse_neighbor(0, 3));

        assert_eq!(Some(0), gwt.access_reverse_neighbor(1, 1));
        assert_eq!(None, gwt.access_reverse_neighbor(1, 2));

        assert_eq!(Some(1), gwt.access_reverse_neighbor(2, 1));
        assert_eq!(Some(3), gwt.access_reverse_neighbor(2, 2));
        assert_eq!(None, gwt.access_reverse_neighbor(2, 3));

        assert_eq!(Some(0), gwt.access_reverse_neighbor(3, 1));
        assert_eq!(Some(1), gwt.access_reverse_neighbor(3, 2));
        assert_eq!(Some(4), gwt.access_reverse_neighbor(3, 3));
        assert_eq!(None, gwt.access_reverse_neighbor(3, 4));

        assert_eq!(None, gwt.access_reverse_neighbor(4, 1));

        assert_eq!(None, gwt.access_reverse_neighbor(5, 1));

    }

    #[test]
    fn test_neighbors() {
        let gwt = GraphWaveletTree::from_graph(example_graph());

        let mut n = gwt.neighbors(0); n.sort();
        assert_eq!(vec![1, 3], n);

        let mut n = gwt.neighbors(1); n.sort();
        assert_eq!(vec![0, 2, 3], n);

        let mut n = gwt.neighbors(2); n.sort();
        assert_eq!(Vec::<usize>::new(), n);

        let mut n = gwt.neighbors(3); n.sort();
        assert_eq!(vec![2], n);

        let mut n = gwt.neighbors(4); n.sort();
        assert_eq!(vec![0, 3], n);

        let mut n = gwt.neighbors(5); n.sort();
        assert_eq!(Vec::<usize>::new(), n);
    }

    #[test]
    fn test_reverse_neighbors() {
        let gwt = GraphWaveletTree::from_graph(example_graph());

        print!("Sequence: ");
        for i in 0..8 {
            print!("{:?}", gwt.tree.access(i).unwrap());
        }
        println!("\nTree: {:?}", gwt.tree);


        let mut n = gwt.reverse_neighbors(0); n.sort();
        assert_eq!(vec![1, 4], n);

        let mut n = gwt.reverse_neighbors(1); n.sort();
        assert_eq!(vec![0], n);

        let mut n = gwt.reverse_neighbors(2); n.sort();
        assert_eq!(vec![1, 3], n);

        let mut n = gwt.reverse_neighbors(3); n.sort();
        assert_eq!(vec![0, 1, 4], n);

        let mut n = gwt.reverse_neighbors(4); n.sort();
        assert_eq!(Vec::<usize>::new(), n);

        let mut n = gwt.reverse_neighbors(5); n.sort();
        assert_eq!(Vec::<usize>::new(), n);
    }

    #[test]
    fn test_edge_exists() {
        let gwt = GraphWaveletTree::from_graph(example_graph());

        assert_eq!(false, gwt.edge_exists(0, 0));
        assert_eq!(true, gwt.edge_exists(0, 1));
        assert_eq!(false, gwt.edge_exists(0, 2));
        assert_eq!(true, gwt.edge_exists(0, 3));
        assert_eq!(false, gwt.edge_exists(0, 4));
        assert_eq!(false, gwt.edge_exists(0, 5));

        assert_eq!(true, gwt.edge_exists(1, 0));
        assert_eq!(false, gwt.edge_exists(1, 1));
        assert_eq!(true, gwt.edge_exists(1, 2));
        assert_eq!(true, gwt.edge_exists(1, 3));
        assert_eq!(false, gwt.edge_exists(1, 4));
        assert_eq!(false, gwt.edge_exists(1, 5));

        assert_eq!(false, gwt.edge_exists(2, 0));
        assert_eq!(false, gwt.edge_exists(2, 1));
        assert_eq!(false, gwt.edge_exists(2, 2));
        assert_eq!(false, gwt.edge_exists(2, 3));
        assert_eq!(false, gwt.edge_exists(2, 4));
        assert_eq!(false, gwt.edge_exists(2, 5));

        assert_eq!(false, gwt.edge_exists(3, 0));
        assert_eq!(false, gwt.edge_exists(3, 1));
        assert_eq!(true, gwt.edge_exists(3, 2));
        assert_eq!(false, gwt.edge_exists(3, 3));
        assert_eq!(false, gwt.edge_exists(3, 4));
        assert_eq!(false, gwt.edge_exists(3, 5));

        assert_eq!(true, gwt.edge_exists(4, 0));
        assert_eq!(false, gwt.edge_exists(4, 1));
        assert_eq!(false, gwt.edge_exists(4, 2));
        assert_eq!(true, gwt.edge_exists(4, 3));
        assert_eq!(false, gwt.edge_exists(4, 4));
        assert_eq!(false, gwt.edge_exists(4, 5));

        assert_eq!(false, gwt.edge_exists(5, 0));
        assert_eq!(false, gwt.edge_exists(5, 1));
        assert_eq!(false, gwt.edge_exists(5, 2));
        assert_eq!(false, gwt.edge_exists(5, 3));
        assert_eq!(false, gwt.edge_exists(5, 4));
        assert_eq!(false, gwt.edge_exists(5, 5));

    }

}