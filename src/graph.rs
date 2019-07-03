use bio::data_structures::rank_select::RankSelect;
use bv::{BitVec, BitsMut};
use petgraph::EdgeType;
use petgraph::graph::{Graph, NodeIndices, NodeIndex, Neighbors, DefaultIx, IndexType};
use serde::{Serialize, Deserialize};

use crate::WaveletTree;
use crate::pointerless::PointerlessWaveletTree;
use core::borrow::BorrowMut;

#[derive(Serialize, Deserialize)]
pub struct LabelledGraphWaveletTree<T> {
    graphTree: GraphWaveletTree,
    labels: Vec<T>
}

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

            // Iterate over all neighbors (i.e. adjacent nodes that appear in the list)
            for n in neighbors {
                // Add the neighbor's index to the sequence
                sequence.push(n.index());

                // '0's in the bit vector represent node indices in the adjacency lists
                lists_bv.push(false);
            }
        }

        // Return a pointerless wavelet tree representing the sequence together with the bitmap.
        GraphWaveletTree {
            tree: PointerlessWaveletTree::from_slice(sequence.as_slice()),
            bitmap: RankSelect::new(lists_bv, 1)
        }
    }

    /// Access the i-th neighbour (i.e. successor) of node v
    pub fn access_neighbor(&self, v: usize, i: usize) -> Option<usize> {
        // First we need to find the start of the adjacency list.
        // Then we need to take into account that the bitmap is longer than the actual sequence
        // because of the additional '1's.
        // Finally we can add the index i to find the correct neighbor's index in the sequence of
        // concatenated adjacency lists.
        let n = self.bitmap.select_1(v as u64)
            .map(|l| l - (v-1) as u64)
            .map(|l| l + i as u64);

        // Now we can look for the correct node index in the wavelet tree.
        n.and_then(|x| self.tree.access(x)).map(|x| x.clone())
    }

    /// Access the i-th reverse neighbor (i.e. predecessor) of node v
    pub fn access_reverse_neighbor(&self, v: usize, i: usize) -> Option<usize> {
        // First we need to find the i-th occurrence of v in the sequence of concatenated adjacency
        // lists, which is represented by the wavelet tree.
        let p = self.tree.select(&v, i as u64);

        // Then we need to find the p-th occurence of a '0' in the bitmap. This corresponds ot the
        // position of v in the bitmap representing the concatenated adjacency lists.
        let p0 = p.and_then(|x| self.bitmap.select_0(x));

        // Now we can count the number of '1's up until position p0 in order to find out which
        // adjacency list the i-th occurence of v is in.
        p0.map(|x| self.bitmap.rank_1(x) as usize)
    }

    pub fn predecessor(&self, index: usize) -> NodeIndex {
        unimplemented!();
    }

    pub fn successor(&self, index: usize) -> NodeIndex {
        unimplemented!();
    }

    pub fn edge_exists(&self, predecessorIndex: usize, successorIndex: usize) -> bool {
        unimplemented!();
    }

}


#[cfg(test)]
mod tests {
    use super::*;
    use super::super::*;

    #[test]
    fn test_from_graph() {
        unimplemented!();
    }

}