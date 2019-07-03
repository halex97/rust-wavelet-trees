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