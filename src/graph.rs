use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;
use petgraph::graph::Graph;
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize)]
pub struct LabelledGraphWaveletTree<T> {
    graphTree: GraphWaveletTree,
    labels: Vec<T>
}

#[derive(Serialize, Deserialize)]
pub struct GraphWaveletTree {
    /// The wavelet tree represents the concatenated adjacency lists of the graph
    tree: WaveletTree<usize>,
    /// The bitmap marks the beginnings of the adjacency lists
    bitmap: RankSelect
}

impl GraphWaveletTree {
    pub fn from_graph(graph: Graph) -> Self {
        unimplemented!();
    }

    pub fn from_adjacency_lists(nodes: NodeIndices, lists: Vec<Neighbours>) -> Self {
        unimplemented!();
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


}