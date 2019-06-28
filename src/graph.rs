use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;
use petgraph::graph::Graph;

pub struct LabelledGraphWaveletTree {
    graphTree: WaveletGraph,
    labels: Vec<T>
}

pub struct GraphWaveletTree {
    /// The wavelet tree represents the concatenated adjacency lists of the graph
    tree: WaveletTree<u64>,
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

    pub fn predecessor(&self, index: u64) -> u64 {
        unimplemented!();
    }

    pub fn successor(&self, index: u64) -> u64 {
        unimplemented!();
    }

    pub fn edge_exists(&self, predecessor: u64, successor: u64) -> bool {
        unimplemented!();
    }

}