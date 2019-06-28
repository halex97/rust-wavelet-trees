use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;
use petgraph::graph::Graph;

pub struct WaveletGraph<T: PartialOrd + Clone> {
    lists: WaveletTree<T>
    seperator: RankSelect
}

impl <T: PartialOrd + Clone> WaveletGraph<T> {
    pub fn from_graph(graph: Graph) -> Self {
        unimplemented!();
    }

    pub fn from_adjacency_list(nodes: &[T], edges: &[T, T]) -> Self {
        unimplemented!();
    }

    pub fn predecessor(&self, index: &T) -> Vec<T> {
        unimplemented!();
    }

    pub fn successor(&self, index: &T) -> Vec<T> {
        unimplemented!();
    }

    pub fn edge_exits(&self, predecessor: &T, successor: &T) -> bool {
        unimplemented!();
    }

}