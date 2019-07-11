pub mod pointers;
pub mod pointerless;
pub mod graph;

pub trait WaveletTree<T> {
    fn from_iterator(sequence: &mut dyn std::iter::Iterator<Item=T>) -> Self;
    fn from_slice(sequence: &[T]) -> Self;

    fn rank(&self, c: &T, i: u64) -> Option<u64>;
    fn select(&self, c: &T, i: u64) -> Option<u64>;
    fn access(&self, i: u64) -> Option<&T>;
}
