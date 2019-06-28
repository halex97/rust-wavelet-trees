pub mod pointers;
pub mod pointerless;

pub trait WaveletTree<T> {
    fn from_slice(sequence: &[T]) -> Self;

    fn rank(&self, c: T, i: u64) -> Option<u64>;
    fn select(&self, c: T, i: u64) -> Option<u64>;
    fn access(&self, i: u64) -> Option<&T>;
}
