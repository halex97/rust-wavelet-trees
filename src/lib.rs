use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;

pub struct PointerWaveletTree<T: PartialOrd + Clone> {
    bitmap: Option<RankSelect>,
    label: Option<T>,
    left_child: Option<Box<PointerWaveletTree<T>>>,
    right_child: Option<Box<PointerWaveletTree<T>>>
}

impl <T: PartialOrd + Clone> PointerWaveletTree<T> {

    /// Create a new wavelet tree
    pub fn new(sequence: Vec<T>, alphabet: Vec<T>) -> PointerWaveletTree<T> {
        if sequence.is_empty() {
            panic!("Wavelet trees cannot be created from empty sequences.");
        }
        if alphabet.is_empty() {
            panic!("Wavelet trees cannot be created from an empty alphabet.");
        }

        // If there is only one symbol in the alphabet of the sequence (i.e. all symbols in the sequence are equal),
        // a leaf should be created.
        // Otherwise we create a bitmap for the sequence and use it to create children.

        if alphabet.len() == 1 {
            // Return a leaf (= a wavelet tree with a label and without children)
            PointerWaveletTree {
                bitmap: None,
                label: Some(alphabet[0].clone()),
                left_child: None,
                right_child: None
            }
        } else {
            // Create the bitmap.
            let center_of_alphabet = alphabet.len()/2;
            let bitmap = create_bitmap(&sequence, &alphabet[center_of_alphabet]);

            // Now we can split up the sequence of this wavelet tree in order to create the left and the right child.
            let mut left_sequence : Vec<T> = Vec::new();
            let mut right_sequence : Vec<T> = Vec::new();

            for i in 0..sequence.len() {
                // If the bitmap contains a 1 at position i, the symbol of the sequence at position i should be
                // in the right child's sequence.
                // Otherwise it should be a part of the left child's sequence.
                if bitmap.get(i as u64) {
                    right_sequence.push(sequence[i].clone());
                } else {
                    left_sequence.push(sequence[i].clone());
                }
            }

            // We can now recursively create the children of this wavelet tree.
            // A child only needs to be created if it's corresponding sequence is not empty.
            let left_child;
            let right_child;

            if left_sequence.is_empty() {
                left_child = None;
            } else {
                left_child = Some(Box::new(PointerWaveletTree::new(left_sequence, alphabet[..center_of_alphabet].to_vec())));
            }

            if right_sequence.is_empty() {
                right_child = None;
            } else {
                right_child = Some(Box::new(PointerWaveletTree::new(right_sequence, alphabet[center_of_alphabet..].to_vec())));
            }

            // Putting everything together, create the wavelet tree that contains the created children
            PointerWaveletTree {
                bitmap: Some(bitmap),
                label: None,
                left_child,
                right_child
            }
        }
    }

    pub fn new_pointer(sequence: &Vec<T>, alphabet: &[T]) -> Option<Box<PointerWaveletTree<T>>> {

        if sequence.is_empty() || (alphabet.len() <= 1) {
            return Option::None;
        }
        
        let center_of_alphabet = alphabet.len()/2;

        Option::Some(Box::new(
            PointerWaveletTree {
                bitmap: Some(create_sequence_bitmap(sequence, alphabet, &alphabet[center_of_alphabet])),
                label: None,
                left_child: PointerWaveletTree::new_pointer(sequence, &alphabet[..center_of_alphabet]),
                right_child: PointerWaveletTree::new_pointer(sequence, &alphabet[center_of_alphabet..])
            }
        ))
    }

    /// Access element at index i
    pub fn access(&self, i: u64) -> Option<&T> {
        if self.is_leaf() {
            self.label.as_ref()
        } else {
            self.bitmap.as_ref().and_then(|bm: &RankSelect| -> Option<&T> {
                if bm.get(i) {
                    self.right_child.as_ref().and_then(|child| child.access(bm.rank_1(i).unwrap() - 1))
                } else {
                    self.left_child.as_ref().and_then(|child| child.access(bm.rank_0(i).unwrap() - 1))
                }
            })
        }
    }

    fn is_leaf(&self) -> bool {
        self.left_child.is_none() && self.right_child.is_none()
    }

}

fn create_bitmap<T: PartialOrd>(sequence: &Vec<T>, mid_symbol: &T) -> RankSelect {
    let mut bits : BitVec<u8> = BitVec::new();

    for symbol in sequence.iter() {
        bits.push(symbol >= mid_symbol);
    }

    RankSelect::new(bits, 1)
}

fn create_sequence_bitmap<T: PartialOrd>(sequence: &Vec<T>, alphabet: &[T], mid_symbol: &T) -> RankSelect {
    let mut bits : BitVec<u8> = BitVec::new();

    for symbol in sequence.iter() {
        if alphabet.contains(symbol) {
            bits.push(symbol >= mid_symbol);
        }
    }

    RankSelect::new(bits, 1)
}


#[cfg(test)]
mod tests {
    use super::*;
    use bv::bit_vec;

    #[test]
    fn test_pointer_wavelet_tree_new_leaf() {
        let pwt : PointerWaveletTree<char> = PointerWaveletTree::new(vec!['a','a','a'], vec!['a']);

        assert!(pwt.is_leaf());
        assert_eq!(pwt.label, Some('a'));
    }

    #[test]
    fn test_create_bitmap() {
        let text = "alabar a la alabarda";
        let sequence : Vec<char> = text.chars().collect();
        let ref_symbol = 'd';

        let bitmap = create_bitmap(&sequence, &ref_symbol);

        let expected_bits = bit_vec![
            false, true, false, false, 
            false, true, false, false, 
            false, true, false, false, 
            false, true, false, false, 
            false, true, true, false];
        assert_eq!(bitmap.bits(), &expected_bits);
    }

    #[test]
    fn test_pointer_wavelet_tree_new_has_children() {
        let text = "alabar a la alabarda";
        let alphabet = vec![' ','a','b','d','l','r'];
        let pwt = PointerWaveletTree::new(text.chars().collect(), alphabet);

        // The wavelet tree should not be a leaf and, therefore, have no label.
        assert!(!pwt.is_leaf());
        assert_eq!(pwt.label, None);

        // The root of the wavelet tree should have a bitmap.
        assert!(pwt.bitmap.is_some());

        // The root of the wavelet tree should have two children.
        assert!(pwt.left_child.is_some());
        assert!(pwt.right_child.is_some());
    }

    #[test]
    fn test_access_leaf() {
        let pwt : PointerWaveletTree<char> = PointerWaveletTree::new(vec!['a','a','a'], vec!['a']);
        assert_eq!(pwt.access(0), Some('a').as_ref());
        assert_eq!(pwt.access(1), Some('a').as_ref());
        assert_eq!(pwt.access(2), Some('a').as_ref());
    }

    #[test]
    fn test_access() {
        let text = "alabar a la alabarda";
        let alphabet = vec![' ','a','b','d','l','r'];
        let pwt = PointerWaveletTree::new(text.chars().collect(), alphabet);

        assert_eq!(pwt.access(0), Some('a').as_ref());
        assert_eq!(pwt.access(1), Some('l').as_ref());
        assert_eq!(pwt.access(2), Some('a').as_ref());
        assert_eq!(pwt.access(3), Some('b').as_ref());
        assert_eq!(pwt.access(4), Some('a').as_ref());
        assert_eq!(pwt.access(5), Some('r').as_ref());
        assert_eq!(pwt.access(6), Some(' ').as_ref());
        assert_eq!(pwt.access(7), Some('a').as_ref());
        assert_eq!(pwt.access(8), Some(' ').as_ref());
        assert_eq!(pwt.access(9), Some('l').as_ref());
        assert_eq!(pwt.access(10), Some('a').as_ref());
        assert_eq!(pwt.access(11), Some(' ').as_ref());
        assert_eq!(pwt.access(12), Some('a').as_ref());
        assert_eq!(pwt.access(13), Some('l').as_ref());
        assert_eq!(pwt.access(14), Some('a').as_ref());
        assert_eq!(pwt.access(15), Some('b').as_ref());
        assert_eq!(pwt.access(16), Some('a').as_ref());
        assert_eq!(pwt.access(17), Some('r').as_ref());
        assert_eq!(pwt.access(18), Some('d').as_ref());
        assert_eq!(pwt.access(19), Some('a').as_ref());
    }
}
