use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;

pub struct PointerWaveletTree<T: PartialOrd + Clone> {
    alphabet: Vec<T>,
    root: WaveletTreeNode
}

impl <T: PartialOrd + Clone> PointerWaveletTree<T> {
    pub fn from_sequence(sequence: &[T]) -> Self {
        // Create a vector for storing the alphabet of the sequence
        let mut alphabet = Vec::new();

        // Add all symbols from the sequence to the new alphabet vector
        for symbol in sequence.iter() {
            if !alphabet.contains(symbol) {
                alphabet.push(symbol.clone());
            }
        }

        // Sort the alphabet
        alphabet.sort_by(|x,y| x.partial_cmp(y).unwrap());
        
        // Create the root of the wavelet tree (which recursively creates all other nodes)
        let root = WaveletTreeNode::new(sequence, &alphabet);

        PointerWaveletTree {
            alphabet,
            root
        }
    }

    
    /// Access element at index i in the sequence
    pub fn access(&self, i: usize) -> Option<&T> {
        if i as u64 >= self.root.bitmap.bits().len() {
            None
        } else {
            let index_in_alphabet = self.root.access(i as u64, 0, self.alphabet.len());

            index_in_alphabet.map(|index| &self.alphabet[index])
        }
    }
}

struct WaveletTreeNode {
    bitmap: RankSelect,
    left_child: Option<Box<WaveletTreeNode>>,
    right_child: Option<Box<WaveletTreeNode>>
}

impl WaveletTreeNode {
    fn new <T: PartialOrd + Clone> (sequence: &[T], alphabet: &[T]) -> Self {
        if sequence.is_empty() {
            panic!("Wavelet trees cannot be created from empty sequences.");
        }
        if alphabet.is_empty() {
            panic!("Wavelet trees cannot be created from an empty alphabet.");
        }

        // Create the bitmap.
        let bitmap = create_bitmap(sequence, alphabet);

        let center_of_alphabet = alphabet.len()/2;
        let left_alphabet = &alphabet[..center_of_alphabet];
        let right_alphabet = &alphabet[center_of_alphabet..];
        
        let left_child = Self::create_boxed_option(sequence, left_alphabet);
        let right_child = Self::create_boxed_option(sequence, right_alphabet);

        WaveletTreeNode {
            bitmap,
            left_child,
            right_child
        }
    }

    fn create_boxed_option <T: PartialOrd + Clone> (sequence: &[T], alphabet: &[T]) -> Option<Box<Self>> {
        if sequence.len() <= 1 || alphabet.len() <= 1 {
            None
        } else {
            Some(Box::new(WaveletTreeNode::new(sequence, alphabet)))
        }
    }

    /// returns the alphabet index of the symbol at position i in the sequence
    /// [a,b) is the subrange of the alphabet that the current node represents
    fn access(&self, i: u64, a: usize, b: usize) -> Option<usize> {
        let bm = &self.bitmap;
        let center = (a+b)/2;

        // If the bitmap contains a 1 at position i, look in the right subtree.
        // Otherwise look in the left subtree.
        if bm.get(i) {
            if b-center <= 1 {
                Some(center)
            } else {
                self.right_child.as_ref().and_then(|child| child.access(bm.rank_1(i).unwrap() - 1, center, b))
            }
        } else {
            if center-a <= 1 {
                Some(a)
            } else {
                self.left_child.as_ref().and_then(|child| child.access(bm.rank_0(i).unwrap() - 1, a, center))
            }
        }
    }
}


fn create_bitmap<T: PartialOrd>(sequence: &[T], alphabet: &[T]) -> RankSelect {
    let ref_symbol = &alphabet[alphabet.len()/2];
    let mut bits : BitVec<u8> = BitVec::new();

    for symbol in sequence.iter() {
        if alphabet.contains(symbol) {
            bits.push(symbol >= ref_symbol);
        }
    }

    RankSelect::new(bits, 1)
}


#[cfg(test)]
mod tests {
    use super::*;
    use bv::bit_vec;

    #[test]
    fn test_create_bitmap() {
        let text = "alabar a la alabarda";
        let sequence : Vec<char> = text.chars().collect();

        let bitmap = create_bitmap(&sequence, &[' ','a','b','d','l','r']);

        let expected_bits = bit_vec![
            false, true, false, false, 
            false, true, false, false, 
            false, true, false, false, 
            false, true, false, false, 
            false, true, true, false];
        assert_eq!(bitmap.bits(), &expected_bits);
    }

    #[test]
    fn test_pointer_wavelet_tree_new() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let pwt = PointerWaveletTree::from_sequence(sequence);

        // The correct alphabet should automatically be created
        assert_eq!(pwt.alphabet, vec![' ','a','b','d','l','r']);

        // The root of the wavelet tree should have two children.
        assert!(pwt.root.left_child.is_some());
        assert!(pwt.root.right_child.is_some());
    }

    // #[test]
    // fn test_access_leaf() {
    //     let pwt : PointerWaveletTree<char> = PointerWaveletTree::new(vec!['a','a','a'], vec!['a']);
    //     assert_eq!(pwt.access(0), Some('a').as_ref());
    //     assert_eq!(pwt.access(1), Some('a').as_ref());
    //     assert_eq!(pwt.access(2), Some('a').as_ref());
    // }

    #[test]
    fn test_access_inside_range() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let pwt = PointerWaveletTree::from_sequence(sequence);

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

    #[test]
    fn test_access_out_of_range() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let pwt = PointerWaveletTree::from_sequence(sequence);

        assert_eq!(pwt.access(20), None);
    }
}
