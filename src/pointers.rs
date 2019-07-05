use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;
use rand::Rng;
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize)]
pub struct PointerWaveletTree<T: PartialOrd + Clone> {
    alphabet: Vec<T>,
    root: WaveletTreeNode
}

impl <T: PartialOrd + Clone> super::WaveletTree<T> for PointerWaveletTree<T> {

    fn from_iterator(sequence: &mut dyn std::iter::Iterator<Item=T>) -> Self {
        let vec_sequence: Vec<T> = sequence.collect();
        return Self::from_slice(&vec_sequence);
    }

    fn from_slice(sequence: &[T]) -> Self {
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

        // Return a new PointerWaveletTree containing the alphabet and the root of the tree
        PointerWaveletTree {
            alphabet,
            root
        }
    }
    
    /// Access element at index i in the sequence
    fn access(&self, i: u64) -> Option<&T> {
        // If the given index is larger than the size of the bitmap (i.e. if it is larger than the amout of symbols in
        // the sequence), no symbol can be returned.
        // Otherwise the access() function of WaveletTreeNode is called.
        if i >= self.root.bitmap.bits().len() {
            None
        } else {
            // Find out the the alphabet index of the symbol, which is at position i in the sequence
            let index_in_alphabet = self.root.access(i, 0, self.alphabet.len());
            // If some index is returned, return the corresponding symbol
            index_in_alphabet.map(|index| &self.alphabet[index])
        }
    }

    /// Returns the number of appearances of symbol c in the sequence up until position i
    fn rank(&self, c: &T, i: u64) -> Option<u64> {
        // If the given index is larger than the size of the bitmap (i.e. if it is larger than the amout of symbols in
        // the sequence), no rank can be returned.
        // Otherwise the rank() function of WaveletTreeNode is called.
        if i as u64 >= self.root.bitmap.bits().len() {
            None
        } else {
            // Find the index of c in the alphabet
            let c_index : Option<usize> = self.alphabet.iter().position(|x| x == c);
            // If c could not be found in the alphabet, return None. Otherwise, return rank_c(i).
            c_index.and_then(|ci| self.root.rank(ci as u64, i, 0, self.alphabet.len()))
        }
    }

    fn select(&self, c: &T, i: u64) -> Option<u64> {
        unimplemented!();
    }
}

#[derive(Serialize, Deserialize)]
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

        // Create the bitmap for this node.
        let bitmap = create_bitmap(sequence, alphabet);

        // Split the alphabet into a left part and a right part.
        let center_of_alphabet = alphabet.len()/2;
        let left_alphabet = &alphabet[..center_of_alphabet];
        let right_alphabet = &alphabet[center_of_alphabet..];

        // Create left and right children to represent the corresponding subranges of the alphabet.
        let left_child = Self::create_boxed_inner_node(sequence, left_alphabet);
        let right_child = Self::create_boxed_inner_node(sequence, right_alphabet);

        // Return a new WaveletTreeNode containing the created bitmap and children.
        WaveletTreeNode {
            bitmap,
            left_child,
            right_child
        }
    }

    fn create_boxed_inner_node <T: PartialOrd + Clone> (sequence: &[T], alphabet: &[T]) -> Option<Box<Self>> {
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

    /// returns the rank of symbol alphabet[c_index] in the sequence
    /// [a,b) is the subrange of the alphabet that the current node represents
    fn rank(&self, c_index: u64, i: u64, a: usize, b: usize) -> Option<u64> {
        let center = (a+b)/2;

        // First, we need to find out if the leaf representing c is in the left or right subtree.
        if c_index < center as u64 {
            // If the left subtree is the leaf, rank_0 of i is the final answer.
            // If rank_0(i)==0, there is no c up to position i (thus, rank_0 of i is the final answer).
            // Otherwise, we call rank on the left subtree with an adjusted index i.
            if center - a <= 1 || self.bitmap.rank_0(i) == Some(0) {
                self.bitmap.rank_0(i)
            } else {
                self.left_child.as_ref()
                               .and_then(|child| child.rank(c_index, self.bitmap.rank_0(i).unwrap()-1, a, center))
            }
        } else {
            // If the right subtree is the leaf, rank_1 of i is the final answer.
            // If rank_1(i)==0, there is no c up to position i (thus, rank_1 of i is the final answer).
            // Otherwise, we call rank on the right subtree with an adjusted index i.
            if b - center <= 1 || self.bitmap.rank_1(i) == Some(0) {
                self.bitmap.rank_1(i)
            } else {
                self.right_child.as_ref()
                                .and_then(|child| child.rank(c_index, self.bitmap.rank_1(i).unwrap()-1, center, b))
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
    use super::super::*;
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
    fn test_from_slice() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let pwt = PointerWaveletTree::from_slice(sequence);

        // The correct alphabet should automatically be created
        assert_eq!(pwt.alphabet, vec![' ','a','b','d','l','r']);

        // The root of the wavelet tree should have two children.
        assert!(pwt.root.left_child.is_some());
        assert!(pwt.root.right_child.is_some());
    }

    #[test]
    fn test_access_inside_range() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let pwt = PointerWaveletTree::from_slice(sequence);

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
        let pwt = PointerWaveletTree::from_slice(sequence);

        assert_eq!(pwt.access(20), None);
    }
  
    #[test]
    fn test_rank_first_appearance() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerWaveletTree::from_slice(sequence);

        // returns 0 before first appearan& ce
        assert_eq!(Some(0), tree.rank(& 'l', 0));
        assert_eq!(Some(0), tree.rank(& 'b', 2));
        assert_eq!(Some(0), tree.rank(& 'r', 4));
        assert_eq!(Some(0), tree.rank(& ' ', 5));
        assert_eq!(Some(0), tree.rank(& 'd', 17));

        // returns 1 at exact position of first appearance
        assert_eq!(Some(1), tree.rank(& 'a', 0));
        assert_eq!(Some(1), tree.rank(& 'l', 1));
        assert_eq!(Some(1), tree.rank(& 'b', 3));
        assert_eq!(Some(1), tree.rank(& 'r', 5));
        assert_eq!(Some(1), tree.rank(& ' ', 6));
        assert_eq!(Some(1), tree.rank(& 'd', 18));
    }

    #[test]
    fn test_rank_all_appearances_of_one_symbol() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerWaveletTree::from_slice(sequence);
        
        assert_eq!(Some(1), tree.rank(& 'a', 0));
        assert_eq!(Some(2), tree.rank(& 'a', 2));
        assert_eq!(Some(3), tree.rank(& 'a', 4));
        assert_eq!(Some(4), tree.rank(& 'a', 7));
        assert_eq!(Some(5), tree.rank(& 'a', 10));
        assert_eq!(Some(6), tree.rank(& 'a', 12));
        assert_eq!(Some(7), tree.rank(& 'a', 14));
        assert_eq!(Some(8), tree.rank(& 'a', 16));
        assert_eq!(Some(9), tree.rank(& 'a', 19));
    }

    #[test]
    fn test_rank_invalid_argument() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerWaveletTree::from_slice(sequence);
        
        // out of range
        assert_eq!(None, tree.rank(& 'a', 20));
        // unknown symbol
        assert_eq!(None, tree.rank(& 'x', 19));
    }

    #[test]
    #[ignore]
    fn test_randomized_access_rank_select() {
        // Time Nedded ~ 3-5 min
        let result = std::panic::catch_unwind(|| {PointerWaveletTree::from_slice(&Vec::<u64>::new());});
        assert!(result.is_err());
        let mut sequence: Vec<u64> = Vec::new();
        sequence.push(1);
        let tree = PointerWaveletTree::from_slice(&sequence);
        assert_eq!(Option::Some(&1), tree.access(0));
        assert_eq!(Option::None, tree.access(1));
        assert_eq!(Option::Some(1), tree.rank(&1, 0));
        assert_eq!(Option::None, tree.rank(&1, 1));
        //assert_eq!(Option::Some(0), tree.select(&1, 1));
        //assert_eq!(Option::None, tree.rank(&1, 2));

        let mut numbergen = rand::thread_rng();
        for size in 1..256 {
            // Build Alphabet
            let mut alphabet: Vec<u64> = Vec::new();
            for i in 0..size {alphabet.push(i);}

            // Build Number Vector
            let mut numbers: Vec<u64> = alphabet.to_vec();
            for _j in 0..(32+size/4-(size*size)/1024) {numbers.push(numbers[numbergen.gen_range(0, size) as usize]);}
            let tree = PointerWaveletTree::from_slice(&numbers);

            // Test Access Valid+Invalid
            for i in 0..numbers.len() {assert_eq!(Option::Some(&numbers[i]), tree.access(i as u64));}
            for i in numbers.len()..numbers.len()+256 {assert_eq!(Option::None, tree.access(i as u64));}

            for symbol in 0..size {
                let mut count = 0;
                // Test Rank Valid
                for index in 0..numbers.len() {
                    if numbers[index] == symbol {count += 1;}
                    assert_eq!(Option::Some(count as u64), tree.rank(&symbol, index as u64))
                }
                // Test Rank Invalid
                for index in numbers.len()..numbers.len()+500 {assert_eq!(Option::None, tree.rank(&symbol, index as u64))}
                // Test Select Valid
                let mut index = 0;
                for j in 1..count+1 {
                    for k in index..numbers.len() {if numbers[k] != symbol {index += 1;} else {break;}}
                    //assert_eq!(index as u64, tree.select(&symbol, j as u64).unwrap());
                    index += 1;
                };
            }
        }
    }
}