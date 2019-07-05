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

        // Return a new PointerWaveletTree containing the alphabet and the root of the tree
        PointerWaveletTree {
            alphabet,
            root
        }
    }

    
    /// Access element at index i in the sequence
    pub fn access(&self, i: usize) -> Option<&T> {
        // If the given index is larger than the size of the bitmap (i.e. if it is larger than the amout of symbols in
        // the sequence), no symbol can be returned.
        // Otherwise the access() function of WaveletTreeNode is called.
        if i as u64 >= self.root.bitmap.bits().len() {
            None
        } else {
            // Find out the the alphabet index of the symbol, which is at position i in the sequence
            let index_in_alphabet = self.root.access(i as u64, 0, self.alphabet.len());
            // If some index is returned, return the corresponding symbol
            index_in_alphabet.map(|index| &self.alphabet[index])
        }
    }
      // returns the position of the  the i-th accurred symbol q in a sequence
     pub fn select(&self, q: T, i: u64) -> Option<u64> {
        // If the given index is larger than the size of the bitmap (i.e. if it is larger than the amout of symbols in
        // the sequence), no result can be returned.
        // Otherwise the select() function of WaveletTreeNode is being called    
        if i as u64 >= self.root.bitmap.bits().len() {
            None
        }

        else {
        // Compute select(q,i) on the root node (bitmap subrange [a, b), representing the whole alphabet)

             // Find the index of q in the alphabet
            let q_index : Option<usize> = self.alphabet.iter().position(|x| x == &q);
            
            self.root.select( q_index.unwrap() as u64, i, 0, self.alphabet.len() as usize)
           
          // If q could not be found in the alphabet, return None. Otherwise, return select_c(i).
         //q_index.and_then(|qi| self.root.select(qi as u64, i, 0, self.alphabet.len() as usize))

        }
     }


    /// Returns the number of appearances of symbol c in the sequence up until position i
    pub fn rank(&self, c: T, i: u64) -> Option<u64> {
        // If the given index is larger than the size of the bitmap (i.e. if it is larger than the amout of symbols in
        // the sequence), no rank can be returned.
        // Otherwise the rank() function of WaveletTreeNode is called.
        if i as u64 >= self.root.bitmap.bits().len() {
            None
        } else {
            // Find the index of c in the alphabet
            let c_index : Option<usize> = self.alphabet.iter().position(|x| x == &c);
            // If c could not be found in the alphabet, return None. Otherwise, return rank_c(i).
            c_index.and_then(|ci| self.root.rank(ci as u64, i, 0, self.alphabet.len()))
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
     // If the bitmap  of q contains a 1 at position i, look in the right subtree.
        // Otherwise look in the left subtree.
      
      
    fn select (&self, q_index: u64, i: u64, a: usize, b: usize) -> Option<u64> {
     //If the leaf is the leftchild of its parent v,then the position i´ corresponding to 
     //i at v is the i-th occurrence of a 0 initsbitmap Bv.
     let btmp = &self.bitmap; 
     let middle = (a+b)/2;
     
              if q_index > middle as u64 {
                   self.right_child.as_ref().and_then(|child| child.select(q_index, btmp.select_1(i).unwrap() - 1, middle as usize, b))
              }
             
              else if q_index < middle as u64 { 
                  self.left_child.as_ref().and_then(|child| child.select(q_index, btmp.select_0(i).unwrap() - 1, a, middle as usize)) 
              }
             else{ Some(middle as u64)  }
           
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

    #[test]
    fn test_rank_first_appearance() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerWaveletTree::from_sequence(sequence);

        // returns 0 before first appearance
        assert_eq!(Some(0), tree.rank('l', 0));
        assert_eq!(Some(0), tree.rank('b', 2));
        assert_eq!(Some(0), tree.rank('r', 4));
        assert_eq!(Some(0), tree.rank(' ', 5));
        assert_eq!(Some(0), tree.rank('d', 17));

        // returns 1 at exact position of first appearance
        assert_eq!(Some(1), tree.rank('a', 0));
        assert_eq!(Some(1), tree.rank('l', 1));
        assert_eq!(Some(1), tree.rank('b', 3));
        assert_eq!(Some(1), tree.rank('r', 5));
        assert_eq!(Some(1), tree.rank(' ', 6));
        assert_eq!(Some(1), tree.rank('d', 18));
    }

    #[test]
    #[ignore]
    fn test_rank_all_appearances_of_one_symbol() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerWaveletTree::from_sequence(sequence);
        
        assert_eq!(Some(1), tree.rank('a', 0));
        assert_eq!(Some(2), tree.rank('a', 2));
        assert_eq!(Some(3), tree.rank('a', 4));
        assert_eq!(Some(4), tree.rank('a', 7));
        assert_eq!(Some(5), tree.rank('a', 10));
        assert_eq!(Some(6), tree.rank('a', 12));
        assert_eq!(Some(7), tree.rank('a', 14));
        assert_eq!(Some(8), tree.rank('a', 16));
        assert_eq!(Some(9), tree.rank('a', 19));
    }

    #[test]
    fn test_rank_invalid_argument() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerWaveletTree::from_sequence(sequence);
        
        // out of range
        assert_eq!(None, tree.rank('a', 20));
        // unknown symbol
        assert_eq!(None, tree.rank('x', 19));
    }


     #[test]

    fn test_pointer_select_first_appearance_of_each_symbol() {

        let text = "alabar a la alabarda";

        let sequence : &Vec<char> = &text.chars().collect();

        let tree = PointerWaveletTree::from_sequence(sequence);



        assert_eq!(Some(0), tree.select('a', 1));

        assert_eq!(Some(1), tree.select('l', 1));

        assert_eq!(Some(3), tree.select('b', 1));

        assert_eq!(Some(5), tree.select('r', 1));

        assert_eq!(Some(6), tree.select(' ', 1));

        assert_eq!(Some(18), tree.select('d', 1));

    }



    #[test]

    fn test_pointer_select_all_appearances_of_one_symbol() {

        let text = "alabar a la alabarda";

        let sequence : &Vec<char> = &text.chars().collect();

        let tree = PointerWaveletTree::from_sequence(sequence);



        assert_eq!(Some(0), tree.select('a', 1));
        assert_eq!(Some(2), tree.select('a', 2));
        assert_eq!(Some(4), tree.select('a', 3));
        assert_eq!(Some(7), tree.select('a', 4));
        assert_eq!(Some(10), tree.select('a', 5));
        assert_eq!(Some(12), tree.select('a', 6));
        assert_eq!(Some(14), tree.select('a', 7));
        assert_eq!(Some(16), tree.select('a', 8));
        assert_eq!(Some(19), tree.select('a', 9));

    }



    #[test]

    fn test_pointer_select_invalid_argument() {

        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerWaveletTree::from_sequence(sequence);



        // Out of range

        assert_eq!(None, tree.select('a', 20));

        // Symbol not in alphabet

        assert_eq!(None, tree.select('x', 1));

        // i too high

        assert_eq!(None, tree.select('d', 2));

        assert_eq!(None, tree.select('l', 4));

        assert_eq!(None, tree.select('r', 3));

    }



    #[test]

    fn test_pointer_select_all_valid_arguments() {

        let text = "The quick brown fox jumps over a lazy dog";

        let sequence : &Vec<char> = &text.chars().collect();

        let tree = PointerWaveletTree::from_sequence(sequence);



        let mut alphabet = Vec::new();

        for symbol in sequence.iter() {

            if !alphabet.contains(symbol) {

                alphabet.push(symbol.clone());

            }

        }

        alphabet.sort_by(|x,y| x.partial_cmp(y).unwrap());



        for i in 0..alphabet.len() {

            let symbol = alphabet[i].clone();

            let mut index = 0;

            for j in 1..sequence.iter().filter(|&n| *n == symbol).count()+1 {

                for k in index..sequence.len() {if sequence[k] != symbol {index += 1;} else {break;}}

                assert_eq!(Option::Some(index as u64), tree.select(symbol, j as u64));

                index += 1;

            };

        }

    }



    #[test]

    fn test_pointer_select_ascending_alphabet_as_sequence() {

        let text = "abcdefghijklmnopqrstuvwxyz";

        let sequence : &Vec<char> = &text.chars().collect();

        let tree = PointerWaveletTree::from_sequence(sequence);



        let mut pos = 0;

        for symbol in tree.alphabet.clone().into_iter() {

            assert_eq!(Some(pos), tree.select(symbol, 1));

            pos = pos + 1;

        }

    }



    #[test]

    fn test_pointer_select_descending_alphabet_as_sequence() {

        let text = "zyxwvutsrqponmlkjihgfedcba";

        let sequence : &Vec<char> = &text.chars().collect();

        let tree = PointerWaveletTree::from_sequence(sequence);



        let mut pos = 0;

        for symbol in tree.alphabet.clone().into_iter().rev() {

            assert_eq!(Some(pos), tree.select(symbol, 1));

            pos = pos + 1;

        }

    }

}