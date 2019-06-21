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
}

pub struct PointerlessWaveletTree<T: PartialOrd + Clone> {
    alphabet: Vec<T>,
    bitmap: RankSelect
}

impl <T: PartialOrd + Clone> PointerlessWaveletTree<T> {

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

        let bitmap = Self::create_bits(sequence, &alphabet);

        // Return a new PointerWaveletTree containing the alphabet and the root of the tree
        PointerlessWaveletTree {
            alphabet,
            bitmap
        }
    }

    /// Create new Pointerless Wavelet Tree
    pub fn create_bits(sequence: &[T], alphabet: &[T]) -> RankSelect {
        if sequence.is_empty() {
            return RankSelect::new(BitVec::new(), 1);
        }
        else {
            let mut bits: BitVec<u8> = BitVec::new();

            // Calculates smallest total number d with alphabet.length <= 2^d
            let bound = Self::alphabet_bound(alphabet.len());

            // Defines how the Alphabet is partitioned
            let partition = Self::partition_alphabet(bound, alphabet.len());

            // index = sub-alphabet length in part of partition
            let mut index;
            // last = start of subsequence, next = end of subsequence, mid = midlle between last and next
            let mut last;
            let mut mid;
            let mut next;

            let mut step = bound;
            // Calculates up to second-greatest depth
            while step > 1 {
                step /= 2;
                index = step;
                last = 0;
                while index <= bound/2 {
                    next = Self::partition_sum(&partition, index);
                    if step == 1 {
                        if partition[index-1] {
                            mid = (last + next) / 2;
                        } else {
                            mid = last;
                        }
                    } else {
                        mid = Self::partition_sum(&partition, index - step / 2);
                    }

                    for symbol in sequence.iter() {
                        if symbol >= &alphabet[last] && symbol <= &alphabet[next-1] {
                            bits.push((mid != last) && (symbol >= &alphabet[mid]));
                        }
                    }

                    last = next;
                    index += step;
                }
            }

            return RankSelect::new(bits, 1);
        }
    }

    pub fn access(&self, i: u64) -> Option<&T> {
        // Upper bound and its log to base 2
        let bound = Self::alphabet_bound(self.alphabet.len());
        let log = Self::bound_log2(bound);

        // Defines how the Alphabet is partitioned
        let partition = Self::partition_alphabet(bound, self.alphabet.len());

        // Returns if Alphabet is empty == Empty String
        if self.alphabet.len() == 0 {return Option::None;}
        // Returns Result for one Symbol in Alphabet
        if self.alphabet.len() == 1 {
            if i < self.bitmap.bits().len() {
                return Option::Some(&self.alphabet[0]);
            } else {
                return Option::None;
            }
        }
        // Returns if Index out of Bounds
        if i >= self.bitmap.bits().len() / (log as u64) {
            return Option::None
        }
        
        // Calculates index, start and end of each new Layer till second to last Layer
        let mut index = i as u64;
        let level_len = self.bitmap.bits().len() / (log as u64);
        let mut depth_start = 0;
        let mut start = 0;
        let mut end = level_len - 1;
        let mut start_index = 0;
        let mut end_index = partition.len();
        while (end_index - start_index) > 1 {
            if self.bitmap.get(depth_start + start + index) {
                start_index = (start_index + end_index) / 2;
                if index > 0 {index = self.bitmap.rank_1(depth_start + start + index).unwrap() - self.bitmap.rank_1(depth_start + start).unwrap() - if !self.bitmap.get(depth_start + start) {1} else {0};}
                start = start + self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + if !self.bitmap.get(depth_start + start) {1} else {0};
            } else {
                end_index = (start_index + end_index) / 2;
                if index > 0 {index = self.bitmap.rank_0(depth_start + start + index).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() - if self.bitmap.get(depth_start + start) {1} else {0};}
                end = start + self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() - if self.bitmap.get(depth_start + start) {1} else {0};
            }
            depth_start += level_len;
        }

        // Returns Result for last Layer
        if partition[start_index] {
            if self.bitmap.get(depth_start + start + index) {
                return Option::Some(&self.alphabet[Self::partition_sum(&partition, start_index)+1]);
            } else {
                return Option::Some(&self.alphabet[Self::partition_sum(&partition, start_index)]);
            }
        } else {
            return Option::Some(&self.alphabet[Self::partition_sum(&partition, start_index)]);
        }
    }

    pub fn rank(&self, symbol: &T, i: u64) -> u64 {
        // Upper bound and its log to base 2
        let bound = Self::alphabet_bound(self.alphabet.len());
        let log = Self::bound_log2(bound);

        // Defines how the Alphabet is partitioned
        let partition = Self::partition_alphabet(bound, self.alphabet.len());
        let symbol_in_alphabet = self.alphabet.iter().position(|x| x == symbol).unwrap_or_default();
        let symbol_index = Self::partition_slice(&partition, symbol_in_alphabet);

        // Returns if Alphabet is empty == Empty String or Alphabet does not contain Symbol
        if (self.alphabet.len() == 0) || !self.alphabet.contains(symbol) {return 0;}
        // Returns Result for one Symbol in Alphabet
        if (self.alphabet.len() == 1) && self.alphabet.contains(symbol) {
            if i < self.bitmap.bits().len() {
                return i+1;
            } else {
                return self.bitmap.bits().len();
            }
        }

        // Calculates index, start and end of each new Layer till second to last Layer, while returning 0 if index smaller than first occurence of symbol
        let mut index = i as u64;
        let level_len = self.bitmap.bits().len() / (log as u64);
        let mut depth_start = 0;
        let mut start = 0;
        let mut end = level_len - 1;
        let mut start_index = 0;
        let mut end_index = partition.len();
        if index >= end {index = end;}
        while (end_index - start_index) > 1 {
            if symbol_index >= ((start_index + end_index) / 2) {
                start_index = (start_index + end_index) / 2;
                if index == 0 && !self.bitmap.get(depth_start + start) {return 0;}
                if index > 0 {index = self.bitmap.rank_1(depth_start + start + index).unwrap() - self.bitmap.rank_1(depth_start + start).unwrap() - if !self.bitmap.get(depth_start + start) {1} else {0};}
                start = start + self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + if !self.bitmap.get(depth_start + start) {1} else {0};
            } else {
                end_index = (start_index + end_index) / 2;
                if index == 0 && self.bitmap.get(depth_start + start) {return 0;}
                if index > 0 {index = self.bitmap.rank_0(depth_start + start + index).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() - if self.bitmap.get(depth_start + start) {1} else {0};}
                end = start + self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() - if self.bitmap.get(depth_start + start) {1} else {0};
            }
            depth_start += level_len;
        }

        // Returns Result for last Layer
        if partition[start_index] {
            if symbol_in_alphabet >= Self::partition_sum(&partition, start_index)+1 {
                return self.bitmap.rank_1(depth_start + start + index).unwrap() - self.bitmap.rank_1(depth_start + start).unwrap() + if self.bitmap.get(depth_start + start) {1} else {0};
            } else {
                return self.bitmap.rank_0(depth_start + start + index).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + if !self.bitmap.get(depth_start + start) {1} else {0};
            }
        } else {
            return self.bitmap.rank_0(depth_start + start + index).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + 1;
        }
    }

    /// Operation SELECT: returns the position of the i-th occurence of symbol c in the sequence represented by this 
    /// wavelet tree
    /// Note that indices start at 0. Therefore, the first occurence (i=1) of the first symbol in the sequence would be
    /// returned as Some(0)!
    pub fn select(&self, c: T, i: u64) -> Option<u64> {
        // If the alphabet contains c, we can execute a SELECT operation.
        // Otherwise, None is returned.
        if self.alphabet.contains(&c) {
            // Calculate n
            let n = self.bitmap.bits().len() / (Self::bound_log2(Self::alphabet_bound(self.alphabet.len())) as u64);
            // Compute select(c,i) on the root node (bitmap subrange 0..n, representing the whole alphabet)
            self.select_on_node(c, i, 0, n, 0, self.alphabet.len() as u64, n)
        } else {
            None
        }

    }

    /// Operation SELECT on the node corresponding to the bitmap given by the range [l,r). The node represents the
    /// subrange [a,b) of the tree's alphabet.
    fn select_on_node(&self, c: T, i: u64, l: u64, r: u64, a: u64, b: u64, n: u64) -> Option<u64> {

        // If we are not on leaf level yet (i.e. if the alphabet represented by this node contains more than 1 symbol),
        // we need to move "downwards" in order to find the leaf corresponding to c.
        if b-a > 1 {
            // position where the alphabet represented by the current node is cut by its children
            let alphabet_cut = a + 2u64.pow((((b-a) as f64).log2().ceil() as u32) -1); 

            // Moving down works similar to ACCESS: we look at the subtree where c can be found, which depends on the
            // position of c in relation to the center of the alphabet represented by the current node.
            if c < self.alphabet[alphabet_cut as usize] {
                println!(" => links.");
                // Compute interval of bitmap and alphabet for the right child
                let new_l = n+l;
                let new_r = n+l + self.bitmap.rank_0(r).unwrap() - self.bitmap.rank_0(l).unwrap();
                let new_a = a;
                let new_b = alphabet_cut;

                // Recursively compute position on lower level
                let p = self.select_on_node(c, i, new_l, new_r, new_a, new_b, n);

                // The variable p indicates the position of the i-th occurence of symbol c on the next level. The 
                // position of the i-th occurence of c on the current level now corresponds to the (p+1)-th occurence of
                // a 0 on the current level, i.e. select_0(p+1) on the current level.
                // However, the SELECT operation needs to be carried out on the bitmap representing the whole tree. 
                // Therefore, we need to add the number of appearances of 0 before the beginning of the [l,r)-part of
                // the bitmap, which is rank_0(l-1), to (p+1); the position returned by SELECT needs to be adjusted as
                // well by subtracting the starting position of the current node's bitmap, i.e. l.

                let offset = if l > 0 {self.bitmap.rank_0(l-1).unwrap()} else {0};

                p.map(|pos| pos + 1)
                 .map(|pos| pos + offset)
                 .and_then(|pos| self.bitmap.select_0(pos))
                 .map(|selected| selected - l)

                // p.map(|pos| self.bitmap.select_0(pos+1).unwrap())
            } else {
                println!(" => rechts.");
                // Compute interval of bitmap and alphabet for the right child
                let new_l = n+l + self.bitmap.rank_0(r).unwrap() - self.bitmap.rank_0(l).unwrap();
                let new_r = n+r;
                let new_a = alphabet_cut;
                let new_b = b;

                // Recursively compute position on lower level
                let p = self.select_on_node(c, i, new_l, new_r, new_a, new_b, n);

                // The variable p indicates the position of the i-th occurence of symbol c on the next level. The 
                // position of the i-th occurence of c on the current level now corresponds to  select_1(p+1) on the 
                // current level. For details see comments above.
                let offset = if l > 0 {self.bitmap.rank_1(l-1).unwrap()} else {0};

                p.map(|pos| pos + 1)
                 .map(|pos| pos + offset)
                 .and_then(|pos| self.bitmap.select_1(pos))
                 .map(|selected| selected - l)
            }
        } else {
            // At leaf level, the i-th occurence of c is at position (i-1).
            Some(i-1)
        }
    }

    // Calculates total log of bound (upper boundary)
    pub fn bound_log2(bound: usize) -> usize {
        let mut log = 0;
        let mut new_bound = 1;
        while new_bound < bound {
            log += 1;
            new_bound *= 2;
        }
        log
    }

    // Calculates smallest total number d with alphabet.length <= 2^d
    pub fn alphabet_bound(alphabetlen: usize) -> usize {
        let mut bound = 1;
        while bound < alphabetlen {
            bound *= 2;
        };
        bound
    }

    // Calculates how the Alphabet is partitioned in Tree, true -> symbol in greatest depth, false -> symbol in second-greatest depth
    pub fn partition_alphabet(bound: usize, alphabetlen: usize) -> Vec<bool> {
        let mut part: Vec<bool> = Vec::new();

        for i in 0..bound/2 {
            part.push(i < (bound/2 + alphabetlen - bound));
        }

        part
    }

    // Calculates Symbol Count in Partition up to index
    pub fn partition_sum(partition: &Vec<bool>, index: usize) -> usize {
        let mut sum = 0;
        for i in 0..index {
            if partition[i] {
                sum += 2;
            }
            else {
                sum += 1;
            }
        }
        sum
    }

    // Calculates index in Partition from index in Alphabet
    pub fn partition_slice(partition: &Vec<bool>, index: usize) -> usize {
        let mut slice = 0;
        for i in 0..partition.len()-1 {
            if partition[i] {
                slice += 2;
            } else {
                slice += 1;
            }
            if slice > index {
                return i;
            } else if slice == index {
                return i+1;
            }
        }
        partition.len()-1
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

    // TESTS FOR POINTERLESS WAVELET TREE START HERE

    #[test]
    fn test_pointerless_wavelet_tree_from_sequence() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_sequence(sequence);

        // The correct alphabet should automatically be created
        assert_eq!(tree.alphabet, vec![' ','a','b','d','l','r']);

        // The bitmap of the pointerless wavelet tree should have a specific format
        // (trailing 0s are needed to compute n)
        // The bitmap is implicitly divided as follows:
        // level 1: 01000100010001000100
        // level 2: 001000000001010 | 01001
        // level 3: 111010101111 | 001
        let expected_bit_string = "010001000100010001000010000000010100100111101010111100100000";

        assert_eq!(expected_bit_string, bit_vec_to_string(tree.bitmap.bits()));
    }

    fn bit_vec_to_string(bit_vector: &BitVec<u8>) -> String {
        let n = bit_vector.len();

        let mut bit_string = String::with_capacity(bit_vector.len() as usize);

        for i in 0..n {
            if bit_vector[i] {
                bit_string.push('1');
            } else {
                bit_string.push('0');
            }
        }

        bit_string
    }

    #[test]
    fn test_pointerless_access_inside_range() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_sequence(sequence);

        assert_eq!(Some('a').as_ref(), tree.access(0));
        assert_eq!(Some('l').as_ref(), tree.access(1));
        assert_eq!(Some('a').as_ref(), tree.access(2));
        assert_eq!(Some('b').as_ref(), tree.access(3));
        assert_eq!(Some('a').as_ref(), tree.access(4));
        assert_eq!(Some('r').as_ref(), tree.access(5));
        assert_eq!(Some(' ').as_ref(), tree.access(6));
        assert_eq!(Some('a').as_ref(), tree.access(7));
        assert_eq!(Some(' ').as_ref(), tree.access(8));
        assert_eq!(Some('l').as_ref(), tree.access(9));
        assert_eq!(Some('a').as_ref(), tree.access(10));
        assert_eq!(Some(' ').as_ref(), tree.access(11));
        assert_eq!(Some('a').as_ref(), tree.access(12));
        assert_eq!(Some('l').as_ref(), tree.access(13));
        assert_eq!(Some('a').as_ref(), tree.access(14));
        assert_eq!(Some('b').as_ref(), tree.access(15));
        assert_eq!(Some('a').as_ref(), tree.access(16));
        assert_eq!(Some('r').as_ref(), tree.access(17));
        assert_eq!(Some('d').as_ref(), tree.access(18));
        assert_eq!(Some('a').as_ref(), tree.access(19));
    }

    #[test]
    fn test_pointerless_access_out_of_range() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_sequence(sequence);

        assert_eq!(None, tree.access(20));
    }

    #[test]
    fn test_pointerless_select_first_appearance_of_each_symbol() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_sequence(sequence);

        assert_eq!(Some(0), tree.select('a', 1));
        assert_eq!(Some(1), tree.select('l', 1));
        assert_eq!(Some(3), tree.select('b', 1));
        assert_eq!(Some(5), tree.select('r', 1));
        assert_eq!(Some(6), tree.select(' ', 1));
        assert_eq!(Some(18), tree.select('d', 1));
}

    #[test]
    fn test_pointerless_select_all_appearances_of_one_symbol() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_sequence(sequence);

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
    fn test_pointerless_select_invalid_argument() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_sequence(sequence);

        // Out of range
        assert_eq!(None, tree.select('a', 20));
        // Symbol not in alphabet
        assert_eq!(None, tree.select('x', 1));
        // i too high
        assert_eq!(None, tree.select('d', 2));
    }
}
