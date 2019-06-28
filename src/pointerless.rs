use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;

pub struct PointerlessWaveletTree<T: PartialOrd + Clone> {
    alphabet: Vec<T>,
    bitmap: RankSelect
}

impl <T: PartialOrd + Clone> std::fmt::Debug for PointerlessWaveletTree<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut bitmap_string = String::new();

        for i in 0..self.bitmap.bits().len() {
            if self.bitmap.get(i) {
                bitmap_string.push_str("1");
            } else {
                bitmap_string.push_str("0");
            }
        }

        write!(f, "PointerlessWaveletTree {{ bitmap: {} }}", bitmap_string)
    }
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

            // Special Case if Alphabet contains one Element
            if alphabet.len() == 1 {
                bits = BitVec::new_fill(false, sequence.len() as u64);
                return RankSelect::new(bits, 1);
            }

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

    pub fn select(&self, symbol: &T, i: u64) -> Option<u64> {
        // Upper bound and its log to base 2
        let bound = Self::alphabet_bound(self.alphabet.len());
        let log = Self::bound_log2(bound);

        // Defines how the Alphabet is partitioned
        let partition = Self::partition_alphabet(bound, self.alphabet.len());
        let symbol_in_alphabet = self.alphabet.iter().position(|x| x == symbol).unwrap_or_default();
        let symbol_index = Self::partition_slice(&partition, symbol_in_alphabet);

        // Returns if Alphabet is empty == Empty String or Alphabet does not contain Symbol or index too small
        if (self.alphabet.len() == 0) || !self.alphabet.contains(symbol) || i == 0 {return Option::None;}
        // Returns Result for one Symbol in Alphabet
        if (self.alphabet.len() == 1) && self.alphabet.contains(symbol) {
            if i < self.bitmap.bits().len() {
                return Option::Some(i-1);
            } else {
                return Option::None;
            }
        }

        // Vectors for traversing tree upwards after traversing downwards
        let mut start_points: Vec<u64> = Vec::new();
        let mut index_points: Vec<bool> = Vec::new();

        // Calculates index, start and end of each new Layer till second to last Layer, while saving data necessary for upwards traversal of Tree
        let mut index = i as u64 - 1;
        let level_len = self.bitmap.bits().len() / (log as u64);
        let mut depth_start = 0;
        let mut start = 0;
        let mut end = level_len - 1;
        let mut start_index = 0;
        let mut end_index = partition.len();
        while (end_index - start_index) > 1 {
            start_points.push(start);
            if symbol_index >= ((start_index + end_index) / 2) {
                start_index = (start_index + end_index) / 2;
                start = start + self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + if !self.bitmap.get(depth_start + start) {1} else {0};
                index_points.push(true);
            } else {
                end_index = (start_index + end_index) / 2;
                end = start + self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() - if self.bitmap.get(depth_start + start) {1} else {0};
                index_points.push(false);
            }

            depth_start += level_len;
        }

        // Variable for fixing bug in dependant library (RankSelect)
        let interim_result;

        // Start Of traversal Upwards in Tree
        if partition[start_index] {
            if symbol_in_alphabet >= Self::partition_sum(&partition, start_index)+1 {
                interim_result = index + self.bitmap.rank_1(depth_start + start).unwrap() - if self.bitmap.get(depth_start + start) {1} else {0} + 1;
                if index >= self.bitmap.rank_1(depth_start + end).unwrap() - self.bitmap.rank_1(depth_start + start).unwrap() + if self.bitmap.get(depth_start + start) {1} else {0} {return Option::None;}
                else {index = self.bitmap.select_1(interim_result).unwrap() - depth_start - start;}
            } else {
                interim_result = index + self.bitmap.rank_0(depth_start + start).unwrap() - if !self.bitmap.get(depth_start + start) {1} else {0} + 1;
                if index >= self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + if !self.bitmap.get(depth_start + start) {1} else {0} {return Option::None;}
                else {index = self.bitmap.select_0(interim_result).unwrap() - depth_start - start;}
            }
        } else {
            if index >= self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + 1 {return Option::None;}
        }

        // Full Upward traversal, unraveling Data while calculating indizes
        while start_points.len() > 0 {
            depth_start -= level_len;
            start = start_points.pop().unwrap();
            index = if index_points.pop().unwrap() {self.bitmap.select_1(index + self.bitmap.rank_1(depth_start + start).unwrap() - if self.bitmap.get(depth_start + start) {1} else {0} + 1).unwrap()}
            else {self.bitmap.select_0(index + self.bitmap.rank_0(depth_start + start).unwrap() - if !self.bitmap.get(depth_start + start) {1} else {0} + 1).unwrap()} - depth_start - start;
        }

        return Option::Some(index);
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

impl <T: PartialOrd + Clone> super::WaveletTree<T> for PointerlessWaveletTree<T> {
    fn from_slice(sequence: &[T]) -> Self {
        unimplemented!();
    }

    fn rank(&self, c: &T, i: u64) -> Option<u64> {
        unimplemented!();
    }
    fn select(&self, c: &T, i: u64) -> Option<u64> {
        unimplemented!();
    }
    fn access(&self, i: u64) -> Option<&T> {
        unimplemented!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

        assert_eq!(Some(0), tree.select(&'a', 1));
        assert_eq!(Some(1), tree.select(&'l', 1));
        assert_eq!(Some(3), tree.select(&'b', 1));
        assert_eq!(Some(5), tree.select(&'r', 1));
        assert_eq!(Some(6), tree.select(&' ', 1));
        assert_eq!(Some(18), tree.select(&'d', 1));
    }

    #[test]
    fn test_pointerless_select_all_appearances_of_one_symbol() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_sequence(sequence);

        assert_eq!(Some(0), tree.select(&'a', 1));
        assert_eq!(Some(2), tree.select(&'a', 2));
        assert_eq!(Some(4), tree.select(&'a', 3));
        assert_eq!(Some(7), tree.select(&'a', 4));
        assert_eq!(Some(10), tree.select(&'a', 5));
        assert_eq!(Some(12), tree.select(&'a', 6));
        assert_eq!(Some(14), tree.select(&'a', 7));
        assert_eq!(Some(16), tree.select(&'a', 8));
        assert_eq!(Some(19), tree.select(&'a', 9));
    }

    #[test]
    fn test_pointerless_select_invalid_argument() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_sequence(sequence);

        // Out of range
        assert_eq!(None, tree.select(&'a', 20));
        // Symbol not in alphabet
        assert_eq!(None, tree.select(&'x', 1));
        // i too high
        assert_eq!(None, tree.select(&'d', 2));
        assert_eq!(None, tree.select(&'l', 4));
        assert_eq!(None, tree.select(&'r', 3));
    }

    #[test]
    fn test_pointerless_select_all_valid_arguments() {
        let text = "The quick brown fox jumps over a lazy dog";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_sequence(sequence);

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
                assert_eq!(Option::Some(index as u64), tree.select(&symbol, j as u64));
                index += 1;
            };
        }
    }

    #[test]
    fn test_pointerless_select_ascending_alphabet_as_sequence() {
        let text = "abcdefghijklmnopqrstuvwxyz";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_sequence(sequence);

        let mut pos = 0;
        for symbol in tree.alphabet.clone().into_iter() {
            assert_eq!(Some(pos), tree.select(&symbol, 1));
            pos = pos + 1;
        }
    }

    #[test]
    fn test_pointerless_select_descending_alphabet_as_sequence() {
        let text = "zyxwvutsrqponmlkjihgfedcba";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_sequence(sequence);

        let mut pos = 0;
        for symbol in tree.alphabet.clone().into_iter().rev() {
            assert_eq!(Some(pos), tree.select(&symbol, 1));
            pos = pos + 1;
        }
    }
}