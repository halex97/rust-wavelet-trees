use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;
use rand::Rng;
use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize)]
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
    // Create new Pointerless Wavelet Tree
    pub fn create_bits(sequence: &[T], alphabet: &[T]) -> RankSelect {
        if sequence.is_empty() {
            // For empty Sequence the bitmap is created empty
            return RankSelect::new(BitVec::new(), 1);
        }
        else {
            let mut bits: BitVec<u8> = BitVec::new();

            // Special Case if Alphabet contains only one Element
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
                    // Calculates mid position depending on Position in Tree
                    if step == 1 {
                        if partition[index-1] {
                            mid = (last + next) / 2;
                        } else {
                            mid = last;
                        }
                    } else {
                        mid = Self::partition_sum(&partition, index - step / 2);
                    }

                    // Calculates and adds Node bitmap
                    for symbol in sequence.iter() {
                        if symbol >= &alphabet[last] && symbol <= &alphabet[next-1] {
                            bits.push((mid != last) && (symbol >= &alphabet[mid]));
                        }
                    }

                    last = next;
                    index += step;
                }
            }

            // Returns calculated Tree bitmap
            return RankSelect::new(bits, 1);
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

impl <T: PartialOrd + Clone> super::WaveletTree<T> for PointerlessWaveletTree<T> {
    
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

        // Create bitmap based on Sequence and Alphabet
        let bitmap = Self::create_bits(sequence, &alphabet);

        // Return a new PointerWaveletTree containing the alphabet and the root of the tree
        PointerlessWaveletTree {
            alphabet,
            bitmap
        }
    }

    fn rank(&self, symbol: &T, i: u64) -> Option<u64> {
        // Upper bound and its log to base 2
        let bound = Self::alphabet_bound(self.alphabet.len());
        let log = Self::bound_log2(bound);

        // Returns if Alphabet is empty or Alphabet does not contain Symbol
        if (self.alphabet.len() == 0) || !self.alphabet.contains(symbol) {return Option::None;}
        // Returns Result for one Symbol in Alphabet
        if (self.alphabet.len() == 1) && self.alphabet.contains(symbol) {
            if i < self.bitmap.bits().len() {
                return Option::Some(i+1);
            } else {
                return Option::None;
            }
        }
        // Returns for index out of bounds
        if i >= self.bitmap.bits().len() / (log as u64) {return Option::None;}

        // Defines how the Alphabet is partitioned
        let partition = Self::partition_alphabet(bound, self.alphabet.len());
        let symbol_in_alphabet = self.alphabet.iter().position(|x| x == symbol).unwrap_or_default();
        let symbol_index = Self::partition_slice(&partition, symbol_in_alphabet);

        // Calculates index, start and end of each new Layer till second to last Layer, while returning 0 if index smaller than first occurence of symbol
        let mut index = i as u64;
        let level_len = self.bitmap.bits().len() / (log as u64);
        let mut depth_start = 0;
        let mut start = 0;
        let mut mid;
        let mut end = level_len - 1;
        let mut start_index = 0;
        let mut mid_index;
        let mut end_index = partition.len();
        while (end_index - start_index) > 1 {
            mid = start + self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + if !self.bitmap.get(depth_start + start) {1} else {0};
            mid_index = (start_index + end_index) / 2;
            if symbol_index >= mid_index {
                if self.bitmap.rank_1(depth_start + start + index).unwrap() == self.bitmap.rank_1(depth_start + start).unwrap() && !self.bitmap.get(depth_start + start) {return Option::Some(0);}
                index = self.bitmap.rank_1(depth_start + start + index).unwrap() - self.bitmap.rank_1(depth_start + start).unwrap() - if !self.bitmap.get(depth_start + start) {1} else {0};
                start = mid;
                start_index = mid_index;
            } else {
                if self.bitmap.rank_0(depth_start + start + index).unwrap() == self.bitmap.rank_0(depth_start + start).unwrap() && self.bitmap.get(depth_start + start) {return Option::Some(0);}
                index = self.bitmap.rank_0(depth_start + start + index).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() - if self.bitmap.get(depth_start + start) {1} else {0};
                end = mid - 1;
                end_index = mid_index;
            }
            depth_start += level_len;
        }

        // Returns Result for last Layer
        if partition[start_index] && (symbol_in_alphabet >= Self::partition_sum(&partition, start_index)+1) {
            return Option::Some(self.bitmap.rank_1(depth_start + start + index).unwrap() - self.bitmap.rank_1(depth_start + start).unwrap() + if self.bitmap.get(depth_start + start) {1} else {0});
        } else {
            return Option::Some(self.bitmap.rank_0(depth_start + start + index).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + if !self.bitmap.get(depth_start + start) {1} else {0});
        }
    }

    fn select(&self, symbol: &T, i: u64) -> Option<u64> {
        // Upper bound and its log to base 2
        let bound = Self::alphabet_bound(self.alphabet.len());
        let log = Self::bound_log2(bound);
        
        // Returns if Alphabet is empty or Alphabet does not contain Symbol or index too small
        if (self.alphabet.len() == 0) || !self.alphabet.contains(symbol) || i == 0 {return Option::None;}
        // Returns Result for one Symbol in Alphabet
        if (self.alphabet.len() == 1) && self.alphabet.contains(symbol) {
            if i <= self.bitmap.bits().len() {
                return Option::Some(i-1);
            } else {
                return Option::None;
            }
        }

        // Defines how the Alphabet is partitioned
        let partition = Self::partition_alphabet(bound, self.alphabet.len());
        let symbol_in_alphabet = self.alphabet.iter().position(|x| x == symbol).unwrap_or_default();
        let symbol_index = Self::partition_slice(&partition, symbol_in_alphabet);

        // Vectors for traversing tree upwards after traversing downwards
        let mut start_points: Vec<u64> = Vec::new();
        let mut index_points: Vec<bool> = Vec::new();

        // Calculates index, start and end of each new Layer till second to last Layer, while saving data necessary for upwards traversal of Tree
        let mut index = i as u64 - 1;
        let level_len = self.bitmap.bits().len() / (log as u64);
        let mut depth_start = 0;
        let mut start = 0;
        let mut mid;
        let mut end = level_len - 1;
        let mut start_index = 0;
        let mut mid_index;
        let mut end_index = partition.len();
        while (end_index - start_index) > 1 {
            start_points.push(start);
            mid = start + self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + if !self.bitmap.get(depth_start + start) {1} else {0};
            mid_index = (start_index + end_index) / 2;
            if symbol_index >= ((start_index + end_index) / 2) {
                index_points.push(true);
                start = mid;
                start_index = mid_index;
            } else {
                index_points.push(false);
                end = mid - 1;
                end_index = mid_index;
            }
            depth_start += level_len;
        }

        // Start Of traversal Upwards in Tree
        if partition[start_index] {
            if symbol_in_alphabet >= Self::partition_sum(&partition, start_index)+1 {
                if index >= self.bitmap.rank_1(depth_start + end).unwrap() - self.bitmap.rank_1(depth_start + start).unwrap() + if self.bitmap.get(depth_start + start) {1} else {0} {return Option::None;}
                else {index = self.bitmap.select_1(index + self.bitmap.rank_1(depth_start + start).unwrap() - if self.bitmap.get(depth_start + start) {1} else {0} + 1).unwrap() - depth_start - start;}
            } else {
                if index >= self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + if !self.bitmap.get(depth_start + start) {1} else {0} {return Option::None;}
                else {index = self.bitmap.select_0(index + self.bitmap.rank_0(depth_start + start).unwrap() - if !self.bitmap.get(depth_start + start) {1} else {0} + 1).unwrap() - depth_start - start;}
            }
        } else {
            if index >= self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + 1 {return Option::None;}
        }

        // Full Upward traversal, unraveling Data while calculating indizes
        while start_points.len() > 0 {
            depth_start -= level_len;
            start = start_points.pop().unwrap();
            index =
            if index_points.pop().unwrap() {self.bitmap.select_1(index + self.bitmap.rank_1(depth_start + start).unwrap() + if self.bitmap.get(depth_start + start) {0} else {1}).unwrap()}
            else {self.bitmap.select_0(index + self.bitmap.rank_0(depth_start + start).unwrap() + if !self.bitmap.get(depth_start + start) {0} else {1}).unwrap()} - depth_start - start;
        }

        return Option::Some(index);
    }

    fn access(&self, i: u64) -> Option<&T> {
        // Upper bound and its log to base 2
        let bound = Self::alphabet_bound(self.alphabet.len());
        let log = Self::bound_log2(bound);

        // Defines how the Alphabet is partitioned
        let partition = Self::partition_alphabet(bound, self.alphabet.len());

        // Returns if Alphabet is empty
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
        let mut mid;
        let mut end = level_len - 1;
        let mut start_index = 0;
        let mut mid_index;
        let mut end_index = partition.len();
        while (end_index - start_index) > 1 {
            mid = start + self.bitmap.rank_0(depth_start + end).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() + if !self.bitmap.get(depth_start + start) {1} else {0};
            mid_index = (start_index + end_index) / 2;
            if self.bitmap.get(depth_start + start + index) {
                if index > 0 {index = self.bitmap.rank_1(depth_start + start + index).unwrap() - self.bitmap.rank_1(depth_start + start).unwrap() - if !self.bitmap.get(depth_start + start) {1} else {0};}
                start = mid;
                start_index = mid_index;
            } else {
                if index > 0 {index = self.bitmap.rank_0(depth_start + start + index).unwrap() - self.bitmap.rank_0(depth_start + start).unwrap() - if self.bitmap.get(depth_start + start) {1} else {0};}
                end = mid - 1;
                end_index = mid_index;
            }
            depth_start += level_len;
        }

        // Returns Result for last Layer
        if partition[start_index] && self.bitmap.get(depth_start + start + index) {
            return Option::Some(&self.alphabet[Self::partition_sum(&partition, start_index)+1]);
        } else {
            return Option::Some(&self.alphabet[Self::partition_sum(&partition, start_index)]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::*;

    #[test]
    fn test_pointerless_wavelet_tree_from_sequence() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_slice(sequence);

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
        let tree = PointerlessWaveletTree::from_slice(sequence);

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
        let tree = PointerlessWaveletTree::from_slice(sequence);

        assert_eq!(None, tree.access(20));
    }

    #[test]
    fn test_pointerless_select_first_appearance_of_each_symbol() {
        let text = "alabar a la alabarda";
        let sequence : &Vec<char> = &text.chars().collect();
        let tree = PointerlessWaveletTree::from_slice(sequence);

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
        let tree = PointerlessWaveletTree::from_slice(sequence);

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
        let tree = PointerlessWaveletTree::from_slice(sequence);

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
        let tree = PointerlessWaveletTree::from_slice(sequence);

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
        let tree = PointerlessWaveletTree::from_slice(sequence);

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
        let tree = PointerlessWaveletTree::from_slice(sequence);

        let mut pos = 0;
        for symbol in tree.alphabet.clone().into_iter().rev() {
            assert_eq!(Some(pos), tree.select(&symbol, 1));
            pos = pos + 1;
        }
    }

    #[test]
    #[ignore]
    fn test_pointerless_randomized_access_rank_select() {
        // Time Nedded ~ 10-15 min
        PointerlessWaveletTree::from_slice(&Vec::<u64>::new());
        let mut sequence: Vec<u64> = Vec::new();
        sequence.push(1);
        let tree = PointerlessWaveletTree::from_slice(&sequence);
        assert_eq!(Option::Some(&1), tree.access(0));
        assert_eq!(Option::None, tree.access(1));
        assert_eq!(Option::Some(1), tree.rank(&1, 0));
        assert_eq!(Option::None, tree.rank(&1, 1));
        assert_eq!(Option::Some(0), tree.select(&1, 1));
        assert_eq!(Option::None, tree.rank(&1, 2));

        let mut numbergen = rand::thread_rng();
        for size in 1..256 {
            // Build Alphabet
            let mut alphabet: Vec<u64> = Vec::new();
            for i in 0..size {alphabet.push(i);}

            // Build Number Vector
            let mut numbers: Vec<u64> = alphabet.to_vec();
            for _j in 0..(32+size/4-(size*size)/1024) {numbers.push(numbers[numbergen.gen_range(0, size) as usize]);}
            let tree = PointerlessWaveletTree::from_slice(&numbers);

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
                    assert_eq!(index as u64, tree.select(&symbol, j as u64).unwrap());
                    index += 1;
                };
                // Test Select Invalid
                for j in count+1..count+256 {assert_eq!(Option::None, tree.select(&symbol, j as u64))}
            }
        }
    }
}