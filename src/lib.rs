use bio::data_structures::rank_select::RankSelect;

pub struct PointerWaveletTree<T: PartialOrd + Clone> {
    bitmap: Option<RankSelect>,
    label: Option<T>,
    left_child: Option<Box<PointerWaveletTree<T>>>,
    right_child: Option<Box<PointerWaveletTree<T>>>
}

impl <T: PartialOrd + Clone> PointerWaveletTree<T> {

    /// Create a new wavelet tree
    pub fn new(sequence: Vec<T>) -> PointerWaveletTree<T> {
        if sequence.is_empty() {
            panic!("Wavelet trees cannot be created from empty sequences.");
        }

        let first_symbol = &sequence[0];

        // If there is only one symbol in the alphabet of the sequence (i.e. all symbols in the sequence are equal),
        // a leaf should be created.
        // Otherwise we create a bitmap for the sequence and use it to create children.

        if sequence.len() == 1 || sequence.iter().fold(true, |acc, x| acc && (x == first_symbol)) {
            // Return a leaf (= a wavelet tree with a label and without children)
            PointerWaveletTree {
                bitmap: None,
                label: Some(sequence[0].clone()),
                left_child: None,
                right_child: None
            }
        } else {
            // Create the bitmap.
            let bitmap = create_bitmap(&sequence);

            // Now we can split up the sequence of this wavelet tree in order to create the left and the right child.
            let mut left_sequence : Vec<T> = Vec::new();
            let mut right_sequence : Vec<T> = Vec::new();

            for i in 0..sequence.len() {
                // If the bitmap contains a 1 at position i, the symbol of the sequence at position i should be
                // in the right child's sequence.
                // Otherwise it should be a part of the left child's sequence.
                if bitmap.get(i as u64) {
                    left_sequence.push(sequence[i].clone());
                } else {
                    right_sequence.push(sequence[i].clone());
                }
            }

            // We can now recursively create the children of this wavelet tree.
            // A child only needs to be created if it's corresponding sequence is not empty.
            let left_child;
            let right_child;

            if left_sequence.is_empty() {
                left_child = None;
            } else {
                left_child = Some(Box::new(PointerWaveletTree::new(left_sequence)));
            }

            if right_sequence.is_empty() {
                right_child = None;
            } else {
                right_child = Some(Box::new(PointerWaveletTree::new(right_sequence)));
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

}

fn create_bitmap<T: PartialOrd>(sequence: &Vec<T>) -> RankSelect {
    unimplemented!()
}


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
