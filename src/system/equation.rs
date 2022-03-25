use crate::symbols::arithmetic;
use crate::symbols::Symbol;
use crate::symbols::SymbolID;
use crate::system::equation::EquationBounds::{Bounds, EmptyBounds};

#[derive(Debug, PartialEq, Eq)]
pub enum EquationBounds {
    Bounds {
        pivot: SymbolID,
        last_nonzero_id: SymbolID,
    },
    EmptyBounds,
}

#[derive(Debug)]
pub struct Equation {
    coefs: Vec<u8>,
    constant_term: Symbol,
    bounds: EquationBounds,
}

impl Equation {
    pub fn new(mut coefs: Vec<u8>, constant_term: Symbol) -> Equation {
        if coefs.len() != constant_term.n_protected_symbols as usize {
            coefs.resize(constant_term.n_protected_symbols as usize, 0u8);
        }
        let mut eq = Equation {
            coefs,
            constant_term,
            bounds: EmptyBounds,
        };
        eq.recompute_bounds();
        eq
    }

    pub fn _new(mut coefs: Vec<u8>, constant_term: Symbol, bounds: EquationBounds) -> Equation {
        if coefs.len() != constant_term.n_protected_symbols as usize {
            coefs.resize(constant_term.n_protected_symbols as usize, 0u8);
        }
        let mut eq = Equation {
            coefs,
            constant_term,
            bounds
        };
        eq
    }

    pub fn constant_term_data(self) -> Vec<u8> {
        return self.constant_term.take_data();
    }

    pub fn add(&mut self, other: &Equation) -> Result<(), ()> {
        assert_eq!(
            self.constant_term.n_protected_symbols as usize,
            self.coefs.len()
        );
        if self.is_zero() && other.is_zero() {
            return Ok(());
        }

        if other.is_zero() {
            // nothing changes
            return Ok(());
        }

        // TODO: equation pivot is not necessarily at index 0 of the coefs ! index 0 is for the repair symbols' first_id field !

        if let (
            Bounds {
                pivot: self_pivot,
                last_nonzero_id: self_last_nonzero_id,
            },
            Bounds {
                pivot: other_pivot,
                last_nonzero_id: other_last_nonzero_id,
            },
        ) = (self.bounds(), other.bounds())
        {
            let (self_pivot, self_last_nonzero_id, other_pivot, other_last_nonzero_id) = (
                *self_pivot,
                *self_last_nonzero_id,
                *other_pivot,
                *other_last_nonzero_id,
            );
            let smallest_pivot = std::cmp::min(self_pivot, other_pivot);
            let largest_nonzero_coef_id =
                std::cmp::max(self_last_nonzero_id, other_last_nonzero_id);

            // realloc the coefficients if needed
            // quite likely to enter in this if
            if self.coefs.len() < ((largest_nonzero_coef_id + 1) - smallest_pivot) as usize {
                self.coefs.resize(
                    ((largest_nonzero_coef_id + 1) - smallest_pivot) as usize,
                    0u8,
                ); // new coefs will have value 0
            }
            // now we're sure than self.coefs is large enough to contain every nonzero coef

            // the coefs still need to be aligned (here they are only aligned for sure if smallest_pivot == other_bounds.pivot
            // Let's first put the new pivot at index 0 of the array
            if self.constant_term.first_id < smallest_pivot {
                self.coefs
                    .rotate_left((smallest_pivot - self.constant_term.first_id) as usize);
            } else if self.constant_term.first_id > smallest_pivot {
                self.coefs
                    .rotate_right((self.constant_term.first_id - smallest_pivot) as usize);
                // rotate will place the trailing zeroes at the beginning so everything is fine
            }
            self.constant_term.first_id = smallest_pivot;
            self.constant_term.n_protected_symbols = largest_nonzero_coef_id + 1 - smallest_pivot;
            // self.coefs.resize(self.constant_term.n_protected_symbols as usize, 0);
            self.constant_term.n_protected_symbols = self.coefs.len() as u32;
            // now we're back at a normal state where self.constant_term.first_id is at index 0 of self.coefs
            // in addition, pivot is at index 0 of self.coefs (which does not necessarily need to be the case)

            // here self.coef is at least (largest_nonzero_coef_id + 1) - smallest_pivot long and at least as large as other.coefs

            self.constant_term.add(&other.constant_term);
            let other_pivot_index_in_self = self._get_coef_index(other_pivot)?;
            let other_pivot_index_in_other = other._get_coef_index(other_pivot)?;
            arithmetic::xor(
                &mut self.coefs[other_pivot_index_in_self..],
                &other.coefs[other_pivot_index_in_other..],
            );
            self._recompute_bounds(smallest_pivot, largest_nonzero_coef_id);
        }
        if self.bounds == EmptyBounds {
            self.coefs.clear();
            self.constant_term.n_protected_symbols = 0;
        }
        assert_eq!(
            self.constant_term.n_protected_symbols as usize,
            self.coefs.len(),
            "n_protected_symbols is not equal to the coefs vector length !"
        );
        Ok(())
    }

    pub fn mul(&mut self, coef: u8) {
        arithmetic::mul(&mut self.coefs[..], coef);
        self.constant_term.mul(coef);
    }

    pub fn div(&mut self, coef: u8) {
        arithmetic::div(&mut self.coefs[..], coef);
        self.constant_term.div(coef);
    }

    pub fn add_mul(&mut self, coef: u8, other: &Equation) {
        assert!(
            self.coefs.len() == other.coefs.len(),
            "equation coefs mismatch !"
        );
        arithmetic::add_mul(&mut self.coefs[..], coef, &other.coefs[..]);
        self.constant_term.add_mul(coef, &other.constant_term);
    }

    pub fn bounds(&self) -> &EquationBounds {
        &self.bounds
    }

    pub fn is_zero(&self) -> bool {
        match self.bounds() {
            Bounds {
                pivot,
                last_nonzero_id,
            } => pivot > last_nonzero_id,
            EmptyBounds => true,
        }
    }

    pub fn normalize_pivot(&mut self) {
        match self.bounds() {
            Bounds {
                pivot,
                last_nonzero_id,
            } => {
                let pivot_coef = self.get_coef(*pivot);
                if pivot_coef != 1 {
                    self.div(pivot_coef);
                }
            },
            EmptyBounds => {},
        }
    }

    fn _recompute_min_bound(&mut self, from: SymbolID) -> bool {
        let from = std::cmp::max(from, self.constant_term.first_id);
        for (i, coef) in self.coefs[(from - self.constant_term.first_id) as usize..]
            .iter()
            .enumerate()
        {
            let idx = i as u32;
            if *coef != 0u8 {
                let new_pivot = idx + self.constant_term.first_id;
                match &mut self.bounds {
                    Bounds { pivot, .. } => {
                        *pivot = new_pivot;
                    }
                    EmptyBounds => {
                        self.bounds = Bounds {
                            pivot: new_pivot,
                            last_nonzero_id: new_pivot,
                        };
                    }
                }
                return true;
            }
        }
        self.bounds = EmptyBounds;
        return false;
    }

    pub fn recompute_min_bound(&mut self) -> bool {
        self._recompute_min_bound(self.constant_term.first_id)
    }

    pub fn shrink_min_bound(&mut self) -> bool {
        if let Bounds {
            pivot: bounds_pivot,
            ..
        } = self.bounds
        {
            let pivot = bounds_pivot;
            self._recompute_min_bound(pivot)
        } else {
            false
        }
    }

    // from included
    fn _recompute_max_bound(&mut self, from: SymbolID) -> bool {
        let from = std::cmp::min(
            from,
            self.constant_term.first_id + self.constant_term.n_protected_symbols - 1,
        );
        for (i, coef) in self.coefs[..(from + 1 - self.constant_term.first_id) as usize]
            .iter()
            .enumerate()
            .rev()
        {
            let idx = i as u32;
            if *coef != 0u8 {
                let last = idx + self.constant_term.first_id;
                match &mut self.bounds {
                    Bounds {
                        last_nonzero_id: bounds_last,
                        ..
                    } => {
                        *bounds_last = last;
                    }
                    EmptyBounds => {
                        self.bounds = Bounds {
                            pivot: last,
                            last_nonzero_id: last,
                        };
                    }
                }
                return true;
            }
        }
        self.bounds = EmptyBounds;
        return false;
    }

    pub fn recompute_max_bound(&mut self) -> bool {
        self._recompute_max_bound(
            self.constant_term.first_id + self.constant_term.n_protected_symbols - 1,
        )
    }

    pub fn shrink_max_bound(&mut self) -> bool {
        if let Bounds {
            last_nonzero_id: bounds_last_nonzero_id,
            ..
        } = &self.bounds
        {
            let last_nonzero_id = *bounds_last_nonzero_id;
            self._recompute_max_bound(last_nonzero_id)
        } else {
            false
        }
    }

    fn _recompute_bounds(
        &mut self,
        lower_bound_included: SymbolID,
        upper_bound_included: SymbolID,
    ) -> bool {
        let result_max = self._recompute_max_bound(upper_bound_included);
        let result_min = self._recompute_min_bound(lower_bound_included);

        assert_eq!(result_max, result_min, "error in bounds adjusting");
        return result_min;
    }

    pub fn recompute_bounds(&mut self) -> bool {
        self._recompute_bounds(
            self.constant_term.first_id,
            self.constant_term.first_id + self.constant_term.n_protected_symbols - 1,
        )
    }

    // pre: the equation is consistent: there is no nonzero coef outside self.bounds
    pub fn shrink_bounds(&mut self) -> bool {
        let result_max = self.shrink_max_bound();
        let result_min = self.shrink_min_bound();

        assert_eq!(result_max, result_min, "error in bounds adjusting");
        return result_min;
    }

    pub fn get_coef(&self, id: SymbolID) -> u8 {
        let index = self._get_coef_index(id);
        match index {
            Ok(index) => self.coefs[index],
            Err(()) => 0,
        }
    }

    pub fn set_coef(&mut self, id: SymbolID, coef: u8) -> Result<(), ()> {
        let index = self._get_coef_index(id);
        match index {
            Ok(index) => {
                self.coefs[index] = coef;
                Ok(())
            }
            Err(()) => Err(()),
        }
    }

    fn _get_coef_index(&self, id: SymbolID) -> Result<usize, ()> {

        if id < self.constant_term.first_id {
            return Err(());
        }
        let index = (id - self.constant_term.first_id) as usize;
        if index >= self.coefs.len() {
            Err(())
        } else {
            Ok(index)
        }
    }

    pub fn has_one_id(&self) -> bool {
        match self.bounds() {
            Bounds {
                pivot,
                last_nonzero_id,
            } => pivot == last_nonzero_id,
            EmptyBounds => false,
        }
    }

    fn clone(&self) -> Equation {
        return Equation {
            coefs: self.coefs.to_vec(),
            constant_term: Symbol::new(
                self.constant_term.first_id,
                self.constant_term.n_protected_symbols,
                self.constant_term.get_data().to_vec(),
            ),
            bounds: match self.bounds {
                Bounds {
                    pivot,
                    last_nonzero_id,
                } => Bounds {
                    pivot,
                    last_nonzero_id,
                },
                EmptyBounds => EmptyBounds,
            },
        };
    }
}

#[cfg(test)]
pub fn eq_clone(eq: &Equation) -> Equation {
    return eq.clone();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::*;
    use rand;
    use rand::distributions::{Distribution, Standard};
    use rand::SeedableRng;
    use rand_pcg::{Lcg64Xsh32, Pcg32};
    use std::cmp::min;

    use crate::tests::*;


    // shifts the array to have the first nonzero value at index 0 and returns a tuple composed of :
    //  - the new upper nonzero bounds (the lower one is always 0 or None if upper is None)
    //  - the number of shifts
    fn get_bounds_and_shift(vec: &mut Vec<u8>) -> (Option<usize>, usize) {
        let lower = (vec).into_iter().position(|x| *x != 0);

        let upper = (vec).into_iter().rposition(|x| *x != 0);

        if let (Some(lower), Some(upper)) = (lower, upper) {
            if lower != 0 {
                vec.rotate_left(lower);
            }

            (Some(upper - lower), lower)
        } else {
            (None, 0)
        }
    }

    fn check_eq_coefs(eq: &Equation, reference_coefs: &Vec<u8>, coefs_first_id: SymbolID) {
        let mut n_assessed_indices = 0;

        // when coefs_first_id < eq.constant_term.first_id
        for id in coefs_first_id..eq.constant_term.first_id {
            // ensure that the additional coefs present in reference_coefs are all 0
            assert_eq!(0, reference_coefs[(id - coefs_first_id) as usize]);
            n_assessed_indices += 1;
        }

        // when coefs_first_id > eq.constant_term.first_id
        for id in eq.constant_term.first_id..coefs_first_id {
            // ensure that the additional coefs present in the equation are all 0
            assert_eq!(0, eq.coefs[(id - eq.constant_term.first_id) as usize]);
            // also check get_coef just in case
            assert_eq!(0, eq.get_coef(id));
            n_assessed_indices += 1;
        }

        let (eq_start_index, reference_start_index) = if eq.constant_term.first_id < coefs_first_id
        {
            (
                coefs_first_id as usize - eq.constant_term.first_id as usize,
                0 as usize,
            )
        } else {
            // we already checked the coefs out of bounds
            (
                0 as usize,
                eq.constant_term.first_id as usize - coefs_first_id as usize,
            )
        };

        let min_len = std::cmp::min(eq.coefs.len(), reference_coefs.len())
            - std::cmp::max(eq_start_index, reference_start_index);
        let max_len = std::cmp::max(eq.coefs.len(), reference_coefs.len())
            - std::cmp::max(eq_start_index, reference_start_index);

        for i in 0..min_len {
            assert_eq!(
                eq.coefs[eq_start_index + i],
                reference_coefs[reference_start_index + i],
                "coefs are not equal at index {}",
                i
            );
            n_assessed_indices += 1;
        }
        for i in min_len..max_len {
            if eq.coefs.len() <= eq_start_index + i {
                assert_eq!(0, reference_coefs[reference_start_index + i]);
            } else {
                assert_eq!(0, eq.coefs[eq_start_index + i]);
            }
            n_assessed_indices += 1;
        }

        // ensure that the test verifies every possible coef (this is like a test for the test...)
        assert_eq!(
            std::cmp::max(eq.coefs.len(), reference_coefs.len()),
            n_assessed_indices
        );
    }

    fn check_eq_bounds(eq: &Equation) {
        match eq.bounds {
            Bounds {
                pivot: eq_pivot,
                last_nonzero_id: eq_last_nonzero_id,
            } => {
                for id in eq.constant_term.first_id..eq_pivot {
                    let index = (id - eq.constant_term.first_id) as usize;
                    assert_eq!(0, eq.get_coef(id));
                    // double check for the same thing, not necessarily useful but we do it anyway
                    // we really want to ensure the integrity of the array because we may
                    // do batch operations with specific CPU instructions so sometimes we
                    // won't use get_coef in the future
                    assert_eq!(0, eq.coefs[index]);
                }
                for id in (eq_last_nonzero_id + 1)
                    ..eq.constant_term.first_id + eq.constant_term.n_protected_symbols
                {
                    let index = (id - eq.constant_term.first_id) as usize;
                    assert_eq!(0, eq.get_coef(id));
                    // double check for the same thing, not necessarily useful but we do it anyway
                    // we really want to ensure the integrity of the array because we may
                    // do batch operations with specific CPU instructions so sometimes we
                    // won't use get_coef in the future
                    assert_eq!(0, eq.coefs[index]);
                }
                // eq.constant_term.n_protected_symbols and eq.coefs.len() should be equal but it may change in the future
                for index in eq.constant_term.n_protected_symbols as usize..eq.coefs.len() {
                    assert_eq!(0, eq.coefs[index]);
                }
            }
            EmptyBounds => {
                for coef in eq.coefs.iter() {
                    assert_eq!(0, *coef);
                }
            }
        }
    }

    fn check_equation(
        eq: &Equation,
        xored_coefs: &mut Vec<u8>,
        coefs_first_id: SymbolID,
        xored_constant_term: &Vec<u8>,
    ) {
        let (computed_last_nonzero_index, nshifts) = get_bounds_and_shift(xored_coefs);
        let new_first_id_and_pivot = coefs_first_id + nshifts as SymbolID;
        if let Some(last) = computed_last_nonzero_index {
            if let Bounds {
                pivot: eq_pivot,
                last_nonzero_id: eq_last_nonzero_id,
            } = eq.bounds
            {
                assert_eq!(eq_pivot, new_first_id_and_pivot);
                assert_eq!(
                    eq_last_nonzero_id,
                    last as SymbolID + new_first_id_and_pivot
                );

                let xored_coefs: &Vec<u8> = xored_coefs.as_ref();
                check_eq_coefs(eq, xored_coefs, coefs_first_id + nshifts as SymbolID);
                assert_eq!(xored_constant_term, eq.constant_term.get_data());
            } else {
                panic!();
            }
        } else {
            assert_eq!(eq.bounds, EmptyBounds);
        }
        check_eq_bounds(eq);
    }

    #[test]
    fn test_aligned_equations_add() {
        let mut rng = get_new_rng();
        let mut eq1 = get_equation(&mut rng, 10, 20, 13, 18);
        let eq2 = get_equation(&mut rng, 10, 20, 13, 18);
        let eq1_clone = eq1.clone();

        let mut xored_coefs = eq1_clone.coefs.to_vec();
        let mut xored_constant_term = eq1_clone.constant_term.get_data().to_vec();
        arithmetic::xor(&mut xored_coefs, &eq2.coefs);
        arithmetic::xor(&mut xored_constant_term, &eq2.constant_term.get_data());

        assert!(matches!(eq1.add(&eq2), Ok(())));

        check_equation(
            &eq1,
            &mut xored_coefs,
            eq1_clone.constant_term.first_id,
            &xored_constant_term,
        );
    }

    #[test]
    fn test_unaligned_inside_equations_add() {
        let mut rng = get_new_rng();
        // eq1's array spans from 10 to 59 included and its nonzero bounds span from 13 to 55
        let mut eq1 = get_equation(&mut rng, 10, 50, 13, 55);
        // eq2's array spans from 28 to 50 included and its nonzero bounds span from 32 to 48
        let eq2 = get_equation(&mut rng, 28, 23, 32, 48);
        let eq1_clone = eq1.clone();

        let mut xored_coefs = eq1_clone.coefs.to_vec();
        let mut xored_constant_term = eq1_clone.constant_term.get_data().to_vec();

        // 0..17 and 51..60 stay unchanged
        arithmetic::xor(&mut xored_coefs[18..41], &eq2.coefs);
        arithmetic::xor(&mut xored_constant_term, &eq2.constant_term.get_data());

        assert!(matches!(eq1.add(&eq2), Ok(())));

        check_equation(
            &eq1,
            &mut xored_coefs,
            eq1_clone.constant_term.first_id,
            &xored_constant_term,
        );
    }

    #[test]
    fn test_unaligned_outside_left_equations_add() {
        let mut rng = get_new_rng();
        // eq1's array spans from 10 to 59 included and its nonzero bounds span from 13 to 55
        let mut eq1 = get_equation(&mut rng, 10, 50, 13, 55);
        // eq2's array spans from 1 to 30 included and its nonzero bounds span from 3 to 16
        let eq2 = get_equation(&mut rng, 1, 30, 3, 16);
        let eq1_clone = eq1.clone();

        let mut xored_coefs = eq1_clone.coefs.to_vec(); // eq2 is larger so eq2's coefs are the base for our xor
        xored_coefs.resize(100, 0);
        xored_coefs.rotate_right(9);
        // 0..8 and 55..150 stay unchanged
        arithmetic::xor(&mut xored_coefs, &eq2.coefs);
        let mut xored_constant_term = eq1_clone.constant_term.get_data().to_vec();

        arithmetic::xor(&mut xored_constant_term, &eq2.constant_term.get_data());

        assert!(matches!(eq1.add(&eq2), Ok(())));

        check_equation(
            &eq1,
            &mut xored_coefs,
            eq2.constant_term.first_id,
            &xored_constant_term,
        );
    }

    #[test]
    fn test_unaligned_outside_equations_add() {
        let mut rng = get_new_rng();
        // eq1's array spans from 10 to 59 included and its nonzero bounds span from 13 to 55
        let mut eq1 = get_equation(&mut rng, 10, 50, 13, 55);
        // eq2's array spans from 1 to 150 included and its nonzero bounds span from 3 to 137
        let eq2 = get_equation(&mut rng, 1, 150, 3, 137);
        let eq1_clone = eq1.clone();

        let mut xored_coefs = eq2.coefs.to_vec(); // eq2 is larger so eq2's coefs are the base for our xor
                                                  // 0..8 and 55..150 stay unchanged
        arithmetic::xor(&mut xored_coefs[9..55], &eq1.coefs);
        let mut xored_constant_term = eq1_clone.constant_term.get_data().to_vec();

        arithmetic::xor(&mut xored_constant_term, &eq2.constant_term.get_data());

        assert!(matches!(eq1.add(&eq2), Ok(())));

        check_equation(
            &eq1,
            &mut xored_coefs,
            eq2.constant_term.first_id,
            &xored_constant_term,
        );
    }

    #[test]
    fn test_unaligned_equations_creating_zeroes() {
        let mut rng = get_new_rng();
        let mut eq1 = get_equation(&mut rng, 10, 150, 13, 140);
        let mut eq2 = get_equation(&mut rng, 4, 200, 13, 140);
        let eq1_clone = eq1.clone();

        if let Bounds { .. } = eq1.bounds() {
            // make the start and end of eq2 cancelling the coefs of eq1
            for id in 13..50 {
                // put the same coef in eq1 and eq2 so that the add cancels them
                assert!(matches!(eq2.set_coef(id, eq1_clone.get_coef(id)), Ok(())));
            }
            for id in 120..=140 {
                // put the same coef in eq1 and eq2 so that the add cancels them
                assert!(matches!(eq2.set_coef(id, eq1_clone.get_coef(id)), Ok(())));
            }
        } else {
            panic!("eq1 should have nonzero bounds !");
        }

        let mut xored_coefs = eq2.coefs.to_vec(); // eq2 is larger so eq2's coefs are the base for our xor
        arithmetic::xor(&mut xored_coefs[6..156], &eq1_clone.coefs);

        let mut xored_constant_term = eq1_clone.constant_term.get_data().to_vec();
        arithmetic::xor(&mut xored_constant_term, &eq2.constant_term.get_data());

        assert!(matches!(eq1.add(&eq2), Ok(())));
        if let Bounds {
            pivot: eq1_pivot,
            last_nonzero_id: eq1_last_nonzero_id,
        } = eq1.bounds()
        {
            assert_eq!(50, *eq1_pivot);
            assert_eq!(119, *eq1_last_nonzero_id);
        } else {
            panic!("eq1 should have nonzero bounds after addition !")
        }

        check_equation(
            &eq1,
            &mut xored_coefs,
            eq2.constant_term.first_id,
            &xored_constant_term,
        );
    }
}
