

use super::*;
use rand;
use rand::distributions::{Distribution, Standard};
use rand::SeedableRng;
use rand_pcg::{Lcg64Xsh32, Pcg32};
use std::cmp::min;

use crate::symbols::{gf_tables, Symbol, SymbolID};
use crate::system::equation::Equation;
use crate::system::equation::EquationBounds::Bounds;

#[cfg(test)]

pub fn get_new_rng() -> Lcg64Xsh32 {
    Pcg32::seed_from_u64(42)
}

pub fn get_random_vec(rng: &mut Lcg64Xsh32, size: usize) -> Vec<u8> {
    Standard.sample_iter(rng).take(size).collect()
}

const SYMBOL_SIZE: usize = 1500;

pub fn get_equation(
    rng: &mut Lcg64Xsh32,
    first_id: SymbolID,
    n_protected_symbols: u32,
    first_nonzero_id: SymbolID,
    last_nonzero_id: SymbolID,
) -> Equation {
    let constant_term = get_random_vec(rng, SYMBOL_SIZE);
    let mut coefs = get_random_vec(rng, n_protected_symbols as usize);

    for id in first_id..first_nonzero_id {
        coefs[(id - first_id) as usize] = 0;
    }
    let first_nonzero_index = (first_nonzero_id - first_id) as usize;
    // ensure the pivot is not randomly zero
    if coefs[first_nonzero_index] == 0 {
        coefs[first_nonzero_index] = 1;
    }

    for id in last_nonzero_id + 1..first_id + n_protected_symbols {
        coefs[(id - first_id) as usize] = 0;
    }

    let last_nonzero_index = (last_nonzero_id - first_id) as usize;
    // ensure the last nonzero coef is not randomly zero
    if coefs[last_nonzero_index] == 0 {
        coefs[last_nonzero_index] = 1;
    }

    Equation::_new(coefs,
                   Symbol::new(first_id, n_protected_symbols, constant_term),
                   Bounds {
                        pivot: first_nonzero_id,
                        last_nonzero_id,
                   })
}