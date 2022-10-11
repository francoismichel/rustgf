use std::cmp;

use crate::symbols::*;

pub fn xor(a: &mut [u8], b: &[u8]) {
    let min_len = cmp::min(a.len(), b.len());
    for i in 0..min_len {
        a[i] = a[i] ^ b[i];
    }
}

pub fn mul(a: &mut [u8], coef: u8) {
    for i in 0..a.len() {
        a[i] = gf_tables::tables::GF256_MUL_TABLE[usize::from(a[i])][usize::from(coef)];
    }
}

pub fn div(a: &mut [u8], coef: u8) {
    mul(a, gf_tables::tables::GF256_INV_TABLE[usize::from(coef)])
}

pub fn add_mul(a: &mut [u8], coef: u8, b: &[u8]) {
    let min_len = cmp::min(a.len(), b.len());
    for i in 0..min_len {
        a[i] = a[i] ^ gf_tables::tables::GF256_MUL_TABLE[usize::from(coef)][usize::from(b[i])];
    }
}

pub fn pow(mut a: u8, mut exp: u8) -> u8 {
    let mut res = 1;
    while exp > 0 {
        if exp & 1 == 0 {
            a = gf_tables::tables::GF256_MUL_TABLE[a as usize][a as usize];
            exp >>= 1;
        } else {
            res = gf_tables::tables::GF256_MUL_TABLE[res as usize][a as usize];
            exp -= 1;
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use crate::symbols::{arithmetic, gf_tables};

    #[test]
    fn test_pow() {
        assert_eq!(arithmetic::pow(42, 4), gf_tables::mul(gf_tables::mul(42, 42), gf_tables::mul(42, 42)));
    }
}
