use std::cmp;

use galois_2p8::Field;

use crate::symbols::*;

pub fn xor(a: &mut [u8], b: &[u8], field: Option<&galois_2p8::PrimitivePolynomialField>) {
    if let Some(field) = field {
        field.add_multiword(a, b);
    } else {
        a.iter_mut()
            .zip(b.iter())
            .for_each(|(x1, x2)| *x1 ^= *x2);
    }
}


pub fn mul(a: &mut [u8], coef: u8, field: Option<&galois_2p8::PrimitivePolynomialField>) {
    if let Some(field) = field {
        field.mult_multiword(a, coef);
    } else {
        for i in 0..a.len() {
            a[i] = gf_tables::tables::GF256_MUL_TABLE[usize::from(a[i])][usize::from(coef)];
        }
    }
}

pub fn div(a: &mut [u8], coef: u8, field: Option<&galois_2p8::PrimitivePolynomialField>) {
    mul(a, gf_tables::tables::GF256_INV_TABLE[usize::from(coef)], field)
}

pub fn add_mul(a: &mut [u8], coef: u8, b: &[u8], field: Option<&galois_2p8::PrimitivePolynomialField>) {
    if let Some(field) = field {
        field.add_scaled_multiword(a, b, coef);
    } else {
        let min_len = cmp::min(a.len(), b.len());
        for i in 0..min_len {
            a[i] ^= gf_tables::tables::GF256_MUL_TABLE[usize::from(coef)][usize::from(b[i])];
        }
    }
}
