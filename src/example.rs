use std::time::Instant;

use rustgf::system::System;
use rand::distributions::{Distribution, Standard};
use rand::SeedableRng;
use rand_pcg::{Lcg64Xsh32, Pcg32};
use rustgf::symbols::{Symbol, SymbolID};
use rustgf::system::equation::{Equation};



fn get_new_rng() -> Lcg64Xsh32 {
    Pcg32::seed_from_u64(42)
}

fn get_random_vec(rng: &mut Lcg64Xsh32, size: usize) -> Vec<u8> {
    Standard.sample_iter(rng).take(size).collect()
}

fn generate_equation_from_payloads(rng: &mut Lcg64Xsh32, first_coef_id: SymbolID, n_coefs: u32, symbols: &Vec<Symbol>, symbol_size: usize) -> Equation {
    let coefs = get_random_vec(rng, n_coefs as usize);
    let mut out_symbol = Symbol::new(first_coef_id, n_coefs, vec![0; symbol_size]);
    for (i, _) in symbols.iter().enumerate() {
        out_symbol.add_mul(coefs[i], &symbols[i]);
    }
    Equation::new(coefs, out_symbol)
}

fn add_several_full_equations() {
    let mut rng = get_new_rng();
    let mut system = System::new();
    let first_id = 40 as SymbolID;
    let n_symbols = 100;
    let symbol_size = 1500;
    let mut symbols = Vec::new();
    for _ in 0..n_symbols {
        symbols.push(Symbol::new(first_id, n_symbols, get_random_vec(&mut rng, symbol_size)));
    }

    let mut equations = Vec::new();

    for _ in 0..n_symbols {
        equations.push(generate_equation_from_payloads(&mut rng, first_id, n_symbols, &symbols, symbol_size));
    }

    let mut decoded_ids: Vec<SymbolID> = Vec::new();


    let start = Instant::now();
    for equation in equations {
        let result = system.add(equation);
        let (_, mut decoded) = result.unwrap();
        decoded_ids.append(&mut decoded);
    }
    let duration = start.elapsed();
    println!("elapsed {:?}", duration);

    assert_eq!(decoded_ids.len(), n_symbols as usize);

    for id in decoded_ids {
        let solved = system.take(id).unwrap();
        assert_eq!(&solved, symbols[id as usize - first_id as usize].get_data());
    }
}

fn main() {
    add_several_full_equations();
}
