pub mod arithmetic;
pub mod gf_tables;

pub type SymbolID = u64;

#[derive(Debug,Clone)]
pub struct Symbol {
    pub first_id: SymbolID,
    pub n_protected_symbols: u64, // states how many symbols have been combined to create this symbol
    data: Vec<u8>,
}

impl Symbol {
    pub fn new(first_id: SymbolID, n_protected_symbols: u64, data: Vec<u8>) -> Symbol {
        Symbol {
            first_id,
            n_protected_symbols,
            data,
        }
    }

    pub fn add(&mut self, other: &Symbol, field: Option<&galois_2p8::PrimitivePolynomialField>) {
        assert!(
            self.data.len() == other.data.len(),
            "symbols size mismatch !"
        );
        arithmetic::xor(&mut self.data, &other.data, field);
    }

    pub fn mul(&mut self, coef: u8, field: Option<&galois_2p8::PrimitivePolynomialField>) {
        arithmetic::mul(&mut self.data, coef, field);
    }

    pub fn div(&mut self, coef: u8, field: Option<&galois_2p8::PrimitivePolynomialField>) {
        arithmetic::div(&mut self.data, coef, field);
    }

    pub fn add_mul(&mut self, coef: u8, other: &Symbol, field: Option<&galois_2p8::PrimitivePolynomialField>) {
        assert!(
            self.data.len() == other.data.len(),
            "symbols size mismatch !"
        );
        arithmetic::add_mul(&mut self.data, coef, &other.data, field);
    }

    pub fn first_id(&self) -> SymbolID {
        self.first_id
    }

    pub fn last_id(&self) -> SymbolID {
        self.first_id + self.n_protected_symbols - 1
    }

    pub fn n_protected_symbols(&self) -> SymbolID {
        self.n_protected_symbols as SymbolID
    }

    pub fn get_data(&self) -> &Vec<u8> {
        &self.data
    }

    pub fn take_data(self) -> Vec<u8> {
        self.data
    }
}
