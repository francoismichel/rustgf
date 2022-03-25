pub mod arithmetic;
pub mod gf_tables;

pub type SymbolID = u32;

#[derive(Debug)]
pub struct Symbol {
    pub first_id: SymbolID,
    pub n_protected_symbols: u32, // states how many symbols have been combined to create this symbol
    data: Vec<u8>,
}

impl Symbol {
    pub fn new(first_id: SymbolID, n_protected_symbols: u32, data: Vec<u8>) -> Symbol {
        Symbol {
            first_id,
            n_protected_symbols,
            data,
        }
    }

    pub fn add(&mut self, other: &Symbol) {
        assert!(
            self.data.len() == other.data.len(),
            "symbols size mismatch !"
        );
        arithmetic::xor(&mut self.data, &other.data);
    }

    pub fn mul(&mut self, coef: u8) {
        arithmetic::mul(&mut self.data, coef);
    }

    pub fn div(&mut self, coef: u8) {
        arithmetic::div(&mut self.data, coef);
    }

    pub fn add_mul(&mut self, coef: u8, other: &Symbol) {
        assert!(
            self.data.len() == other.data.len(),
            "symbols size mismatch !"
        );
        arithmetic::add_mul(&mut self.data, coef, &other.data);
    }

    pub fn first_id(&self) -> SymbolID {
        self.first_id
    }

    pub fn last_id(&self) -> SymbolID {
        self.first_id + self.n_protected_symbols - 1
    }

    pub fn n_protected_symbols(&self) -> SymbolID {
        self.n_protected_symbols
    }

    pub fn get_data(&self) -> &Vec<u8> {
        &self.data
    }

    pub fn take_data(self) -> Vec<u8> {
        self.data
    }
}
