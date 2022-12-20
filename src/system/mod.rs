pub mod equation;

use crate::symbols::{gf_tables};
use crate::symbols::SymbolID;

use crate::system::equation::EquationBounds;
use crate::system::SystemBounds::EmptyBounds;
use equation::Equation;
use crate::system::SystemError::{InternalError, UnusedEquation};

#[derive(Debug, PartialEq, Eq)]
pub enum SystemBounds {
    Bounds {
        first_equation_pivot_id: SymbolID,  // this first equation may not be present in the system
        last_present_equation_pivot_id: SymbolID,
        largest_nonzero_id: SymbolID,
    },
    EmptyBounds,
}

#[derive(Debug)]
pub enum SystemError {
    UnusedEquation,
    InternalError(String),
}

#[derive(Debug)]
pub struct System {
    bounds: SystemBounds,
    max_equations: u64,
    n_equations: u64,
    equations: Vec<Option<Equation>>,
}

impl System {
    pub fn new() -> System {
        let max_equations = 10000u64;
        let mut equations = Vec::new();
        equations.resize_with(max_equations as usize, || None);
        System {
            bounds: SystemBounds::EmptyBounds,
            n_equations: 0,
            max_equations,
            equations,
        }
    }

    pub fn bounds(&self) -> &SystemBounds {
        &self.bounds
    }

    pub fn is_empty(&self) -> bool {
        self.bounds == SystemBounds::EmptyBounds
    }

    // the pivot must be valid !
    fn _adjust_last_nonzero_id_bound_from(&mut self, from: SymbolID) {
        match &mut self.bounds {
            SystemBounds::Bounds {
                first_equation_pivot_id,
                last_present_equation_pivot_id,
                largest_nonzero_id: _,
            } => {
                for id in (*first_equation_pivot_id..from).rev() {
                    if let Some(_) = self.equations[(id - *first_equation_pivot_id) as usize] {
                        *last_present_equation_pivot_id = id;
                        return;
                    }
                }
            }
            EmptyBounds => {
                return;
            }
        }
    }

    fn _compute_largest_symbol_id(&self) -> Result<SymbolID, ()> {
        let mut retval = None;
        let mut equations_checked = 0u64;

        for maybe_eq in &self.equations {
            if let Some(eq) = maybe_eq {
                if let EquationBounds::Bounds {
                    last_nonzero_id: eq_last_nonzero_id,
                    ..
                } = eq.bounds()
                {
                    retval = match retval {
                        Some(val) => Some(std::cmp::max(*eq_last_nonzero_id, val)),
                        Option::None => Some(*eq_last_nonzero_id),
                    };
                    equations_checked += 1;
                    if equations_checked == self.n_equations {
                        break;
                    }
                }
            }
        }

        match retval {
            Some(val) => Ok(val),
            Option::None => Err(()),
        }
    }

    // @brief Add an equation to the system.
    //
    // Returns Some(eq) will be returned if `eq` was removed fom the system as it was occupying the slot of `equation`
    fn add_pivot_equation(&mut self, equation: Equation) -> Result<Option<Equation>, ()> {
        let mut removed = None;
        if let EquationBounds::Bounds {
            pivot: equation_pivot,
            last_nonzero_id: equation_last_nonzero_id,
        } = equation.bounds()
        {
            let candidate_pivot = *equation_pivot;
            let candidate_last_nonzero_id = *equation_last_nonzero_id;
            /* The added symbol has first nonzero index outside the system */
            match self.bounds {
                SystemBounds::Bounds {
                    first_equation_pivot_id: system_first_equation_pivot_id,
                    last_present_equation_pivot_id: system_last_present_equation_pivot_id,
                    largest_nonzero_id: system_largest_nonzero_id,
                } => {
                    let system_pivot = system_first_equation_pivot_id;
                    if candidate_pivot < system_pivot {
                        let mut new_last_nonzero_equation_pivot = candidate_pivot; // starts at the new pivot, but will likely grow
                        let mut new_largest_nonzero_id =
                            std::cmp::max(candidate_last_nonzero_id, system_largest_nonzero_id);
                        let new_pivot = candidate_pivot;
                        if system_pivot - candidate_pivot < self.max_equations {
                            // self.equations.insert(0, Some(equation));
                            // FIXME: inefficient but should be rare
                            // there are system_pivot - candidate_pivot places to create between the new pivot and the old one
                            // we delete the old ones if needed
                            let places_to_create = system_pivot - candidate_pivot;

                            if places_to_create
                                > self.max_equations
                                    - (system_last_present_equation_pivot_id
                                        - system_first_equation_pivot_id)
                            {
                                // there are existing equations that will be shifted out to make new place, so remove them
                                for maybe_eq in self.equations[((self.max_equations
                                    - places_to_create)
                                    as usize)
                                    ..(self.max_equations as usize)]
                                    .iter_mut()
                                {
                                    if let Some(_) = maybe_eq {
                                        // we remove this equation
                                        *maybe_eq = Option::None;
                                        self.n_equations -= 1;
                                    }
                                }
                                // update the largest_nonzero_id if needed
                                let res = self._compute_largest_symbol_id();
                                new_largest_nonzero_id = match res {
                                    Ok(val) => std::cmp::max(candidate_last_nonzero_id, val),
                                    Err(_) => candidate_last_nonzero_id,
                                }
                            } else {
                                new_largest_nonzero_id = std::cmp::max(
                                    system_largest_nonzero_id,
                                    candidate_last_nonzero_id,
                                );
                            }
                            // now the area between self.max_equations - places_to_create and self.max_equations is only composed of None
                            // now move the previously present equations
                            let last_possibly_nonzero_index = std::cmp::min(
                                system_last_present_equation_pivot_id
                                    - system_first_equation_pivot_id,
                                self.max_equations - places_to_create,
                            );

                            // let's shift all equations places_to_create places to the right
                            let mut moved_equations = 0;
                            for i in (0..last_possibly_nonzero_index).rev() {
                                if let Some(_) = self.equations[i as usize] {
                                    let new_index = i + places_to_create;
                                    // "rust way" to shift the elem at i to new_index
                                    self.equations[new_index as usize] =
                                        self.equations[i as usize].take();
                                    if (new_last_nonzero_equation_pivot)
                                        < (candidate_pivot) + (new_index)
                                    {
                                        new_last_nonzero_equation_pivot =
                                            candidate_pivot + new_index;
                                    }
                                    moved_equations -= 1;
                                    if moved_equations == self.n_equations {
                                        // no remaining equation
                                        break;
                                    }
                                }
                            }
                            // everything has shifted, now insert the new pivot
                        } else {
                            // wipe the whole system and add the equation as the only equation
                            for val in self.equations.iter_mut() {
                                *val = None
                            }
                            self.n_equations = 0;
                            new_last_nonzero_equation_pivot = new_pivot;
                        }
                        if let SystemBounds::Bounds {
                            first_equation_pivot_id,
                            last_present_equation_pivot_id,
                            largest_nonzero_id,
                        } = &mut self.bounds
                        {
                            *first_equation_pivot_id = new_pivot;
                            *last_present_equation_pivot_id = new_last_nonzero_equation_pivot;
                            *largest_nonzero_id = new_largest_nonzero_id;
                        }
                        self.n_equations += 1;
                        self.equations[0] = Some(equation);
                    } else {
                        let idx_pos = candidate_pivot - system_first_equation_pivot_id;
                        if idx_pos >= self.max_equations {
                            return Err(());
                        }
                        let eq_to_insert = if let Some(eq_already_present) =
                            self.equations[idx_pos as usize].take()
                        {
                            self.n_equations -= 1;
                            let (to_remove, to_insert) = match eq_already_present.bounds() {
                                EquationBounds::Bounds {
                                    last_nonzero_id: eq_already_present_last_nonzero_id,
                                    pivot: eq_already_present_pivot,
                                } => {
                                    if *eq_already_present_last_nonzero_id
                                        - *eq_already_present_pivot
                                        <= candidate_last_nonzero_id - candidate_pivot
                                    {
                                        // don't need to insert the new equation, the one already present is better
                                        (equation, eq_already_present)
                                    } else {
                                        (eq_already_present, equation)
                                    }
                                }
                                EquationBounds::EmptyBounds => (eq_already_present, equation),
                            };
                            removed = Some(to_remove);
                            to_insert
                        } else {
                            equation
                        };

                        if let SystemBounds::Bounds {
                            largest_nonzero_id,
                            last_present_equation_pivot_id,
                            ..
                        } = &mut self.bounds
                        {
                            *last_present_equation_pivot_id =
                                std::cmp::max(*last_present_equation_pivot_id, candidate_pivot);
                            *largest_nonzero_id =
                                std::cmp::max(*largest_nonzero_id, candidate_last_nonzero_id);
                        }
                        self.equations[idx_pos as usize] = Some(eq_to_insert);
                        self.n_equations += 1;
                    }
                }
                SystemBounds::EmptyBounds => {
                    self.bounds = SystemBounds::Bounds {
                        first_equation_pivot_id: candidate_pivot,
                        last_present_equation_pivot_id: candidate_pivot,
                        largest_nonzero_id: candidate_last_nonzero_id,
                    };
                    self.equations[0] = Some(equation);
                    self.n_equations = 1;
                }
            }
        }
        self.recompute_bounds();
        Ok(removed)
    }

    fn reduce_equation(&self, new_equation: &mut Equation) -> Result<(), ()> {
        let mut n_non_null_equations = 0;
        for eq_opt in &self.equations {
            if n_non_null_equations >= self.n_equations {
                break;
            }

            if let Some(eq) = eq_opt {
                n_non_null_equations += 1;

                if let EquationBounds::Bounds {
                    pivot: eq_pivot,
                    last_nonzero_id: _,
                } = eq.bounds() {
                    let mut new_equation_pivot_value = None;
                    if let EquationBounds::Bounds {
                        pivot: new_equation_pivot,
                        last_nonzero_id: _,
                    } = new_equation.bounds() {
                        new_equation_pivot_value = Some(*new_equation_pivot);
                    }
                    if let Some(_) = new_equation_pivot_value {
                        let coef = new_equation.get_coef(*eq_pivot);
                        if coef != 0 {
                            // reduce this coef in the new eq
                            new_equation.mul(gf_tables::mul(gf_tables::inv(coef), eq.get_coef(*eq_pivot)));
                            new_equation.add(eq)?;
                        }
                    }
                }
            }
        }
        Ok(())

    }

    fn recompute_bounds(&mut self) {
        if self.n_equations == 0 {
            self.bounds = SystemBounds::EmptyBounds;
        }
    }

    /// adds eq to the system
    /// without error, a tuple (removed_eq, decoded_symbol_ids) is returned. removed_eq is an equation that
    /// has been removed from the system if it was already occupying new_equation's pivot position.
    /// decoded_symbol_ids are the ids of all the newly decoded equations if they exist
    pub fn add(
        &mut self,
        mut new_equation: Equation,
    ) -> Result<(Option<Equation>, Vec<SymbolID>), SystemError> {
        let mut decoded_symbol_ids = Vec::new();

        new_equation.recompute_bounds();

        match self.reduce_equation(&mut new_equation) {
            Err(()) => {
                return Err(InternalError("could not reduce the new equation".to_string()));
            }
            Ok(()) => ()
        }

        if let EquationBounds::EmptyBounds = new_equation.bounds() {
            return Err(UnusedEquation);
        }

        if let EquationBounds::Bounds {
            pivot: equation_pivot,
            last_nonzero_id: _,
        } = new_equation.bounds()
        {
            let first_id = *equation_pivot;

            let mut n_non_null_equations = 0;
            for eq_opt in &mut self.equations {
                if n_non_null_equations >= self.n_equations {
                    break;
                }

                if let Some(eq) = eq_opt {
                    n_non_null_equations += 1;
                    let coef = eq.get_coef(first_id);
                    if coef != 0 {
                        // the current equation has a non-null coef at the index of the pivot of the new equation
                        // let's reduce it
                        let new_equation_coef = new_equation.get_coef(first_id);
                        new_equation.mul(gf_tables::mul(coef, gf_tables::inv(new_equation_coef)));
                        let has_one_id_before_add = eq.has_one_id();
                        match eq.add(&new_equation){
                            Err(()) => {
                                return Err(InternalError("could not add the equations".to_string()));
                            }
                            Ok(()) => ()
                        }
                        match eq.bounds() {
                            EquationBounds::Bounds {
                                pivot: eq_pivot, ..
                            } => {
                                let is_decoded = !has_one_id_before_add && eq.has_one_id();
                                if is_decoded {
                                    let decoded_id = *eq_pivot;
                                    if eq.get_coef(decoded_id) != 1 {
                                        eq.div(eq.get_coef(decoded_id));
                                    }
                                    decoded_symbol_ids.push(decoded_id);
                                }
                            }
                            EquationBounds::EmptyBounds => return Err(InternalError("a newly added equation totally zeroed the coefficients of another one".to_string())),
                        }
                    }
                }
            }
        }


        if new_equation.has_one_id() {
            new_equation.normalize_pivot();
            if let EquationBounds::Bounds {
                pivot: decoded_pivot,
                ..
            } = new_equation.bounds() {
                let decoded_pivot_value = *decoded_pivot;
                decoded_symbol_ids.push(decoded_pivot_value);
            }
        }

        self.recompute_bounds();

        match self.add_pivot_equation(new_equation) {
            Err(()) => {
                Err(InternalError("error when adding new pivot equation".to_string()))
            }
            Ok(removed) => {
                self.recompute_bounds();

                Ok((removed, decoded_symbol_ids))
            }
        }
    }

    pub fn take(&mut self, id: SymbolID) -> Option<Vec<u8>> {
        if let SystemBounds::Bounds {
            first_equation_pivot_id, last_present_equation_pivot_id, largest_nonzero_id: _
        } = self.bounds {
            if id < first_equation_pivot_id || id > last_present_equation_pivot_id {
                return None;
            }
            let mut retval = self.equations[(id - first_equation_pivot_id) as usize].take()?;
            self.n_equations -= 1;
            if !retval.has_one_id() {
                return None;
            }
            retval.normalize_pivot();
            self.recompute_bounds();
            return Some(retval.constant_term_data());
        }
        return None;
    }

    fn _recompute_bounds_heavy(&mut self) {
        if let SystemBounds::Bounds { first_equation_pivot_id, .. } = self.bounds() {
            let mut new_last_present_equation_pivot_id = None;
            let mut new_largest_nonzero_id = None;
            let mut n_non_null_equations = 0;
            let first_pivot_id = *first_equation_pivot_id;
            for eq_opt in self.equations.iter_mut().rev() {
                if let Some(eq) = eq_opt {
                    if let EquationBounds::Bounds { pivot, last_nonzero_id } = eq.bounds() {
                        if let None = new_last_present_equation_pivot_id {
                            new_last_present_equation_pivot_id = Some(*pivot);
                        }
                        new_largest_nonzero_id = match new_largest_nonzero_id {
                            None => Some(*last_nonzero_id),
                            Some(v) => Some(std::cmp::max(v, *last_nonzero_id)),
                        };
                        n_non_null_equations += 1;
                    } else {
                        eq_opt.take();
                    }
                } 
            }
    
            self.bounds = match (new_largest_nonzero_id, new_last_present_equation_pivot_id) {
                (Some(largest_nonzero), Some(last_present)) => {
                    SystemBounds::Bounds { first_equation_pivot_id: first_pivot_id, last_present_equation_pivot_id: last_present, largest_nonzero_id: largest_nonzero }
                }
                _ => {
                    SystemBounds::EmptyBounds
                }
            };
            self.n_equations = n_non_null_equations;
        }


        self.recompute_bounds();
    }

    /// removes all the equations whose coef for the given id is not zero
    /// this is useful if we decide to not recover a specific ID, then we 
    /// just drop the equations that need this id to be solved
    pub fn drop_id(&mut self, id: SymbolID) {
        let mut need_to_recompute_bounds = false;
        if let SystemBounds::Bounds { first_equation_pivot_id, last_present_equation_pivot_id, largest_nonzero_id } = self.bounds() {
            let (first_equation_pivot_id, last_present_equation_pivot_id, largest_nonzero_id) = (*first_equation_pivot_id, *last_present_equation_pivot_id, *largest_nonzero_id);
            if id < first_equation_pivot_id {
                return;
            }
            if id > largest_nonzero_id {
                return;
            }
            let mut n_non_null_equations = 0;
            let n_equations = self.n_equations;
            for eq_opt in &mut self.equations {
                if let Some(eq) = eq_opt {
                    n_non_null_equations += 1;
                    if eq.get_coef(id) != 0 {
                        if let EquationBounds::Bounds { pivot, last_nonzero_id } = eq.bounds() {
                            if *pivot == last_present_equation_pivot_id || *last_nonzero_id == largest_nonzero_id {
                                need_to_recompute_bounds = true;
                            }
                        }
                        // drop the equation
                        eq_opt.take();
                        self.n_equations -= 1;
                    }
                }
                if n_non_null_equations >= n_equations {
                    if need_to_recompute_bounds {
                        self._recompute_bounds_heavy();
                    } else {
                        self.recompute_bounds();
                    }
                    return;
                }
            }
        }
    }

    pub fn clear(&mut self) {
        for eq_opt in self.equations.iter_mut() {
            eq_opt.take();
        } 
        self.n_equations = 0;
        self.bounds = SystemBounds::EmptyBounds;
    }

    pub fn print(&self) {
        let mut n_non_null_equations = 0;
        println!("System {:?}", self.bounds);
        if let SystemBounds::Bounds {
            first_equation_pivot_id, last_present_equation_pivot_id: _last_present_equation_pivot_id, largest_nonzero_id
        } = self.bounds {
            for eq in &self.equations {
                if let Some(eq) = eq {
                    n_non_null_equations += 1;
                    for id in first_equation_pivot_id..=largest_nonzero_id {
                        print!("{} ", eq.get_coef(id));
                    }
                    println!();
                }
                if n_non_null_equations == self.n_equations {
                    break;
                }
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use crate::system::System;
    use rand::distributions::{Distribution, Standard};
    use rand::{Rng, SeedableRng};
    use rand_pcg::{Lcg64Xsh32, Pcg32};
    use crate::symbols::{Symbol, SymbolID};
    use crate::system::equation::{eq_clone, Equation};

    fn get_new_rng() -> Lcg64Xsh32 {
        Pcg32::seed_from_u64(42)
    }

    fn get_random_vec(rng: &mut Lcg64Xsh32, size: usize) -> Vec<u8> {
        Standard.sample_iter(rng).take(size).collect()
    }

    fn _generate_equation_from_payloads(rng: &mut Lcg64Xsh32, first_coef_id: SymbolID, coefs: Vec<u8>, symbols: &Vec<Symbol>, symbol_size: usize) -> Equation {
        let mut out_symbol = Symbol::new(first_coef_id, coefs.len() as u64, vec![0; symbol_size]);
        for (i, symbol) in symbols.iter().enumerate() {
            out_symbol.add_mul(coefs[i], &symbols[i]);
        }
        Equation::new(coefs, out_symbol)
    }

    fn generate_equation_from_payloads(rng: &mut Lcg64Xsh32, first_coef_id: SymbolID, n_coefs: u64, symbols: &Vec<Symbol>, symbol_size: usize) -> Equation {
        let coefs = get_random_vec(rng, n_coefs as usize);
        _generate_equation_from_payloads(rng, first_coef_id, coefs, symbols, symbol_size)
    }

    fn generate_equation_from_payloads_with_zeroes(rng: &mut Lcg64Xsh32, first_coef_id: SymbolID, n_coefs: u64, symbols: &Vec<Symbol>, symbol_size: usize, zeroes_rate: f64, already_solved_rate: f64) -> Equation {
        let mut coefs = get_random_vec(rng, n_coefs as usize);
        let already_solved = rng.gen_bool(already_solved_rate);
        let alone_index = rng.gen_range(0..n_coefs) as usize;
        if already_solved {
        }
        for (i, coef) in coefs.iter_mut().enumerate() {
            if (already_solved && i != alone_index) || (!already_solved && rng.gen_bool(zeroes_rate)) {
                *coef = 0;
            }
        }
        _generate_equation_from_payloads(rng, first_coef_id, coefs, symbols, symbol_size)
    }

    #[test]
    fn add_solved_equation() {
        let mut rng = get_new_rng();
        let mut system = System::new();
        let n_protected_symbols = 10;
        let first_id = 40 as SymbolID;
        let lost_id = 44 as SymbolID;
        let mut coefs = vec![0 as u8; n_protected_symbols];
        coefs[(lost_id - first_id) as usize] = 42;
        let data = get_random_vec(&mut rng, 1500);
        let mut solved_eq = Equation::new(coefs,
                                      Symbol::new(first_id,
                                                  n_protected_symbols as u64,
                                                  data));
        let mut solved_eq_copy = eq_clone(&solved_eq);
        solved_eq_copy.normalize_pivot();

        let result = system.add(solved_eq);

        assert!(matches!(result, Ok((None, _))));

        if let Ok((None, vec)) = result {
            assert_eq!(vec.len(), 1);
            assert_eq!(vec[0], lost_id);
            let solved = system.take(lost_id).unwrap();
            assert_eq!(solved, solved_eq_copy.constant_term_data());

        }

    }

    #[test]
    fn add_several_full_equations() {

        let mut rng = get_new_rng();
        let mut system = System::new();
        let first_id = 40 as SymbolID;
        let n_symbols = 50;
        let symbol_size = 1500;
        let mut symbols = Vec::new();
        for i in 0..n_symbols {
            symbols.push(Symbol::new(first_id, n_symbols, get_random_vec(&mut rng, symbol_size)));
        }
        let mut decoded_ids: Vec<SymbolID> = Vec::new();
        for i in 0..n_symbols {
            let result = system.add(generate_equation_from_payloads(&mut rng, first_id, n_symbols, &symbols, symbol_size));
            let (removed, mut decoded) = result.unwrap();
            decoded_ids.append(&mut decoded);
        }

        assert_eq!(decoded_ids.len(), n_symbols as usize);

        for id in decoded_ids {
            let solved = system.take(id).unwrap();
            assert_eq!(&solved, symbols[id as usize - first_id as usize].get_data());
        }
    }

    fn generate_symbols_for_ranges(rng: &mut Lcg64Xsh32, first_coef_id: SymbolID, n_coefs: u32, symbols: &Vec<Symbol>, symbol_size: usize, zeroes_rate: f64, already_solved_rate: f64) {

    }

    #[test]
    fn add_several_full_equations_with_zeroes_in_coefs() {

        let mut rng = get_new_rng();
        let mut system = System::new();
        let first_id = 40 as SymbolID;
        let n_symbols = 50;
        let symbol_size = 1500;
        let system_oversize_ratio = 1.1; // we add 20% more equations than symbols as in this tests we have many coefs set to 0
        let mut symbols = Vec::new();
        let mut n_unused = 0;
        for i in 0..n_symbols {
            symbols.push(Symbol::new(first_id, n_symbols, get_random_vec(&mut rng, symbol_size)));
        }
        let mut decoded_ids: Vec<SymbolID> = Vec::new();
        for i in 0..(((n_symbols as f64)*system_oversize_ratio) as u32) {
            let result = system.add(generate_equation_from_payloads_with_zeroes(&mut rng, first_id, n_symbols, &symbols, symbol_size, 0.3, 0.1));
            match result {
                Ok((removed, mut decoded)) => {
                    decoded_ids.append(&mut decoded);
                }
                Err(crate::system::SystemError::UnusedEquation) => {
                    n_unused += 1;
                }
                Err(_) => panic!("unexpected error"),
            }
        }

        assert_eq!(decoded_ids.len(), n_symbols as usize);

        for id in decoded_ids {
            let solved = system.take(id).unwrap();
            assert_eq!(&solved, symbols[id as usize - first_id as usize].get_data());
        }
    }
}
