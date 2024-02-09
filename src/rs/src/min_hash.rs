use rand::RngCore;

pub struct MinHashFunction {
    shift: u32,
    mult: u32,
}

impl MinHashFunction {
    pub fn new(rng: &mut rand::rngs::ThreadRng) -> Self {
        Self {
            shift: (rng.next_u64()%(4_294_967_291-1)) as u32,
            mult: (rng.next_u64()%(4_294_967_291-2)) as u32+1,
        }
    }

    pub fn hash(&self, value: u64) -> u32 {
        let hashed = (self.mult as u64*value+self.shift as u64)%4_294_967_291;
        return hashed as u32;
    }
}

#[derive(Debug)]
pub struct MinHashSignature {
    pub bands: Vec<MinHashBandSignature>
}

impl MinHashSignature {
    fn new(num_bands: usize) -> Self {
        Self {
            bands: (0..num_bands).map(|_| MinHashBandSignature::new()).collect()
        }
    }
}

#[derive(Debug)]
pub struct MinHashBandSignature {
    pub signature: u32
}

impl MinHashBandSignature {
    pub fn new() -> Self {
        Self {
            signature: u32::MAX
        }
    }
}

pub struct MinHashBand {
    fnc: MinHashFunction,
}

impl MinHashBand {
    pub fn new(rng: &mut rand::rngs::ThreadRng) -> Self {
        Self {
            fnc: MinHashFunction::new(rng)
        }
    }
}

pub struct MinHash {
    num_bands: usize,
    bands: Vec<MinHashBand>
}

impl MinHash {
    pub fn new(num_bands: usize, rng: &mut rand::rngs::ThreadRng) -> Self {
        let mut bands = Vec::new();
        for _ in 0..num_bands {
            bands.push(MinHashBand::new(rng));
        }
        Self {
            bands,
            num_bands
        }
    }

    pub fn generate_signature(&self, value: u64) -> MinHashSignature {
        let mut sig = MinHashSignature::new(self.num_bands);
        self.add(value,&mut sig);
        sig
    }

    pub fn add(&self, value: u64, sig: &mut MinHashSignature) {
        for i in 0..self.num_bands {
            let band = &self.bands[i];
            let val = band.fnc.hash(value);
            sig.bands[i].signature = sig.bands[i].signature.min(val);
        }
    }
}
