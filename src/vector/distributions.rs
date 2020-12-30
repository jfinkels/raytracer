use crate::vector::vec3::Vec3;
use rand::distributions::Distribution;
use rand::Rng;
use rand::distributions::Uniform;

struct EntrywiseRandomVec3 {
    min: f64,
    max: f64,
}

impl EntrywiseRandomVec3 {
    fn new(min: f64, max: f64) -> EntrywiseRandomVec3 {
        EntrywiseRandomVec3 {
            min,
            max,
        }
    }
}


impl Distribution<Vec3> for EntrywiseRandomVec3 {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
        let min = self.min;
        let max = self.max;
        let unif = Uniform::from(min..max);
        let x = unif.sample(rng);
        let y = unif.sample(rng);
        let z = unif.sample(rng);
        Vec3::new(x, y, z)
    }
}


struct RandomInUnitCubeVec3 {
    distribution: EntrywiseRandomVec3,
}

impl RandomInUnitCubeVec3 {
    fn new() -> RandomInUnitCubeVec3 {
        let distribution = EntrywiseRandomVec3::new(-1., 1.);
        RandomInUnitCubeVec3 {
            distribution,
        }
    }
}


impl Distribution<Vec3> for RandomInUnitCubeVec3 {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
        self.distribution.sample(rng)
    }
}


pub struct RandomInUnitBallVec3 {
    distribution: RandomInUnitCubeVec3,
}

impl RandomInUnitBallVec3 {
    pub fn new() -> RandomInUnitBallVec3 {
        let distribution = RandomInUnitCubeVec3::new();
        RandomInUnitBallVec3 {
            distribution,
        }
    }
}


impl Distribution<Vec3> for RandomInUnitBallVec3 {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
        loop {
            let p = self.distribution.sample(rng);
            if p.normsquared() < 1. {
                return p;
            }
        }
    }
}


pub struct RandomUnitVec3 {
    distribution: RandomInUnitBallVec3,
}

impl RandomUnitVec3 {
    pub fn new() -> RandomUnitVec3 {
        let distribution = RandomInUnitBallVec3::new();
        RandomUnitVec3 {
            distribution,
        }
    }
}


impl Distribution<Vec3> for RandomUnitVec3 {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
        self.distribution.sample(rng).unit()
    }
}


struct EntrywiseRandomVec2 {
    min: f64,
    max: f64,
}

impl EntrywiseRandomVec2 {
    fn new(min: f64, max: f64) -> EntrywiseRandomVec2 {
        EntrywiseRandomVec2 {
            min,
            max,
        }
    }
}


impl Distribution<(f64, f64)> for EntrywiseRandomVec2 {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> (f64, f64) {
        let min = self.min;
        let max = self.max;
        let unif = Uniform::from(min..max);
        let x = unif.sample(rng);
        let y = unif.sample(rng);
        (x, y)
    }
}


struct RandomInUnitSquareVec2 {
    distribution: EntrywiseRandomVec2,
}

impl RandomInUnitSquareVec2 {
    fn new() -> RandomInUnitSquareVec2 {
        let distribution = EntrywiseRandomVec2::new(-1., 1.);
        RandomInUnitSquareVec2 {
            distribution,
        }
    }
}


impl Distribution<(f64, f64)> for RandomInUnitSquareVec2 {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> (f64, f64) {
        self.distribution.sample(rng)
    }
}



struct RandomInUnitDiskVec2 {
    distribution: RandomInUnitSquareVec2,
}

impl RandomInUnitDiskVec2 {
    fn new() -> RandomInUnitDiskVec2 {
        let distribution = RandomInUnitSquareVec2::new();
        RandomInUnitDiskVec2 {
            distribution,
        }
    }
}


impl Distribution<(f64, f64)> for RandomInUnitDiskVec2 {
    // TODO Refactor this with the nearly identical code from `impl
    // Distribution for RandomInUnitSphereVec3`.
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> (f64, f64) {

        fn normsquared((x, y): (f64, f64)) -> f64 {
            x * x + y * y
        }

        loop {
            let p = self.distribution.sample(rng);
            if normsquared(p) < 1. {
                return p;
            }
        }
    }
}


pub struct RandomInUnitDiskVec3 {
    distribution: RandomInUnitDiskVec2,
}

impl RandomInUnitDiskVec3 {
    pub fn new() -> RandomInUnitDiskVec3 {
        let distribution = RandomInUnitDiskVec2::new();
        RandomInUnitDiskVec3 {
            distribution,
        }
    }
}


impl Distribution<Vec3> for RandomInUnitDiskVec3 {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
        let (x, y) = self.distribution.sample(rng);
        Vec3::new(x, y, 0.)
    }
}
