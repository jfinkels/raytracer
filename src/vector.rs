use rand::distributions::Uniform;
use rand::distributions::Distribution;

use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;


#[derive(Debug, Copy, Clone)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {

    pub fn new(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3 { x, y, z }
    }

    pub fn cross(self, other: Vec3) -> Vec3 {
        Vec3::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    pub fn dot(self, other: Vec3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn entrywise_mul(self, other: Vec3) -> Vec3 {
        Vec3::new(
            self.x * other.x,
            self.y * other.y,
            self.z * other.z,
        )
    }

    pub fn interpolate(u: Vec3, v: Vec3, t: f64) -> Vec3 {
        u * t + v * (1. - t)
    }

    pub fn near_zero(&self) -> bool {
        const EPS: f64 = 1e-8;
        return self.x.abs() < EPS && self.y.abs() < EPS && self.z.abs() < EPS;
    }

    pub fn normsquared(self) -> f64 {
        self.dot(self)
    }

    pub fn norm(self) -> f64 {
        self.normsquared().sqrt()
    }

    pub fn unit(self: Vec3) -> Vec3 {
        self / self.norm()
    }

    pub fn zero() -> Vec3 {
        Vec3::new(0., 0., 0.)
    }
}


impl Div<f64> for Vec3 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}


impl Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}


impl Sub for Vec3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}


impl Add for Vec3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}


impl AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}


impl Neg for Vec3 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}


fn random_vec3(min: f64, max: f64) -> Vec3 {
    let unif = Uniform::from(min..max);
    let mut rng = rand::thread_rng();
    let x = unif.sample(&mut rng);
    let y = unif.sample(&mut rng);
    let z = unif.sample(&mut rng);
    Vec3::new(x, y, z)
}


pub fn random_unit_vector() -> Vec3 {
    random_in_unit_sphere().unit()
}


pub fn random_in_unit_sphere() -> Vec3 {
    loop {
        let p = random_vec3(-1., 1.);
        if p.normsquared() < 1. {
            return p;
        }
    }
}


pub fn random_in_unit_disk() -> Vec3 {
    let min = -1.;
    let max = 1.;
    let unif = Uniform::from(min..max);
    let mut rng = rand::thread_rng();
    loop {
        let x = unif.sample(&mut rng);
        let y = unif.sample(&mut rng);
        let z = 0.;
        let p = Vec3::new(x, y, z);
        if p.normsquared() < 1. {
            return p;
        }
    }
}

pub fn random_in_hemisphere(normal: Vec3) -> Vec3 {
    let v = random_in_unit_sphere();
    if v.dot(normal) > 0. {
        v
    } else {
        -v
    }
}
