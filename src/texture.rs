use crate::vector::Vec3;
use std::rc::Rc;

pub trait Texture {
    fn value(&self, point: &Vec3, surface_coords: (f64, f64)) -> Vec3;
}

pub struct SolidColor {
    color: Vec3,
}

impl SolidColor {
    pub fn new(color: Vec3) -> SolidColor {
        SolidColor { color }
    }
}

impl Texture for SolidColor {
    fn value(&self, _point: &Vec3, _surface_coords: (f64, f64)) -> Vec3 {
        self.color
    }
}

pub struct Checker {
    even: Rc<dyn Texture>,
    odd: Rc<dyn Texture>,
}

impl Checker {
    pub fn new(even: Rc<dyn Texture>, odd: Rc<dyn Texture>) -> Checker {
        Checker { even, odd }
    }
}

impl Texture for Checker {
    fn value(&self, point: &Vec3, surface_coords: (f64, f64)) -> Vec3 {
        let f = 10.;
        let Vec3 { x, y, z } = point;
        if (f * x).sin() * (f * y).sin() * (f * z).sin() < 0. {
            self.odd.value(point, surface_coords)
        } else {
            self.even.value(point, surface_coords)
        }
    }
}
