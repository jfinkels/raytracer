use crate::vector::Vec3;

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
