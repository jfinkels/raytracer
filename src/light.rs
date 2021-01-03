use crate::hittable::AttenuatedRay;
use crate::hittable::HitRecord;
use crate::hittable::Material;
use crate::ray::Ray;
use crate::texture::Texture;
use crate::vector::Vec3;
use std::rc::Rc;

pub struct DiffuseLight {
    texture: Rc<dyn Texture>,
}

impl DiffuseLight {
    pub fn new(texture: Rc<dyn Texture>) -> DiffuseLight {
        DiffuseLight { texture }
    }
}

impl Material for DiffuseLight {
    fn emit(&self, point: &Vec3, surface_coords: (f64, f64)) -> Vec3 {
        self.texture.value(point, surface_coords)
    }
    fn scatter(&self, ray: Ray, hit_record: HitRecord) -> Option<AttenuatedRay> {
        None
    }
}
