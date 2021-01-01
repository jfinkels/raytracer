use crate::ray::Ray;
use crate::vector::Vec3;
use std::rc::Rc;

pub type AttenuatedRay = (Ray, Vec3);

pub trait Material {
    fn scatter(&self, ray: Ray, hit_record: HitRecord) -> Option<AttenuatedRay>;
}

pub struct HitRecord {
    pub point: Vec3,
    pub normal: Vec3,
    pub time: f64,
    pub front_face: bool,
    pub material: Rc<dyn Material>,
}

impl HitRecord {
    pub fn new(
        ray: &Ray,
        time: f64,
        outward_normal: Vec3,
        material: Rc<dyn Material>,
    ) -> HitRecord {
        let point = ray.at(time);
        let front_face = ray.direction.dot(outward_normal) < 0.;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };
        HitRecord {
            point,
            normal,
            time,
            front_face,
            material,
        }
    }
}

pub trait Hittable {
    fn hits(&self, ray: &Ray, time_bounds: (f64, f64)) -> Option<HitRecord>;
}

// TODO This could be more generic, but I don't know how to do it in
// Rust. We could implement `Hittable` for any iterator over `Box<dyn
// Hittable>` structs, not just `Vec<Box<dyn Hittable>>`.
impl Hittable for Vec<Box<dyn Hittable>> {
    fn hits(&self, ray: &Ray, time_bounds: (f64, f64)) -> Option<HitRecord> {
        let mut maybe_closest = None;
        for hittable in self.into_iter() {
            if let Some(hit_record) = hittable.hits(ray, time_bounds) {
                match maybe_closest {
                    None => {
                        maybe_closest = Some(hit_record);
                    }
                    Some(ref closest) => {
                        if !hit_record.time.is_nan() {
                            if hit_record.time < closest.time {
                                maybe_closest = Some(hit_record);
                            }
                        }
                    }
                };
            };
        }
        maybe_closest
    }
}
