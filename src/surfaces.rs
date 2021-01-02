use crate::hittable::HitRecord;
use crate::hittable::Hittable;
use crate::hittable::Material;
use crate::ray::Ray;
use crate::time::Duration;
use crate::vector::Vec3;
use std::rc::Rc;

pub struct Sphere {
    center_start: Vec3,
    center_end: Vec3,
    duration: Duration,
    radius: f64,
    material: Rc<dyn Material>,
}

impl Sphere {
    pub fn new(center: Vec3, radius: f64, material: Rc<dyn Material>) -> Sphere {
        let center_start = center;
        let center_end = center;
        let duration = Duration::max();
        Sphere::new_moving(center_start, center_end, duration, radius, material)
    }

    pub fn new_moving(
        center_start: Vec3,
        center_end: Vec3,
        duration: Duration,
        radius: f64,
        material: Rc<dyn Material>,
    ) -> Sphere {
        Sphere {
            center_start,
            center_end,
            duration,
            radius,
            material,
        }
    }

    fn center(&self, time: f64) -> Vec3 {
        let t = (time - self.duration.start) / self.duration.len();
        let a = self.center_start;
        let b = self.direction();
        a + b * t
    }

    fn direction(&self) -> Vec3 {
        self.center_end - self.center_start
    }
}

impl Hittable for Sphere {
    fn hits(&self, ray: &Ray, time_bounds: (f64, f64)) -> Option<HitRecord> {
        let center = self.center(ray.time);
        let oc = ray.origin - center;
        let a = ray.direction.normsquared();
        let b = ray.direction.dot(oc) * 2.;
        let c = oc.normsquared() - self.radius.powf(2.);
        let discriminant = b * b - 4. * a * c;
        if discriminant < 0. {
            None
        } else {
            let sqrtd = discriminant.sqrt();

            let (t_min, t_max) = time_bounds;

            // There may be two roots, we want the nearest one.
            let root1 = (-b - sqrtd) / (2. * a);
            let root2 = (-b + sqrtd) / (2. * a);

            let t = if t_min <= root1 && root1 <= t_max {
                root1
            } else {
                if t_min <= root2 && root2 <= t_max {
                    root2
                } else {
                    return None;
                }
            };

            let point = ray.at(t);
            let normal = (point - center) / self.radius;
            // Because ownership of a material might be shared among
            // multiple surfaces *and* among multiple `HitRecord`
            // structs, we use a reference counting smart pointer.
            let material = Rc::clone(&self.material);
            Some(HitRecord::new(ray, t, normal, material))
        }
    }
}
