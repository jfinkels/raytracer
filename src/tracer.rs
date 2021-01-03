use crate::hittable::Hittable;
use crate::ray::Ray;
use crate::vector::Vec3;
use std::rc::Rc;

pub struct Tracer {
    world: Vec<Box<dyn Hittable>>,
    max_depth: u8,
    background_color: Vec3,
}

impl Tracer {
    pub fn new(world: Vec<Box<dyn Hittable>>, max_depth: u8, background_color: Vec3) -> Tracer {
        Tracer {
            world,
            max_depth,
            background_color,
        }
    }

    fn ray_color(&self, ray: Ray, depth: u8) -> Vec3 {
        if depth <= 0 {
            return Vec3::zero();
        }
        let time_bounds = (0.001, f64::INFINITY);
        let maybe_hit_record = self.world.hits(&ray, time_bounds);
        match maybe_hit_record {
            Option::None => self.background_color,
            Option::Some(hit_record) => {
                let material = Rc::clone(&hit_record.material);
                let point = hit_record.point;
                let surface_coords = hit_record.surface_coords;
                let emitted = material.emit(&point, surface_coords);
                match material.scatter(ray, hit_record) {
                    Some((scattered_ray, attenuation)) => {
                        let recursive_color = self.ray_color(scattered_ray, depth - 1);
                        emitted + attenuation.entrywise_mul(recursive_color)
                    }
                    None => emitted,
                }
            }
        }
    }

    pub fn trace(&self, ray: Ray) -> Vec3 {
        self.ray_color(ray, self.max_depth)
    }
}
