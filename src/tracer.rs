use crate::hittable::Hittable;
use crate::ray::Ray;
use crate::vector::Vec3;

pub struct Tracer {
    world: Vec<Box<dyn Hittable>>,
    max_depth: u8,
}

impl Tracer {
    pub fn new(world: Vec<Box<dyn Hittable>>, max_depth: u8) -> Tracer {
        Tracer { world, max_depth }
    }

    fn background_color(ray: Ray) -> Vec3 {
        let u = ray.direction.unit();
        let t = 0.5 * (u.y + 1.);
        let a = Vec3::new(0.5, 0.7, 1.0);
        let b = Vec3::new(1., 1., 1.);
        Vec3::interpolate(a, b, t)
    }

    fn ray_color(&self, ray: Ray, depth: u8) -> Vec3 {
        if depth <= 0 {
            return Vec3::zero();
        }
        let time_bounds = (0.001, f64::INFINITY);
        let maybe_hit_record = self.world.hits(&ray, time_bounds);
        match maybe_hit_record {
            Option::None => Tracer::background_color(ray),
            Option::Some(hit_record) => {
                // FIXME I don't understand why the clone() is necessary.
                let material = hit_record.material.clone();
                match material.scatter(ray, hit_record) {
                    Some((scattered_ray, attenuation)) => {
                        let recursive_color = self.ray_color(scattered_ray, depth - 1);
                        attenuation.entrywise_mul(recursive_color)
                    }
                    None => Vec3::zero(),
                }
            }
        }
    }

    pub fn trace(&self, ray: Ray) -> Vec3 {
        self.ray_color(ray, self.max_depth)
    }
}
