use crate::camera::Camera;
use crate::color::Color;
use crate::hittable::Hittable;
use crate::image::Image;
use crate::ray::Ray;
use crate::vector::Vec3;

pub struct Renderer {
    camera: Camera,
    world: Vec<Box<dyn Hittable>>,
    samples_per_pixel: u8,
    max_depth: u8,
}

impl Renderer {
    pub fn new(
        camera: Camera,
        world: Vec<Box<dyn Hittable>>,
        samples_per_pixel: u8,
        max_depth: u8,
    ) -> Renderer {
        Renderer {
            camera,
            world,
            samples_per_pixel,
            max_depth,
        }
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
            Option::None => Renderer::background_color(ray),
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
    pub fn render(&self, image: Image) -> Vec<(usize, usize, Color)> {
        // Iterate over each pixel in the image.
        //
        // Following
        // https://raytracing.github.io/books/RayTracingInOneWeekend.html
        //
        // For some reason, this author decided to use `i` as the column
        // number and `j` as the row number. I think this is so that `i`
        // matches `x` and `j` matches `y`.
        //
        // The number `i` ranges from 0 on the left to `IMAGE_WIDTH` on
        // the right. The number `j` ranges from 0 on the bottom to
        // `IMAGE_HEIGHT` on the top.
        let mut result = Vec::new();
        for j in (0..image.height).rev() {
            for i in 0..image.width {
                eprintln!("tracing pixel {}, {}", i, j);
                // We are supersampling at each pixel and taking the
                // average color in order to anti-alias.
                let mut sub_color_sum = Vec3::zero();
                for _ in 0..self.samples_per_pixel {
                    // Compute the direction of the vector from the camera to
                    // the viewport. The camera is at the origin.
                    //
                    // The number `u` determines how far the pixel is along
                    // the horizontal axis, 0 meaning all the way left and 1
                    // meaning all the way right.
                    //
                    // The number `v` determines how far the pixel is along
                    // the vertical axis, 0 meaning all the way at the bottom
                    // and 1 meaning all the way at the top.
                    let u = (i as f64 + rand::random::<f64>()) / (image.width as f64 - 1.);
                    let v = (j as f64 + rand::random::<f64>()) / (image.height as f64 - 1.);
                    let ray = self.camera.ray_through(u, v);

                    sub_color_sum += self.ray_color(ray, self.max_depth);
                }
                let color = Color::from_vec3(sub_color_sum / self.samples_per_pixel as f64);
                result.push((j, i, color));
            }
        }
        result
    }
}
