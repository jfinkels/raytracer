use crate::ray::Ray;
use crate::vector::random_in_unit_disk;
use crate::vector::Vec3;

pub struct Camera {
    origin: Vec3,
    lower_left_corner: Vec3,
    horizontal: Vec3,
    vertical: Vec3,
    u: Vec3,
    v: Vec3,
    w: Vec3,
    lens_radius: f64,
}

impl Camera {
    // vertical field of view in degrees
    pub fn new(
        lookfrom: Vec3,
        lookat: Vec3,
        vup: Vec3,
        vfov: f64,
        aspect_ratio: f64,
        aperture: f64,
        focus_dist: f64,
    ) -> Camera {
        let theta = vfov.to_radians();
        let h = (theta / 2.).tan();
        let viewport_height: f64 = 2. * h;
        let viewport_width: f64 = aspect_ratio * viewport_height;

        let w = (lookfrom - lookat).unit();
        let u = vup.cross(w).unit();
        let v = w.cross(u);

        let origin = lookfrom;
        let horizontal = u * viewport_width * focus_dist;
        let vertical = v * viewport_height * focus_dist;

        let lower_left_corner: Vec3 =
            origin - (horizontal / 2.) - (vertical / 2.) - (w * focus_dist);
        let lens_radius = aperture / 2.;

        Camera {
            origin,
            lower_left_corner,
            horizontal,
            vertical,
            u,
            v,
            w,
            lens_radius,
        }
    }

    // `s` and `t` are numbers between 0 and 1, representing how far
    // along the viewport axes to generate the ray.
    pub fn ray_through(&self, s: f64, t: f64) -> Ray {
        let rd = random_in_unit_disk() * self.lens_radius;
        let offset = self.u * rd.x + self.v * rd.y;
        let origin = self.origin + offset;
        let direction =
            self.lower_left_corner + self.horizontal * s + self.vertical * t - self.origin - offset;
        Ray::new(origin, direction)
    }
}
