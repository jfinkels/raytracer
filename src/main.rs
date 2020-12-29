use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;

#[derive(Debug)]
struct Color {
    red: u8,
    green: u8,
    blue: u8,
}

fn clamp(x: f64, min: f64, max: f64) -> f64 {
    x.min(max).max(min)
}


impl Color {

    fn new(red: u8, green: u8, blue: u8) -> Color {
        Color {red, green, blue}
    }

    fn from_ratios(red: f64, green: f64, blue: f64) -> Color {
        let red_u8 = (256. * clamp(red, 0., 0.999)) as u8;
        let green_u8 = (256. * clamp(green, 0., 0.999)) as u8;
        let blue_u8 = (256. * clamp(blue, 0., 0.999)) as u8;
        Color::new(red_u8, green_u8, blue_u8)
    }

    fn from_vec3(v: Vec3) -> Color {
        Color::from_ratios(v.x, v.y, v.z)
    }

    fn to_ppm(&self) -> String {
        format!("{} {} {}", self.red, self.green, self.blue)
    }
}


#[derive(Debug, Copy, Clone)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3 {

    fn new(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3 { x, y, z }
    }

    fn dot(self, other: Vec3) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn interpolate(u: Vec3, v: Vec3, t: f64) -> Vec3 {
        u * t + v * (1. - t)
    }

    fn normsquared(self) -> f64 {
        self.dot(self)
    }

    fn norm(self) -> f64 {
        self.normsquared().sqrt()
    }

    fn unit(self: Vec3) -> Vec3 {
        self / self.norm()
    }

    fn zero() -> Vec3 {
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


#[derive(Debug)]
struct Ray {
    origin: Vec3,
    direction: Vec3,
}


impl Ray {

    fn new(origin: Vec3, direction: Vec3) -> Ray {
        Ray { origin, direction }
    }

    fn at(&self, t: f64) -> Vec3 {
        self.origin + self.direction * t
    }

}


struct HitRecord {
    point: Vec3,
    normal: Vec3,
    time: f64,
    front_face: bool,
}


impl HitRecord {
    fn new(ray: &Ray, time: f64, outward_normal: Vec3) -> HitRecord {
        let point = ray.at(time);
        let front_face = ray.direction.dot(outward_normal) < 0.;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };
        HitRecord { point, normal, time, front_face }
    }
}


trait Hittable {
    fn hits(&self, ray: &Ray, time_bounds: (f64, f64)) -> Option<HitRecord>;
}


struct Sphere {
    center: Vec3,
    radius: f64,
}

impl Sphere {

    fn new(center: Vec3, radius: f64) -> Sphere {
        Sphere { center, radius }
    }

}

impl Hittable for Sphere {

    fn hits(&self, ray: &Ray, time_bounds: (f64, f64)) -> Option<HitRecord> {
        let oc = ray.origin - self.center;
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
                    return None
                }
            };

            let point = ray.at(t);
            let normal = (point - self.center) / self.radius;
            Some(HitRecord::new(ray, t, normal))
        }
    }

}


impl Hittable for Vec<Box<dyn Hittable>> {

    fn hits(&self, ray: &Ray, time_bounds: (f64, f64)) -> Option<HitRecord> {
        let mut maybe_closest = None;
        for hittable in self.into_iter() {
            if let Some(hit_record) = hittable.hits(ray, time_bounds) {
                match maybe_closest {
                    None => {
                        maybe_closest = Some(hit_record);
                    },
                    Some(ref closest) => {
                        if !hit_record.time.is_nan() {
                            if hit_record.time < closest.time {
                                maybe_closest = Some(hit_record);
                            }
                        }
                    },
                };
            };
        };
        maybe_closest

    }
}



fn ray_color(world: &Vec<Box<dyn Hittable>>, ray: Ray) -> Vec3 {

    fn background_color(ray: Ray) -> Vec3 {
        let u = ray.direction.unit();
        let t = 0.5 * (u.y + 1.);
        let a = Vec3::new(0.5, 0.7, 1.0);
        let b = Vec3::new(1., 1., 1.);
        Vec3::interpolate(a, b, t)
    }

    let time_bounds = (0., f64::INFINITY);
    let maybe_hit_record = world.hits(&ray, time_bounds);
    match maybe_hit_record {
        Option::None => background_color(ray),
        Option::Some(hit_record) => (hit_record.normal + Vec3::new(1., 1., 1.)) * 0.5,
    }
}


struct Camera {}

impl Camera {

    // `u` and `v` are numbers between 0 and 1, representing how far
    // along the viewport axes to generate the ray.
    fn ray_through(&self, u: f64, v: f64) -> Ray {

        const ASPECT_RATIO: f64 = 16.0 / 9.0;
        const VIEWPORT_HEIGHT: f64 = 2.;
        const VIEWPORT_WIDTH: f64 = ASPECT_RATIO * VIEWPORT_HEIGHT;
        const FOCAL_LENGTH: f64 = 1.;

        let origin: Vec3 = Vec3::new(0., 0., 0.);
        let horizontal: Vec3 = Vec3::new(VIEWPORT_WIDTH, 0., 0.);
        let vertical: Vec3 = Vec3::new(0., VIEWPORT_HEIGHT, 0.);
        let depth: Vec3 = Vec3::new(0., 0., FOCAL_LENGTH);

        let lower_left_corner: Vec3 = origin - horizontal / 2. - vertical / 2. - depth;

        let direction = lower_left_corner + horizontal * u + vertical * v - origin;
        Ray::new(origin, direction)
    }

}



fn main() {

    // Image
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: usize = 400;
    const IMAGE_HEIGHT: usize = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as usize;
    const SAMPLES_PER_PIXEL: u8 = 100;

    // World
    let center1 = Vec3::new(0., 0., -1.);
    let radius1 = 0.5;
    let sphere1 = Sphere::new(center1, radius1);
    let center2 = Vec3::new(0., -100.5, -1.);
    let radius2 = 100.;
    let sphere2 = Sphere::new(center2, radius2);
    let world: Vec<Box<dyn Hittable>> = vec![Box::new(sphere1), Box::new(sphere2)];

    // Camera
    //
    // - x is positive to the right,
    // - y is positive going up,
    // - z is positive *coming out of the screen*.
    //
    let camera = Camera {};

    // Render
    println!("P3");
    println!("{} {}", IMAGE_WIDTH, IMAGE_HEIGHT);
    println!("255");

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
    for j in (0..IMAGE_HEIGHT).rev() {
        for i in 0..IMAGE_WIDTH {

            // We are supersampling at each pixel and taking the
            // average color in order to anti-alias.
            let mut sub_color_sum = Vec3::zero();
            for s in 0..SAMPLES_PER_PIXEL {

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
                let u = (i as f64 + rand::random::<f64>()) / (IMAGE_WIDTH as f64 - 1.);
                let v = (j as f64 + rand::random::<f64>()) / (IMAGE_HEIGHT as f64 - 1.);
                let ray = camera.ray_through(u, v);

                let sub_color = ray_color(&world, ray);
                sub_color_sum += sub_color;
            }
            let color = Color::from_vec3(sub_color_sum / SAMPLES_PER_PIXEL as f64);
            println!("{}", color.to_ppm());
        }
    }

}
