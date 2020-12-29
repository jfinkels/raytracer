use std::ops::Add;
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

impl Color {

    fn new(red: u8, green: u8, blue: u8) -> Color {
        Color {red, green, blue}
    }

    fn from_ratios(red: f64, green: f64, blue: f64) -> Color {
        Color::new((red * 255.) as u8, (green * 255.) as u8, (blue * 255.) as u8)
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



fn ray_color(world: &Vec<Box<dyn Hittable>>, ray: Ray) -> Color {

    fn background_color(ray: Ray) -> Color {
        let u = ray.direction.unit();
        let t = 0.5 * (u.y + 1.);
        let a = Vec3::new(0.5, 0.7, 1.0);
        let b = Vec3::new(1., 1., 1.);
        let ratios = Vec3::interpolate(a, b, t);

        let red = ratios.x;
        let green = ratios.y;
        let blue = ratios.z;

        Color::from_ratios(red, green, blue)
    }

    let time_bounds = (0., f64::INFINITY);
    let maybe_hit_record = world.hits(&ray, time_bounds);
    match maybe_hit_record {
        Option::None => background_color(ray),
        Option::Some(hit_record) => {

            // Compute the color for this hit, represented as
            // percentages for each color.
            let normal = hit_record.normal;
            let color_as_ratios = (normal + Vec3::new(1., 1., 1.)) * 0.5;

            // Convert the percentages to a `Color` struct.
            let red = color_as_ratios.x;
            let green = color_as_ratios.y;
            let blue = color_as_ratios.z;
            Color::from_ratios(red, green, blue)

        },
    }
}


fn main() {

    // Image
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: usize = 400;
    const IMAGE_HEIGHT: usize = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as usize;

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
    const VIEWPORT_HEIGHT: f64 = 2.;
    const VIEWPORT_WIDTH: f64 = ASPECT_RATIO * VIEWPORT_HEIGHT;
    const FOCAL_LENGTH: f64 = 1.;

    let origin: Vec3 = Vec3::new(0., 0., 0.);
    let horizontal: Vec3 = Vec3::new(VIEWPORT_WIDTH, 0., 0.);
    let vertical: Vec3 = Vec3::new(0., VIEWPORT_HEIGHT, 0.);
    let depth: Vec3 = Vec3::new(0., 0., FOCAL_LENGTH);

    let lower_left_corner: Vec3 = origin - horizontal / 2. - vertical / 2. - depth;

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
            let u = i as f64 / (IMAGE_WIDTH as f64 - 1.);
            let v = j as f64 / (IMAGE_HEIGHT as f64 - 1.);
            let direction = lower_left_corner + horizontal * u + vertical * v - origin;
            let ray = Ray::new(origin, direction);
            let color = ray_color(&world, ray);
            println!("{}", color.to_ppm());
        }
    }

}
