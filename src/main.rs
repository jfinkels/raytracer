use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;
use rand::distributions::Uniform;
use rand::distributions::Distribution;

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

        // Gamma correction with gamma = 2.
        let red = red.sqrt();
        let green = green.sqrt();
        let blue = blue.sqrt();

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

    fn entrywise_mul(self, other: Vec3) -> Vec3 {
        Vec3::new(
            self.x * other.x,
            self.y * other.y,
            self.z * other.z,
        )
    }

    fn interpolate(u: Vec3, v: Vec3, t: f64) -> Vec3 {
        u * t + v * (1. - t)
    }

    fn near_zero(&self) -> bool {
        const EPS: f64 = 1e-8;
        return self.x.abs() < EPS && self.y.abs() < EPS && self.z.abs() < EPS;
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
    material: Box<dyn Material>,
}


impl HitRecord {
    fn new(ray: &Ray, time: f64, outward_normal: Vec3, material: Box<dyn Material>) -> HitRecord {
        let point = ray.at(time);
        let front_face = ray.direction.dot(outward_normal) < 0.;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };
        HitRecord { point, normal, time, front_face, material }
    }
}


trait Hittable {
    fn hits(&self, ray: &Ray, time_bounds: (f64, f64)) -> Option<HitRecord>;
}


struct Sphere {
    center: Vec3,
    radius: f64,
    material: Box<dyn Material>,
}

impl Sphere {

    fn new(center: Vec3, radius: f64, material: Box<dyn Material>) -> Sphere {
        Sphere { center, radius, material }
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
            // FIXME I don't understand how to avoid cloning here.
            let material = self.material.clone();
            Some(HitRecord::new(ray, t, normal, material))
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


fn background_color(ray: Ray) -> Vec3 {
    let u = ray.direction.unit();
    let t = 0.5 * (u.y + 1.);
    let a = Vec3::new(0.5, 0.7, 1.0);
    let b = Vec3::new(1., 1., 1.);
    Vec3::interpolate(a, b, t)
}

fn random_vec3(min: f64, max: f64) -> Vec3 {
    let unif = Uniform::from(min..max);
    let mut rng = rand::thread_rng();
    let x = unif.sample(&mut rng);
    let y = unif.sample(&mut rng);
    let z = unif.sample(&mut rng);
    Vec3::new(x, y, z)
}


fn random_unit_vector() -> Vec3 {
    random_in_unit_sphere().unit()
}


fn random_in_unit_sphere() -> Vec3 {
    loop {
        let p = random_vec3(-1., 1.);
        if p.normsquared() < 1. {
            return p;
        }
    }
}

fn random_in_hemisphere(normal: Vec3) -> Vec3 {
    let v = random_in_unit_sphere();
    if v.dot(normal) > 0. {
        v
    } else {
        -v
    }
}



fn ray_color(world: &Vec<Box<dyn Hittable>>, ray: Ray, depth: u8) -> Vec3 {
    if depth <= 0 {
        return Vec3::zero();
    }
    let time_bounds = (0.001, f64::INFINITY);
    let maybe_hit_record = world.hits(&ray, time_bounds);
    match maybe_hit_record {
        Option::None => background_color(ray),
        Option::Some(hit_record) => {

            // FIXME I don't understand why the clone() is necessary.
            let material = hit_record.material.clone();
            match material.scatter(ray, hit_record) {
                Some((scattered_ray, attenuation)) => {
                    let recursive_color = ray_color(world, scattered_ray, depth - 1);
                    attenuation.entrywise_mul(recursive_color)
                },
                None => Vec3::zero(),
            }

            // let point = hit_record.point;
            // let normal = hit_record.normal;

            // let origin = point;
            // let target = point + normal + random_in_hemisphere(normal);
            // let direction = target - point;

            // ray_color(world, Ray::new(origin, direction), depth - 1) * 0.5
        }
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



type AttenuatedRay = (Ray, Vec3);


// All this is to allow for cloning `Box<dyn Material>`. See
// https://stackoverflow.com/questions/30353462/how-to-clone-a-struct-storing-a-boxed-trait-object

trait Material : MaterialClone {
    fn scatter(&self, ray: Ray, hit_record: HitRecord) -> Option<AttenuatedRay>;
}


trait MaterialClone {
    fn clone_box(&self) -> Box<dyn Material>;
}

impl<T: 'static + Material + Clone> MaterialClone for T {
    fn clone_box(&self) -> Box<dyn Material> {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn Material> {
    fn clone(&self) -> Box<dyn Material> {
        self.clone_box()
    }
}

#[derive(Clone)]
struct Lambertian {
    albedo: Vec3,
}

impl Lambertian {
    fn new(albedo: Vec3) -> Lambertian {
        Lambertian { albedo }
    }
}

impl Material for Lambertian {
    fn scatter(&self, _ray: Ray, hit_record: HitRecord) -> Option<AttenuatedRay> {
        let origin = hit_record.point;
        let direction = hit_record.normal + random_unit_vector();
        let direction = if direction.near_zero() {
            hit_record.normal
        } else {
            direction
        };
        Some((Ray::new(origin, direction), self.albedo))
    }
}

#[derive(Clone)]
struct Metal {
    albedo: Vec3,
}

impl Metal {
    fn new(albedo: Vec3) -> Metal {
        Metal { albedo }
    }
}

fn reflect(v: Vec3, normal: Vec3) -> Vec3 {
    v - normal * (v.dot(normal) * 2.)
}

impl Material for Metal {
    fn scatter(&self, ray: Ray, hit_record: HitRecord) -> Option<AttenuatedRay> {
        let origin = hit_record.point;
        let direction = reflect(ray.direction.unit(), hit_record.normal);
        let scattered_ray = Ray::new(origin, direction);
        let attenuation = self.albedo;
        Some((scattered_ray, attenuation))
    }
}


fn main() {

    // Image
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: usize = 400;
    const IMAGE_HEIGHT: usize = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as usize;
    const SAMPLES_PER_PIXEL: u8 = 100;
    const MAX_DEPTH: u8 = 10;

    // World

    let center1 = Vec3::new(0., -100.5, -1.);
    let radius1 = 100.;
    let albedo1 = Vec3::new(0.8, 0.8, 0.);
    let material1 = Lambertian::new(albedo1);
    let sphere1 = Sphere::new(center1, radius1, Box::new(material1));

    let center2 = Vec3::new(0., 0., -1.);
    let radius2 = 0.5;
    let albedo2 = Vec3::new(0.7, 0.3, 0.3);
    let material2 = Lambertian::new(albedo2);
    let sphere2 = Sphere::new(center2, radius2, Box::new(material2));

    let center3 = Vec3::new(-1., 0., -1.);
    let radius3 = 0.5;
    let albedo3 = Vec3::new(0.8, 0.8, 0.8);
    let material3 = Metal::new(albedo3);
    let sphere3 = Sphere::new(center3, radius3, Box::new(material3));

    let center4 = Vec3::new(1., 0., -1.);
    let radius4 = 0.5;
    let albedo4 = Vec3::new(0.8, 0.6, 0.2);
    let material4 = Metal::new(albedo4);
    let sphere4 = Sphere::new(center4, radius4, Box::new(material4));

    let world: Vec<Box<dyn Hittable>> = vec![
        Box::new(sphere1),
        Box::new(sphere2),
        Box::new(sphere3),
        Box::new(sphere4),
    ];

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
            eprintln!("tracing pixel {}, {}", i, j);
            // We are supersampling at each pixel and taking the
            // average color in order to anti-alias.
            let mut sub_color_sum = Vec3::zero();
            for _ in 0..SAMPLES_PER_PIXEL {

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

                sub_color_sum += ray_color(&world, ray, MAX_DEPTH);
            }
            let color = Color::from_vec3(sub_color_sum / SAMPLES_PER_PIXEL as f64);
            println!("{}", color.to_ppm());
        }
    }

}
