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

    fn cross(self, other: Vec3) -> Vec3 {
        Vec3::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
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


fn random_in_unit_disk() -> Vec3 {
    let min = -1.;
    let max = 1.;
    let unif = Uniform::from(min..max);
    let mut rng = rand::thread_rng();
    loop {
        let x = unif.sample(&mut rng);
        let y = unif.sample(&mut rng);
        let z = 0.;
        let p = Vec3::new(x, y, z);
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


struct Camera {
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
    fn new(lookfrom: Vec3, lookat: Vec3, vup: Vec3, vfov: f64, aspect_ratio: f64, aperture: f64, focus_dist: f64) -> Camera {
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

        let lower_left_corner: Vec3 = origin - (horizontal / 2.) - (vertical / 2.) - (w * focus_dist);
        let lens_radius = aperture / 2.;

        Camera { origin, lower_left_corner, horizontal, vertical, u, v, w, lens_radius }
        
    }

    // `s` and `t` are numbers between 0 and 1, representing how far
    // along the viewport axes to generate the ray.
    fn ray_through(&self, s: f64, t: f64) -> Ray {
        let rd = random_in_unit_disk() * self.lens_radius;
        let offset = self.u * rd.x + self.v * rd.y;
        let origin = self.origin + offset;
        let direction = self.lower_left_corner + self.horizontal * s + self.vertical * t - self.origin - offset;
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
    fuzz: f64,
}

impl Metal {
    fn new(albedo: Vec3, fuzz: f64) -> Metal {
        Metal { albedo, fuzz }
    }
}

fn reflect(v: Vec3, normal: Vec3) -> Vec3 {
    v - normal * (v.dot(normal) * 2.)
}

impl Material for Metal {
    fn scatter(&self, ray: Ray, hit_record: HitRecord) -> Option<AttenuatedRay> {
        let origin = hit_record.point;
        let reflected = reflect(ray.direction.unit(), hit_record.normal);
        let direction = reflected + random_in_unit_sphere() * self.fuzz;
        let scattered_ray = Ray::new(origin, direction);
        let attenuation = self.albedo;
        Some((scattered_ray, attenuation))
    }
}


#[derive(Clone)]
struct Dielectric {
    refraction_index: f64,
}

impl Dielectric {
    fn new(refraction_index: f64) -> Dielectric {
        Dielectric { refraction_index }
    }
}


// Schlick's approximation for reflectance.
fn reflectance(cosine: f64, refraction_index: f64) -> f64 {
    let r0 = ((1. - refraction_index) / (1. + refraction_index)).powf(2.);
    r0 + (1. - r0) * (1. - cosine).powf(5.)
}

fn is_reflectance_large(cos_theta: f64, refraction_ratio: f64) -> bool {
    reflectance(cos_theta, refraction_ratio) > rand::random::<f64>()
}


fn refract(unit_vec: Vec3, normal: Vec3, eta_i_over_eta_t: f64) -> Vec3 {
    let cos_theta = (-unit_vec).dot(normal).min(1.);
    let r_out_perp = (unit_vec + normal * cos_theta) * eta_i_over_eta_t;
    let r_out_parallel = -normal * (1. - r_out_perp.normsquared()).abs().sqrt();
    r_out_perp + r_out_parallel
}

impl Material for Dielectric {

    fn scatter(&self, ray: Ray, hit_record: HitRecord) -> Option<AttenuatedRay> {
        // Invert the index of refraction if the ray is entering
        // versus exiting the material.
        let refraction_ratio = if hit_record.front_face {
            self.refraction_index.recip()
        } else {
            self.refraction_index
        };

        let origin = hit_record.point;

        // Determine if angle of incidence admits refraction. If not, reflect.
        let unit_direction = ray.direction.unit();
        let cos_theta = (-unit_direction).dot(hit_record.normal).min(1.);
        let sin_theta = (1. - cos_theta.powf(2.)).sqrt();
        let cannot_refract = refraction_ratio * sin_theta > 1.0;
        let direction = if cannot_refract || is_reflectance_large(cos_theta, refraction_ratio) {
            reflect(unit_direction, hit_record.normal)
        } else {
            refract(unit_direction, hit_record.normal, refraction_ratio)
        };

        let scattered_ray = Ray::new(origin, direction);
        let attenuation = Vec3::new(1., 1., 1.);
        Some((scattered_ray, attenuation))
    }
}

fn main() {

    // Image
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: usize = 400;
    const IMAGE_HEIGHT: usize = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as usize;
    const SAMPLES_PER_PIXEL: u8 = 200;
    const MAX_DEPTH: u8 = 20;

    // World

    let material_ground = Box::new(Lambertian::new(Vec3::new(0.8, 0.8, 0.0)));
    let material_center = Box::new(Lambertian::new(Vec3::new(0.1, 0.2, 0.5)));
    let material_left = Box::new(Dielectric::new(1.5));
    let material_right = Box::new(Metal::new(Vec3::new(0.8, 0.6, 0.2), 0.));

    let world: Vec<Box<dyn Hittable>> = vec![
        Box::new(Sphere::new(Vec3::new(0., -100.5, -1.), 100.0, material_ground)),
        Box::new(Sphere::new(Vec3::new(0., 0., -1.), 0.5, material_center)),
        Box::new(Sphere::new(Vec3::new(-1., 0., -1.), 0.5, material_left.clone())),
        Box::new(Sphere::new(Vec3::new(-1., 0., -1.), -0.45, material_left)),
        Box::new(Sphere::new(Vec3::new(1., 0., -1.), 0.5, material_right)),
    ];

    // Camera
    //
    // - x is positive to the right,
    // - y is positive going up,
    // - z is positive *coming out of the screen*.
    //
    let lookfrom = Vec3::new(3., 3., 2.);
    let lookat = Vec3::new(0., 0., -1.);
    let vup = Vec3::new(0., 1., 0.);
    let dist_to_focus = (lookfrom - lookat).norm();
    const APERTURE: f64 = 1.;
    const VFOV: f64 = 20.;
    let camera = Camera::new(lookfrom, lookat, vup, VFOV, ASPECT_RATIO, APERTURE, dist_to_focus);

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
