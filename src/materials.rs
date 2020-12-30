use crate::hittable::AttenuatedRay;
use crate::hittable::HitRecord;
use crate::hittable::Material;
use crate::ray::Ray;
use crate::vector::random_in_unit_sphere;
use crate::vector::random_unit_vector;
use crate::vector::Vec3;

#[derive(Clone)]
pub struct Lambertian {
    albedo: Vec3,
}

impl Lambertian {
    pub fn new(albedo: Vec3) -> Lambertian {
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
pub struct Metal {
    albedo: Vec3,
    fuzz: f64,
}

impl Metal {
    pub fn new(albedo: Vec3, fuzz: f64) -> Metal {
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
pub struct Dielectric {
    refraction_index: f64,
}

impl Dielectric {
    pub fn new(refraction_index: f64) -> Dielectric {
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
