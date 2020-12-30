use rand::distributions::Distribution;
use crate::hittable::AttenuatedRay;
use crate::hittable::HitRecord;
use crate::hittable::Material;
use crate::ray::Ray;
use crate::vector::RandomInUnitBallVec3;
use crate::vector::RandomUnitVec3;
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

        // The incoming ray is ignored, and we treat the surface of
        // the material as emitting its own color.
        let emitted = hit_record.normal;

        // Generate some random noise for the emitted vector.
        let mut rng = rand::thread_rng();
        let distribution = RandomUnitVec3::new();
        let noise = distribution.sample(&mut rng);

        // The scattered ray originates at the hit point and goes in
        // the noisy emission direction.
        //
        // If the noise would put the vector too close to the zero
        // vector, then we just eliminate the noise altogether.
        let origin = hit_record.point;
        let noisy_emitted = emitted + noise;
        let direction = if noisy_emitted.near_zero() {
            emitted
        } else {
            noisy_emitted
        };
        let scattered_ray = Ray::new(origin, direction);


        // The attenuation is the color of the material itself.
        let attenuation = self.albedo;

        Some((scattered_ray, attenuation))
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

        // Reflect the incoming ray off the plane defined by the
        // surface normal.
        let unit_direction = ray.direction.unit();
        let reflected = reflect(unit_direction, hit_record.normal);

        // Generate some random noise for the reflected vector.
        //
        // When the noise factor `self.fuzz` is zero, the metal is a
        // perfect mirror. As the noise factor increases, the metal
        // becomes less shiny and more matte.
        let mut rng = rand::thread_rng();
        let distribution = RandomInUnitBallVec3::new();
        let random_in_unit_ball = distribution.sample(&mut rng);
        let noise = random_in_unit_ball * self.fuzz;

        // The scattered ray originates at the hit point and goes in
        // the noisy reflected direction.
        let origin = hit_record.point;
        let direction = reflected + noise;
        let scattered_ray = Ray::new(origin, direction);

        // The attenuation is the color of the material itself.
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
fn reflectance(cos_theta: f64, refraction_index: f64) -> f64 {
    let r0 = ((1. - refraction_index) / (1. + refraction_index)).powf(2.);
    r0 + (1. - r0) * (1. - cos_theta).powf(5.)
}


fn refract(unit_vec: Vec3, normal: Vec3, refraction_index: f64) -> Vec3 {
    let cos_theta = -unit_vec.dot(normal).min(1.);
    let r_out_perp = (unit_vec + normal * cos_theta) * refraction_index;
    let r_out_parallel = -normal * (1. - r_out_perp.normsquared()).abs().sqrt();
    r_out_perp + r_out_parallel
}


fn should_reflect(unit_vec: Vec3, normal: Vec3, refraction_index: f64) -> bool {
    let cos_theta = -unit_vec.dot(normal).min(1.);
    let sin_theta = (1. - cos_theta.powf(2.)).sqrt();
    let cannot_refract = refraction_index * sin_theta > 1.0;
    let is_reflectance_large = reflectance(cos_theta, refraction_index) > rand::random::<f64>();
    cannot_refract || is_reflectance_large
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

        // Decide whether to reflect or refract the ray.
        //
        // Determine if angle of incidence admits refraction. If not, reflect.
        let unit_direction = ray.direction.unit();
        let normal = hit_record.normal;
        let should_reflect_ = should_reflect(unit_direction, normal, refraction_ratio);
        let new_direction = if should_reflect_ {
            reflect(unit_direction, normal)
        } else {
            refract(unit_direction, normal, refraction_ratio)
        };

        // Generate some random noise for the reflected/refracted vector.
        let noise = Vec3::zero();

        // The scattered ray originates at the hit point and goes in
        // the noisy reflected/refracted direction.
        let origin = hit_record.point;
        let direction = new_direction + noise;
        let scattered_ray = Ray::new(origin, direction);

        let attenuation = Vec3::new(1., 1., 1.);

        Some((scattered_ray, attenuation))
    }
}
