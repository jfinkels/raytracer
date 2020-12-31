use raytracer::Camera;
use raytracer::Dielectric;
use raytracer::Hittable;
use raytracer::Image;
use raytracer::Lambertian;
use raytracer::Lens;
use raytracer::Metal;
use raytracer::Orientation;
use raytracer::Renderer;
use raytracer::Sphere;
use raytracer::Vec3;
use raytracer::Viewport;

const ASPECT_RATIO: f64 = 16.0 / 9.0;

fn make_image() -> Image {
    const IMAGE_WIDTH: usize = 400;
    const IMAGE_HEIGHT: usize = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as usize;
    Image::new(IMAGE_WIDTH, IMAGE_HEIGHT)
}

fn make_world() -> Vec<Box<dyn Hittable>> {
    let material_ground = Box::new(Lambertian::new(Vec3::new(0.8, 0.8, 0.0)));
    let material_center = Box::new(Lambertian::new(Vec3::new(0.1, 0.2, 0.5)));
    let material_left = Box::new(Dielectric::new(1.5));
    let material_right = Box::new(Metal::new(Vec3::new(0.8, 0.6, 0.2), 0.));

    vec![
        Box::new(Sphere::new(
            Vec3::new(0., -100.5, -1.),
            100.0,
            material_ground,
        )),
        Box::new(Sphere::new(Vec3::new(0., 0., -1.), 0.5, material_center)),
        Box::new(Sphere::new(
            Vec3::new(-1., 0., -1.),
            0.5,
            material_left.clone(),
        )),
        Box::new(Sphere::new(Vec3::new(-1., 0., -1.), -0.45, material_left)),
        Box::new(Sphere::new(Vec3::new(1., 0., -1.), 0.5, material_right)),
    ]
}

fn make_camera() -> Camera {

    // Specify the orientation of the camera.
    let lookfrom = Vec3::new(3., 3., 2.);
    let lookat = Vec3::new(0., 0., -1.);
    let vup = Vec3::new(0., 1., 0.);
    let orientation = Orientation::new(lookfrom, lookat, vup);

    // Create the viewport.
    const VFOV: f64 = 20.;
    let viewport = Viewport::new(VFOV, ASPECT_RATIO);

    // Create the lens.
    const APERTURE: f64 = 0.25;
    let lens = Lens::new(APERTURE);

    let dist_to_focus = (lookfrom - lookat).norm();

    Camera::new(
        orientation,
        viewport,
        lens,
        dist_to_focus,
    )
}

fn make_renderer(camera: Camera, world: Vec<Box<dyn Hittable>>) -> Renderer {
    const SAMPLES_PER_PIXEL: u8 = 50;
    const MAX_DEPTH: u8 = 10;
    Renderer::new(camera, world, SAMPLES_PER_PIXEL, MAX_DEPTH)
}

fn main() {
    let image = make_image();
    let world = make_world();
    let camera = make_camera();
    let renderer = make_renderer(camera, world);

    println!("P3");
    println!("{} {}", image.width, image.height);
    println!("255");
    for (_, _, color) in renderer.render(image) {
        println!("{}", color.to_ppm())
    }
}
