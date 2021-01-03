use raytracer::AveragingPixelRenderer;
use raytracer::Camera;
use raytracer::Checker;
use raytracer::ConstantMedium;
use raytracer::Dielectric;
use raytracer::DiffuseLight;
use raytracer::Duration;
use raytracer::Hittable;
use raytracer::Image;
use raytracer::Lambertian;
use raytracer::Lens;
use raytracer::Metal;
use raytracer::Orientation;
use raytracer::Rectangle;
use raytracer::Renderer;
use raytracer::SolidColor;
use raytracer::Sphere;
use raytracer::Tracer;
use raytracer::Vec3;
use raytracer::Viewport;
use std::rc::Rc;

const ASPECT_RATIO: f64 = 16.0 / 9.0;

fn make_image() -> Image {
    const IMAGE_WIDTH: usize = 400;
    const IMAGE_HEIGHT: usize = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as usize;
    Image::new(IMAGE_WIDTH, IMAGE_HEIGHT)
}

fn make_world() -> Vec<Box<dyn Hittable>> {
    let material_ground = Rc::new(Lambertian::new(Rc::new(Checker::new(
        Rc::new(SolidColor::new(Vec3::new(0.8, 0.8, 0.0))),
        Rc::new(SolidColor::new(Vec3::new(0.5, 0.5, 0.0))),
    ))));
    let material_center = Rc::new(Lambertian::new(Rc::new(SolidColor::new(Vec3::new(
        0.1, 0.2, 0.5,
    )))));
    let material_left = Rc::new(Dielectric::new(1.5));
    let material_right = Rc::new(Metal::new(Vec3::new(0.8, 0.6, 0.2), 0.));

    let material_center2 = Rc::new(Lambertian::new(Rc::new(SolidColor::new(Vec3::new(
        0.4, 0.1, 0.1,
    )))));
    let material_left2 = Rc::new(Dielectric::new(2.5));
    let material_right2 = Rc::new(Metal::new(Vec3::new(0.7, 0.7, 0.7), 0.1));

    let material_light = Rc::new(DiffuseLight::new(Rc::new(SolidColor::new(Vec3::new(
        4., 4., 4.,
    )))));

    vec![
        Box::new(Sphere::new(
            Vec3::new(0., -100.5, -1.),
            100.0,
            material_ground,
        )),
        Box::new(Sphere::new_moving(
            Vec3::new(0., 0., 0.),
            Vec3::new(0., 0., -1.),
            Duration::new(0., 1.),
            0.5,
            material_center.clone(),
        )),
        Box::new(Sphere::new_moving(
            Vec3::new(-1., 0., -1.),
            Vec3::new(-1.25, 0., -1.),
            Duration::new(0., 1.),
            0.5,
            material_left.clone(),
        )),
        Box::new(Sphere::new_moving(
            Vec3::new(-1., 0., -1.),
            Vec3::new(-1.25, 0., -1.),
            Duration::new(0., 1.),
            -0.45,
            material_left.clone(),
        )),
        Box::new(Sphere::new_moving(
            Vec3::new(1., 0., -1.),
            Vec3::new(1.2, 0., -1.),
            Duration::new(0., 1.),
            0.5,
            material_right.clone(),
        )),
        Box::new(Sphere::new_moving(
            Vec3::new(-0.5, 0., -2.),
            Vec3::new(-0.6, 0., -2.2),
            Duration::new(0., 1.),
            0.5,
            material_left2.clone(),
        )),
        Box::new(Sphere::new_moving(
            Vec3::new(-0.5, 0., -2.),
            Vec3::new(-0.6, 0., -2.2),
            Duration::new(0., 1.),
            -0.45,
            material_left2.clone(),
        )),
        Box::new(Sphere::new(Vec3::new(-1.5, 0., -2.), 0.5, material_right2)),
        Box::new(Sphere::new(Vec3::new(0.5, 0., -2.), 0.5, material_center2)),
        Box::new(ConstantMedium::new(
            Box::new(Sphere::new(
                Vec3::new(1., 1., -3.),
                0.5,
                Rc::new(Lambertian::new(Rc::new(SolidColor::new(Vec3::new(
                    1., 1., 1.,
                ))))),
            )),
            Rc::new(SolidColor::new(Vec3::new(1., 1., 1.))),
            0.2,
        )),
        // lights
        Box::new(Rectangle::new(
            (-0.5, 0.),
            (0.5, 1.),
            -4.,
            material_light.clone(),
        )),
        Box::new(Sphere::new(
            Vec3::new(0., 5., -1.),
            3.,
            material_light.clone(),
        )),
    ]
}

fn make_camera() -> Camera {
    // Specify the orientation of the camera.
    let lookfrom = Vec3::new(3., 3., 2.);
    let lookat = Vec3::new(0., 0., -1.);
    let vup = Vec3::new(0., 1., 0.);
    let orientation = Orientation::new(lookfrom, lookat, vup);

    // Create the viewport.
    const VFOV: f64 = 35.;
    let viewport = Viewport::new(VFOV, ASPECT_RATIO);

    // Create the lens.
    const APERTURE: f64 = 0.25;
    let lens = Lens::new(APERTURE);

    // Define the time at which the shutter opens and closes.
    let duration = Duration::new(0., 1.);

    // Define the focus distance.
    let dist_to_focus = (lookfrom - lookat).norm();

    Camera::new(orientation, viewport, lens, dist_to_focus, duration)
}

fn make_renderer(camera: Camera, world: Vec<Box<dyn Hittable>>) -> Renderer {
    const MAX_DEPTH: u8 = 10;
    // let background_color = Vec3::new(0.7, 0.8, 1.);
    let background_color = Vec3::zero();
    let tracer = Tracer::new(world, MAX_DEPTH, background_color);

    const SAMPLES_PER_PIXEL: usize = 100;
    let pixel_renderer = Box::new(AveragingPixelRenderer::new(
        camera,
        tracer,
        SAMPLES_PER_PIXEL,
    ));

    Renderer::new(pixel_renderer)
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
