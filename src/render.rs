use crate::camera::Camera;
use crate::color::Color;
use crate::image::Image;
use crate::tracer::Tracer;
use crate::vector::Vec3;

pub trait PixelRenderer {
    fn render(&self, image: &Image, j: usize, i: usize) -> Vec3;
}

pub struct NoisyPixelRenderer {
    camera: Camera,
    tracer: Tracer,
}

impl NoisyPixelRenderer {
    pub fn new(camera: Camera, tracer: Tracer) -> NoisyPixelRenderer {
        NoisyPixelRenderer { camera, tracer }
    }
}

impl PixelRenderer for NoisyPixelRenderer {
    fn render(&self, image: &Image, j: usize, i: usize) -> Vec3 {
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
        self.tracer.trace(ray)
    }
}

pub struct AveragingPixelRenderer {
    subrenderer: NoisyPixelRenderer,
    samples_per_pixel: usize,
}

impl AveragingPixelRenderer {
    pub fn new(camera: Camera, tracer: Tracer, samples_per_pixel: usize) -> AveragingPixelRenderer {
        let subrenderer = NoisyPixelRenderer::new(camera, tracer);
        AveragingPixelRenderer {
            subrenderer,
            samples_per_pixel,
        }
    }
}

impl PixelRenderer for AveragingPixelRenderer {
    fn render(&self, image: &Image, j: usize, i: usize) -> Vec3 {
        // We are supersampling at each pixel and taking the
        // average color in order to anti-alias.
        let mut sub_color_sum = Vec3::zero();
        for _ in 0..self.samples_per_pixel {
            sub_color_sum += self.subrenderer.render(image, j, i);
        }
        sub_color_sum / (self.samples_per_pixel as f64)
    }
}

pub struct Renderer {
    pixel_renderer: Box<dyn PixelRenderer>,
}

impl Renderer {
    pub fn new(pixel_renderer: Box<dyn PixelRenderer>) -> Renderer {
        Renderer { pixel_renderer }
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
                let color_as_vec3 = self.pixel_renderer.render(&image, j, i);
                let color = Color::from_vec3(color_as_vec3);
                let item = (j, i, color);
                result.push(item);
            }
        }
        result
    }
}
