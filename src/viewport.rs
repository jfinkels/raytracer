pub struct Viewport {
    pub width: f64,
    pub height: f64,
}

impl Viewport {
    pub fn new(vfov: f64, aspect_ratio: f64) -> Viewport {
        let height = 2. * (vfov.to_radians() / 2.).tan();
        let width = aspect_ratio * height;
        Viewport {
            width,
            height,
        }
    }
}
