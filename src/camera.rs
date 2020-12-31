use crate::ray::Ray;
use crate::vector::RandomInUnitDiskVec3;
use crate::vector::Vec3;
use crate::viewport::Viewport;
use rand::distributions::Distribution;

pub struct Orientation {
    lookfrom: Vec3,
    lookat: Vec3,
    vup: Vec3,
}

impl Orientation {
    pub fn new(lookfrom: Vec3, lookat: Vec3, vup: Vec3) -> Orientation {
        Orientation {
            lookfrom,
            lookat,
            vup,
        }
    }
}

pub struct Lens {
    radius: f64,
}

impl Lens {
    pub fn new(aperture: f64) -> Lens {
        let radius = aperture / 2.;
        Lens { radius }
    }
}

struct Basis {
    u: Vec3,
    v: Vec3,
    w: Vec3,
}

impl Basis {
    fn new(u: Vec3, v: Vec3, w: Vec3) -> Basis {
        Basis { u, v, w }
    }
}

struct ViewportDimensions {
    horizontal: Vec3,
    vertical: Vec3,
    depth: Vec3,
}

impl ViewportDimensions {
    fn new(basis: &Basis, viewport: Viewport, focus_dist: f64) -> ViewportDimensions {
        // Get the orthonormal basis vectors of the viewport.
        let u = basis.u;
        let v = basis.v;
        let w = basis.w;

        // The vectors that define the three dimensions of the
        // viewport in world coordinates.
        //
        // The `depth` is the vector that points from the camera to
        // the viewport plane.
        let horizontal = u * viewport.width * focus_dist;
        let vertical = v * viewport.height * focus_dist;
        let depth = w * focus_dist;

        ViewportDimensions {
            horizontal,
            vertical,
            depth,
        }
    }

    fn lower_left_corner(&self) -> Vec3 {
        (self.horizontal / 2.) + (self.vertical / 2.) + self.depth
    }
}

pub struct Camera {
    origin: Vec3,
    lower_left_corner: Vec3,
    viewport_dimensions: ViewportDimensions,
    basis: Basis,
    lens: Lens,
}

impl Camera {
    // vertical field of view in degrees
    pub fn new(
        orientation: Orientation,
        viewport: Viewport,
        lens: Lens,
        focus_dist: f64,
    ) -> Camera {
        // The location of the camera in world coordinates.
        let origin = orientation.lookfrom;

        // The orthonormal basis for the viewport in world coordinates.
        //
        // - `u` is directed left-to-right,
        // - `v` is directed bottom-to-top,
        // - `w` is directed close-to-far.
        //
        let w = (orientation.lookfrom - orientation.lookat).unit();
        let u = orientation.vup.cross(w).unit();
        let v = w.cross(u);
        let basis = Basis::new(u, v, w);

        // Compute the coordinates of the lower left corner of the viewport.
        //
        // Storing this vector is an optimization to make it easier to
        // compute the direction of the rays generated in the
        // `ray_through` method.
        let viewport_dimensions = ViewportDimensions::new(&basis, viewport, focus_dist);
        let lower_left_corner: Vec3 = origin - viewport_dimensions.lower_left_corner();

        Camera {
            origin,
            lower_left_corner,
            viewport_dimensions,
            basis,
            lens,
        }
    }

    // `s` and `t` are numbers between 0 and 1, representing how far
    // along the viewport axes to generate the ray.
    pub fn ray_through(&self, s: f64, t: f64) -> Ray {
        // Generate random noise in the location in the camera at
        // which the ray is generated. This random noise is scaled by
        // the lens radius. If the radius is zero, there is no noise.
        let mut rng = rand::thread_rng();
        let distribution = RandomInUnitDiskVec3::new();
        let rd = distribution.sample(&mut rng) * self.lens.radius;
        let noise = self.basis.u * rd.x + self.basis.v * rd.y;

        // The ray originates at the location of the camera with some
        // noise. The noise gives the camera a "depth of field"
        // effect.
        let origin = self.origin + noise;

        // Compute the direction in which the ray is travelling.
        //
        // The ray is directed from the origin point through the point
        // in the viewport indicated by the `s` and `t` percentages.
        let horizontal = self.viewport_dimensions.horizontal;
        let vertical = self.viewport_dimensions.vertical;
        let direction =
            self.lower_left_corner + horizontal * s + vertical * t - self.origin - noise;

        Ray::new(origin, direction)
    }
}
