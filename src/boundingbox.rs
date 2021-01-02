use crate::ray::Ray;
use crate::vector::Vec3;

#[derive(Copy, Clone)]
pub enum BoundingBox {
    Empty,
    Nonempty { minimum: Vec3, maximum: Vec3 },
}

impl BoundingBox {
    // TODO Can this be refactored with `Hittable::hits()`?
    pub fn intersects(&self, ray: &Ray, time_bounds: (f64, f64)) -> bool {
        match self {
            BoundingBox::Empty => false,
            BoundingBox::Nonempty { minimum, maximum } => {
                let (mut earliest, mut latest) = time_bounds;

                let origin = ray.origin;
                let direction = ray.direction;

                let x_t0 = ((minimum.x - origin.x) / direction.x)
                    .min((maximum.x - origin.x) / direction.x);
                let x_t1 = ((minimum.x - origin.x) / direction.x)
                    .max((maximum.x - origin.x) / direction.x);
                earliest = earliest.min(x_t0);
                latest = latest.max(x_t1);
                if latest <= earliest {
                    return false;
                }

                let y_t0 = ((minimum.y - origin.y) / direction.y)
                    .min((maximum.y - origin.y) / direction.y);
                let y_t1 = ((minimum.y - origin.y) / direction.y)
                    .max((maximum.y - origin.y) / direction.y);
                earliest = earliest.min(y_t0);
                latest = latest.max(y_t1);
                if latest <= earliest {
                    return false;
                }

                let z_t0 = ((minimum.z - origin.z) / direction.z)
                    .min((maximum.z - origin.z) / direction.z);
                let z_t1 = ((minimum.z - origin.z) / direction.z)
                    .max((maximum.z - origin.z) / direction.z);
                earliest = earliest.min(z_t0);
                latest = latest.max(z_t1);
                if latest <= earliest {
                    return false;
                }

                true
            }
        }
    }

    pub fn surrounding(box0: BoundingBox, box1: BoundingBox) -> BoundingBox {
        match (box0, box1) {
            (BoundingBox::Empty, BoundingBox::Empty) => BoundingBox::Empty,
            (BoundingBox::Empty, BoundingBox::Nonempty { minimum, maximum }) => {
                BoundingBox::Nonempty { minimum, maximum }
            }
            (BoundingBox::Nonempty { minimum, maximum }, BoundingBox::Empty) => {
                BoundingBox::Nonempty { minimum, maximum }
            }
            (
                BoundingBox::Nonempty {
                    minimum: min0,
                    maximum: max0,
                },
                BoundingBox::Nonempty {
                    minimum: min1,
                    maximum: max1,
                },
            ) => {
                let min_x = min0.x.min(min1.x);
                let min_y = min0.y.min(min1.y);
                let min_z = min0.z.min(min1.z);
                let minimum = Vec3::new(min_x, min_y, min_z);

                let max_x = max0.x.max(max1.x);
                let max_y = max0.y.max(max1.y);
                let max_z = max0.z.max(max1.z);
                let maximum = Vec3::new(max_x, max_y, max_z);

                BoundingBox::Nonempty { minimum, maximum }
            }
        }
    }
}
