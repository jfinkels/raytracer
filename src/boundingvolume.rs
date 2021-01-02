use crate::boundingbox::BoundingBox;
use crate::hittable::HitRecord;
use crate::hittable::Hittable;
use crate::ray::Ray;

struct BVHNode {
    left: Box<BVHNode>,
    right: Box<BVHNode>,
    bounding_box: BoundingBox,
}

impl Hittable for BVHNode {
    fn bounding_box(&self, _time_bounds: (f64, f64)) -> BoundingBox {
        self.bounding_box
    }
    fn hits(&self, ray: &Ray, time_bounds: (f64, f64)) -> Option<HitRecord> {
        if self.bounding_box.intersects(ray, time_bounds) {
            None
        } else {
            match self.left.hits(ray, time_bounds) {
                // In this case, see if the right child has an earlier hit.
                Some(left_hit_record) => {
                    let new_time_bounds = (time_bounds.0, left_hit_record.time);
                    match self.right.hits(ray, new_time_bounds) {
                        Some(right_hit_record) => Some(right_hit_record),
                        None => Some(left_hit_record),
                    }
                }
                None => self.right.hits(ray, time_bounds),
            }
        }
    }
}
