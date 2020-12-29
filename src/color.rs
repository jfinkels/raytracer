use crate::vector::Vec3;

#[derive(Debug)]
pub struct Color {
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

    pub fn from_vec3(v: Vec3) -> Color {
        Color::from_ratios(v.x, v.y, v.z)
    }

    pub fn to_ppm(&self) -> String {
        format!("{} {} {}", self.red, self.green, self.blue)
    }
}
