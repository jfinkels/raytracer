#[derive(Debug)]
struct Color {
    red: usize,
    green: usize,
    blue: usize,
}

impl Color {

    fn new(red: usize, green: usize, blue: usize) -> Color {
        Color {red, green, blue}
    }

    fn to_ppm(&self) -> String {
        format!("{} {} {}", self.red, self.green, self.blue)
    }
}


fn main() {

    // Image
    const IMAGE_WIDTH: usize = 0xff;
    const IMAGE_HEIGHT: usize = 0xff;

    // Render
    println!("P3");
    println!("{} {}", IMAGE_WIDTH, IMAGE_HEIGHT);
    println!("255");

    for i in 0..IMAGE_HEIGHT {
        for j in 0..IMAGE_WIDTH {
            let color = Color::new(i, j, 0x7f);
            println!("{}", color.to_ppm());
        }
    }
}
