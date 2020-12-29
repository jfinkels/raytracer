pub struct Image {
    pub width: usize,
    pub height: usize,
}

impl Image {
    pub fn new(width: usize, height: usize) -> Image {
        Image { width, height }
    }
}
