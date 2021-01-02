pub struct Duration {
    pub start: f64,
    pub end: f64,
}

impl Duration {
    pub fn new(start: f64, end: f64) -> Duration {
        Duration { start, end }
    }

    pub fn len(&self) -> f64 {
        self.end - self.start
    }

    pub fn max() -> Duration {
        // FIXME This is not good.
        Duration {
            start: -10000.,
            end: 10000.,
        }
    }
}
