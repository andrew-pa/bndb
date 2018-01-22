extern crate cgmath;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate bincode;
extern crate rand;
extern crate zip;

use cgmath::prelude::*;
use cgmath::{dot};

pub type Real = f64;
pub type Vec3 = cgmath::Vector3<Real>;

pub fn pi() -> Real {
    3.14159265358979323846264338327950288
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Body {
    pub index: usize,
    pub position: Vec3,
    pub velocity: Vec3,
    pub mass: Real
}

pub const G: Real = 6.67408e-11;
