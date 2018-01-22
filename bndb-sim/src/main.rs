extern crate cgmath;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate bincode;
extern crate rand;
extern crate zip;
extern crate bndb_core;

use cgmath::prelude::*;
use cgmath::{dot};
use serde::ser::*;
use bincode::*;
use rand::{random, Rng};
use bndb_core::*;

#[derive(Debug)]
enum Tree<'b> {
    Internal {
        center: Vec3,
        extents: Vec3,
        mass_center: Vec3,
        total_mass: Real,
        bodies: Vec<&'b Body>,
        children: Vec<Option<Box<Tree<'b>>>>
    },
    External {
        body: &'b Body
    }
}

fn new_octree_center(i: usize, center: Vec3, extents: Vec3) -> Vec3 {
    let d = match i {
        0 => Vec3::new( 1.0, 1.0, 1.0),
        1 => Vec3::new(-1.0, 1.0, 1.0),
        2 => Vec3::new(-1.0,-1.0, 1.0),
        3 => Vec3::new( 1.0,-1.0, 1.0),

        4 => Vec3::new( 1.0, 1.0,-1.0),
        5 => Vec3::new(-1.0, 1.0,-1.0),
        6 => Vec3::new(-1.0,-1.0,-1.0),
        7 => Vec3::new( 1.0,-1.0,-1.0),
        _ => panic!("invalid octree index")
    };
    center + d.mul_element_wise(extents/2.0)
}

fn select_octant_for_point(wp: Vec3, center: Vec3, extents: Vec3) -> usize {
    let p = wp-center;
    match (p.x > 0.0, p.y > 0.0, p.z > 0.0) {
        (true,  true,  true) => 0,
        (false, true,  true) => 1,
        (false, false, true) => 2,
        (true,  false, true) => 3,
        (true,  true,  false) => 4,
        (false, true,  false) => 5,
        (false, false, false) => 6,
        (true,  false, false) => 7,
    }
}

fn gravity_between(body: &Body, pos: Vec3, mass: Real) -> Vec3 {
    let rv = pos-body.position;
    let ir2 = 1.0 / dot(rv,rv);
    let rv = rv * ir2.sqrt();
    rv * (G * body.mass * mass * ir2) 
}

impl<'b> Tree<'b> {
    fn construct(bodies: &Vec<&'b Body>, center: Vec3, extents: Vec3) -> Tree<'b> {
        if bodies.len() == 1 {
            Tree::External { body: &bodies[0] }
        } else {
            let mut mass_center = Vec3::zero();
            let mut total_mass = 0.0;
            let mut nodeb : Vec<&'b Body> = Vec::new();
            let mut chlb : [Option<Vec<&'b Body>>; 8] = [ None, None, None, None, None, None, None, None ];
            for b in bodies.iter() {
                mass_center += b.position * b.mass;
                total_mass += b.mass;
                nodeb.push(b);
                chlb[select_octant_for_point(b.position, center, extents)].get_or_insert_with(|| Vec::new()).push(b);
            }
            mass_center /= total_mass;
            Tree::Internal {
                mass_center, total_mass, center, extents, bodies: nodeb,
                children: chlb.iter().enumerate().map(|(i, x)| x.as_ref().map(|bs| Box::new(Tree::construct(bs, new_octree_center(i, center, extents), extents/2.0)))).collect()
            }
        }
    }

    fn calculate_force(&self, body: &Body) -> Vec3 {
        match self {
            &Tree::External { body: obody } => {
                if body.index == obody.index { Vec3::zero() } else {
                    gravity_between(body, obody.position, obody.mass)
                }
            },
            &Tree::Internal { mass_center, total_mass, ref children, extents, .. } => {
                if extents.x.max(extents.y.max(extents.z)) / body.position.distance(mass_center) > 0.5 {
                    children.iter().map(|oc| oc.as_ref().map_or_else(|| Vec3::zero(), |c| c.calculate_force(body))).fold(Vec3::zero(), |a,b| a+b)
                } else {
                    gravity_between(body, mass_center, total_mass)
                }
            }
        }
    }
}

// tree N-bodies simulation
//      for each 𝚫t
//          build tree
//          calculate forces
//          move bodies
//          write simulation data into a file for playback



fn main() {
    
    /*let bodies = vec![
        Body { position: Vec3::new(1.0, 0.0, 1.0), mass: 1.0, index: 0 },
        Body { position: Vec3::new(-1.0, 0.0, 1.0), mass: 1.0, index: 1 },
    ];*/
    let mut bodies = Vec::new();
    let mut rnd = rand::thread_rng();
    for i in 0..100 {
        let theta: Real = rnd.gen::<Real>()*2.0*pi();
        let r: Real = rnd.gen::<Real>()*6.0 + 2.0;
        let y: Real = rnd.gen::<Real>()*2.0 - 1.0;
        let p = Vec3::new(r*theta.cos(), r*theta.sin(), y);
        bodies.push(Body {
            index: i,
            position: p,
            velocity: Vec3::zero(),
            mass: 1.0
        });
    }
    let mut next_bodies = Vec::new();

    let mut out = zip::write::ZipWriter::new(std::fs::OpenOptions::new()
                                             .write(true).create(true).truncate(true).open("out.zip").expect("open output file"));

    // simulate time
    let dt: Real = 0.00001;
    for step in 0..50000 {
        if step % 5000 == 0 { println!("step {}", step); }
        let (cb, nb) = if step % 2 == 0 { (&bodies, &mut next_bodies) } else { (&next_bodies, &mut bodies) };
        out.start_file(format!("{}", step), zip::write::FileOptions::default()).expect("start step file");
        bincode::serialize_into(&mut out, cb, bincode::Infinite).expect("serialize step data");
        let tree = Tree::construct(&cb.iter().map(|x| x).collect(), Vec3::zero(), Vec3::new(200.0, 200.0, 200.0));
        *nb = cb.iter().map(|body| {
            let f = tree.calculate_force(body);
            let a = f / body.mass;
            Body {
                position: body.position + body.velocity * dt,
                velocity: body.velocity + a * dt,
                .. *body
            }
        }).collect::<Vec<_>>();
        if step % 5000 == 0 {
            println!("system momentum = {:?}", nb.iter().map(|b| b.velocity*b.mass).fold(Vec3::zero(), |a,b| a+b));
        }
    }
}
