extern crate cgmath;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate bincode;
extern crate rand;
extern crate tar;
extern crate bndb_core;
extern crate rayon;

use cgmath::prelude::*;
use cgmath::{dot};
use serde::ser::*;
use bincode::*;
use rand::{random, Rng};
use bndb_core::*;
use rayon::prelude::*;

use std::time::{Instant, Duration};

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
    let rv = body.position-pos;
    let ir2 = 1.0 / dot(rv,rv);
    let rv = rv * ir2.sqrt();
    -rv * (G * body.mass * mass * ir2) 
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
                let bpa = (b.position-center);
                if bpa.x.abs() > extents.x || bpa.y.abs() > extents.y || bpa.z.abs() > extents.z { continue; }
                mass_center += b.position * b.mass;
                total_mass += b.mass;
                //nodeb.push(b);
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
                    children.iter().map(|oc| oc.as_ref().map_or_else(|| Vec3::zero(), |c| c.calculate_force(body))).sum()
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
    rayon::initialize(rayon::Configuration::new().num_threads(16)).unwrap();
    println!("?");
    let start_init = Instant::now();
    let mut bodies = Vec::new();
    let mut rnd = rand::thread_rng();
    for i in 0..1600 {
        let theta: Real = rnd.gen::<Real>()*2.0*pi();
        let r: Real = rnd.gen::<Real>()*8.0 + 8.0;
        let y: Real = (rnd.gen::<Real>()*2.0 - 1.0)*0.3;
        let p = Vec3::new(r*theta.cos(), y, r*theta.sin());
        bodies.push(Body {
            index: i,
            position: p,
            velocity: p.normalize().cross(Vec3::unit_y()) * 0.1,
            mass: 50.0
        });
    }
    let bc = bodies.len();
    for i in 0..1600 {
        let theta: Real = rnd.gen::<Real>()*2.0*pi();
        let r: Real = rnd.gen::<Real>()*8.0 + 8.0;
        let y: Real = (rnd.gen::<Real>()*2.0 - 1.0)*0.3;
        let p = Vec3::new(r*theta.cos(), r*theta.sin(), y);
        bodies.push(Body {
            index: bc+i,
            position: p,
            velocity: p.normalize().cross(Vec3::unit_z()) * 0.1,
            mass: 50.0
        });
    }
    let bc = bodies.len();
    for i in 0..1600 {
        let theta: Real = rnd.gen::<Real>()*2.0*pi();
        let r: Real = rnd.gen::<Real>()*8.0 + 8.0;
        let y: Real = (rnd.gen::<Real>()*2.0 - 1.0)*0.3;
        let p = Vec3::new(y, r*theta.cos(), r*theta.sin());
        bodies.push(Body {
            index: bc+i,
            position: p,
            velocity: p.normalize().cross(Vec3::unit_x()) * 0.1,
            mass: 50.0
        });
    }
    let bc = bodies.len();
    bodies.push(Body {
            index: bc+1,
            position: Vec3::zero(),
            velocity: Vec3::zero(),
            mass: 5e9
        });
    let end_init = Instant::now();
    println!("init took: {:?}", (end_init - start_init));
    let mut next_bodies = Vec::new();

    let mut out = zip::write::ZipWriter::new(std::fs::OpenOptions::new()
                                             .write(true).create(true).truncate(true).open("out.zip").expect("open output file"));

    let zipopts = zip::write::FileOptions::default().compression_method(zip::CompressionMethod::Bzip2);
    out.start_file("metadata", zipopts).expect("start metadata file");

    let simm = SimulationMetadata {
        delta_time: 0.05,
        steps: 10_000
    };

    println!("serializing metadata");
    bincode::serialize_into(&mut out, &simm, bincode::Infinite).expect("serialize metadata");

    let mut last_momentum = Vec3::zero();
    let mut last_time = Instant::now();

    // simulate time
    let substeps = 5;
    let dt = simm.delta_time / (substeps as Real);
    println!("Δt = {}", dt);
    for step in 0..(simm.steps * substeps) {
        if step % 500 == 0 {
            let time = Instant::now();
            println!("step {} of {} ({}), Δ(real time)={:?}", step, simm.steps*substeps, step/substeps, (time - last_time));
            last_time = time;
        }
        let (cb, nb) = if step % 2 == 0 { (&bodies, &mut next_bodies) } else { (&next_bodies, &mut bodies) };
        if step % substeps == 0 {
            out.start_file(format!("{}", step / substeps), zipopts).expect("start step file");
            bincode::serialize_into(&mut out, &cb.iter().map(|x| x.position).collect::<Vec<_>>(), bincode::Infinite).expect("serialize step data");
        }
        let tree = Tree::construct(&cb.iter().map(|x| x).collect(), Vec3::zero(), Vec3::new(500.0, 500.0, 500.0));
        *nb = cb.par_iter().map(|body| {
            let f = tree.calculate_force(body);
            let a = f / body.mass;
            Body {
                position: body.position + body.velocity * dt,
                velocity: body.velocity + a * dt,
                .. *body
            }
        }).collect::<Vec<_>>();
        if step % 500 == 0 {
            let p = nb.par_iter().map(|b| b.velocity*b.mass).sum();
            println!("system momentum = {:?},\n\tΔp = {:?}", p, last_momentum-p);
            last_momentum = p;
        }
    }
}
