#![feature(vec_resize_default)]
extern crate cgmath;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate bincode;
extern crate rand;
extern crate tar;
extern crate bndb_core;
extern crate rayon;
extern crate image;

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
                children: chlb.par_iter().enumerate().map(|(i, x)| x.as_ref().map(|bs| Box::new(Tree::construct(bs, new_octree_center(i, center, extents), extents/2.0)))).collect()
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
                    /*let mut s = Vec3::zero();
                    for c in children {
                        if let &Some(ref ch) = c {
                            s += ch.calculate_force(body);
                        }
                    }
                    s*/
                } else {
                    gravity_between(body, mass_center, total_mass)
                }
            }
        }
    }
}

fn hsv_to_rgb((h,s,v) : (Real,Real,Real)) -> (Real,Real,Real) {
    let r = (((h * 6.0 - 3.0).abs() - 1.0).min(1.0).max(0.0) - 1.0) * s + 1.0;
    let g = ((2.0 - (h * 6.0 - 2.0).abs()).min(1.0).max(0.0) - 1.0) * s + 1.0;
    let b = ((2.0 - (h * 6.0 - 4.0).abs()).min(1.0).max(0.0) - 1.0) * s + 1.0;
    (r*v,g*v,b*v)
}

fn render_image(bodies: &Vec<Body>) -> image::RgbImage {
    use cgmath::*;
    let width:i32 = 640;
    let height:i32 = 480;
    let size = 8;
    let T =
        Matrix4::from_translation(Vector3::new(320.0, 240.0, 0.0))*
        Matrix4::from_nonuniform_scale(320.0, 240.0, 1.0)*
        perspective(cgmath::Deg(60.0), (width as Real) / (height as Real), 0.1, 100.0)*
        Matrix4::look_at(cgmath::Point3::new(15.0, 10.0, -30.0), cgmath::Point3::new(0.0, 0.0, 0.0), cgmath::Vector3::unit_y()); 

    let mut img = image::ImageBuffer::new(width as u32, height as u32);
    let mut pxs: Vec<(Real,Real,Real)> = Vec::new();
    pxs.resize_default((width * height) as usize);

    for b in bodies {
        let mut p = T * cgmath::Vector4::new(b.position.x, b.position.y, b.position.z, 1.0);
        let x = (p.x / p.w) as i32;
        let y = (p.y / p.w) as i32;
        if x < 0 || x >= width || y < 0 || y >= height { continue; }
        for iy in 0..(size+1) {
            let dy = (iy as i32 * 2 - size as i32)/2;
            let v = ((dy as Real) / (size as Real));
            if (y+dy) < 0 || (y+dy) >= height { continue; }
            for ix in 0..(size+1) {
                let dx = (ix as i32 * 2 - size as i32)/2;
                let u = ((dx as Real) / (size as Real));
                let idx = ((x+dx) + (y+dy)*(width as i32)) as usize;
                if idx >= pxs.len() { continue; }
                let sat = (1.0-(u*u + v*v)) * 0.0015;
                let (r,g,b) = hsv_to_rgb(( (dot(b.velocity,b.velocity)*10.0).sin()*0.5+0.5 ,0.9,sat));
                pxs[idx].0 += r;
                pxs[idx].1 += g;
                pxs[idx].2 += b;
            }
        }
    }

    for (x, y, px) in img.enumerate_pixels_mut() {
        let r: u8 = (pxs[(x + y*(width as u32)) as usize].0 * 255.0).min(255.0) as u8;
        let g: u8 = (pxs[(x + y*(width as u32)) as usize].1 * 255.0).min(255.0) as u8;
        let b: u8 = (pxs[(x + y*(width as u32)) as usize].2 * 255.0).min(255.0) as u8;
        *px = image::Rgb([r, g, b]);
    }
    img
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
    /*for i in 0..32_000 {
        let theta: Real = rnd.gen::<Real>()*2.0*pi();
        let r: Real = rnd.gen::<Real>()*8.0 + 8.0;
        let y: Real = (rnd.gen::<Real>()*2.0 - 1.0)*0.3;
        let p = Vec3::new(r*theta.cos(), y, r*theta.sin());
        bodies.push(Body {
            index: i,
            position: p,
            velocity: p.normalize().cross(Vec3::unit_y()) * 0.1,
            mass: 25.0
        });
    }
    let bc = bodies.len();
    for i in 0..32_000 {
        let theta: Real = rnd.gen::<Real>()*2.0*pi();
        let r: Real = rnd.gen::<Real>()*8.0 + 8.0;
        let y: Real = (rnd.gen::<Real>()*2.0 - 1.0)*0.3;
        let p = Vec3::new(r*theta.cos(), r*theta.sin(), y);
        bodies.push(Body {
            index: bc+i,
            position: p,
            velocity: p.normalize().cross(Vec3::unit_z()) * 0.1,
            mass: 25.0
        });
    }
    let bc = bodies.len();
    for i in 0..32_000 {
        let theta: Real = rnd.gen::<Real>()*2.0*pi();
        let r: Real = rnd.gen::<Real>()*8.0 + 8.0;
        let y: Real = (rnd.gen::<Real>()*2.0 - 1.0)*0.3;
        let p = Vec3::new(y, r*theta.cos(), r*theta.sin());
        bodies.push(Body {
            index: bc+i,
            position: p,
            velocity: p.normalize().cross(Vec3::unit_x()) * 0.1,
            mass: 25.0
        });
    }*/
    for i in 0..200_000 {
        let theta: Real = rnd.gen::<Real>()*2.0*pi();
        let phi: Real = rnd.gen::<Real>()*2.0*pi();
        let r: Real = rnd.gen::<Real>()*2.0 + 15.0;
        let p = Vec3::new(r*theta.sin()*phi.cos(), r*theta.sin()*phi.sin(), r*theta.cos());
        bodies.push(Body {
            index: i,
            position: p,
            velocity: p.normalize().cross(Vec3::unit_y()) * 0.1,
            mass: 5.0
        });
    }
    let bc = bodies.len();
    for i in 0..32_000 {
        let theta: Real = rnd.gen::<Real>()*2.0*pi();
        let r: Real = rnd.gen::<Real>()*8.0 + 13.0;
        let y: Real = (rnd.gen::<Real>()*2.0 - 1.0)*0.4;
        let p = Vec3::new(r*theta.cos(), y, r*theta.sin());
        bodies.push(Body {
            index: bc+i,
            position: p,
            velocity: p.normalize().cross(Vec3::unit_y()) * -0.2,
            mass: 8.0
        });
    }
    let bc = bodies.len();
    bodies.push(Body {
        index: bc+1,
        position: Vec3::zero(),
        velocity: Vec3::zero(),
        mass: 6e9
    });
    let end_init = Instant::now();
    println!("init took: {:?}", (end_init - start_init));
    let mut next_bodies = Vec::new();

    let mut dest_path = std::path::PathBuf::from(std::env::args().skip(1).next().expect("destination directory as first argument"));
    //let mut out = tar::Builder::new(std::fs::OpenOptions::new()
    //                                        .write(true).create(true).truncate(true).open(filename).expect("open output file"));
    //let mut out = zip::write::ZipWriter::new(std::fs::OpenOptions::new()
    //                                         .write(true).create(true).truncate(true).open(filename).expect("open output file"));

    //let zipopts = zip::write::FileOptions::default().compression_method(zip::CompressionMethod::Bzip2);
    //out.start_file("metadata", zipopts).expect("start metadata file");

    let simm = SimulationMetadata {
        delta_time: 0.5,
        steps: 100_000
    };

    /*println!("serializing metadata");
    let simm_size = bincode::serialized_size(&simm);
    let mut header = tar::Header::new_gnu();
    header.set_size(simm_size);
    header.set_cksum();
    let mut staging_data = Vec::new();
    staging_data.resize_default(simm_size as usize);
    bincode::serialize_into(&mut staging_data, &simm, bincode::Infinite).expect("serialize metadata");
    println!("data = {:?}", staging_data);
    out.append_data(&mut header, "metadata", &staging_data[..]);*/

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
            //out.start_file(format!("{}", step / substeps), zipopts).expect("start step file");
            /*let positions = cb.iter().map(|x| x.position).collect::<Vec<_>>();
            let size = bincode::serialized_size(&positions);
            header.set_size(size);
            header.set_cksum();
            staging_data.resize_default(size as usize);
            bincode::serialize_into(&mut staging_data, &positions, bincode::Infinite).expect("serialize step data");
            out.append_data(&mut header, format!("{}", step/substeps), &staging_data[..]);*/
            let img = image::ImageRgb8(render_image(&cb));
            dest_path.push(format!("step{}.png", step/substeps));
            img.save(&mut std::fs::File::create(&dest_path).unwrap(), image::PNG).unwrap();
            dest_path.pop();
        }
//        let start_tree_build = Instant::now();
        let tree = Tree::construct(&cb.iter().map(|x| x).collect(), Vec3::zero(), Vec3::new(500.0, 500.0, 500.0));
 //       let end_tree_build = Instant::now();
  //      println!("\ttree build took: {:?}", end_tree_build-start_tree_build);
        cb.par_iter().map(|body| {
            let f = tree.calculate_force(body);
            let a = f / body.mass;
            Body {
                position: body.position + body.velocity * dt,
                velocity: body.velocity + a * dt,
                .. *body
            }
        }).collect_into(nb);
        //nb.clear();
        /*cb.par_iter().map(|body| {
            let f = tree.calculate_force(body);
            let a = f / body.mass;
            Body {
                position: body.position + body.velocity * dt,
                velocity: body.velocity + a * dt,
                .. *body
            }
        }).for_each(|b| nb.push(b));*/
        /*for body in cb {
            let f = tree.calculate_force(body);
            let a = f / body.mass;
            nb.push(Body {
                position: body.position + body.velocity * dt,
                velocity: body.velocity + a * dt,
                .. *body
            });
        }*/
   //     let end_update = Instant::now();
    //    println!("\tupdate took: {:?}", end_update-end_tree_build);
        if step % 500 == 0 {
            let p = nb.par_iter().map(|b| b.velocity*b.mass).sum();
            println!("system momentum = {:?},\n\tΔp = {:?}", p, last_momentum-p);
            last_momentum = p;
        }
    }

    //out.finish().expect("finish archive");
}
