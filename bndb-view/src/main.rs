extern crate gl;
extern crate glutin;
extern crate cgmath;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate bincode;
extern crate zip;
extern crate bndb_core;

use cgmath::prelude::*;
use cgmath::{dot, Point3, Vector3, Matrix4, Rad};
use serde::de::*;
use bincode::*;
use bndb_core::*;
use glutin::GlContext;

use gl::types::*;

use std::{mem, ptr, str};
use std::ffi::CString;

// Vertex data
static VERTEX_DATA: [GLfloat; 6] = [0.0, 0.5, 0.5, -0.5, -0.5, -0.5];

// Shader sources
static VS_SRC: &'static str = "
#version 150
in vec3 position;
uniform mat4 vp;
void main() {
    gl_Position = vp*vec4(position, 1.0);
}";

static FS_SRC: &'static str = "
#version 150
out vec4 out_color;

void main() {
    out_color = vec4(0.8, 0.8, 1.0, 0.5);
}";

fn compile_shader(src: &str, ty: GLenum) -> GLuint {
    let shader;
    unsafe {
        shader = gl::CreateShader(ty);
        // Attempt to compile the shader
        let c_str = CString::new(src.as_bytes()).unwrap();
        gl::ShaderSource(shader, 1, &c_str.as_ptr(), ptr::null());
        gl::CompileShader(shader);

        // Get the compile status
        let mut status = gl::FALSE as GLint;
        gl::GetShaderiv(shader, gl::COMPILE_STATUS, &mut status);

        // Fail on error
        if status != (gl::TRUE as GLint) {
            let mut len = 0;
            gl::GetShaderiv(shader, gl::INFO_LOG_LENGTH, &mut len);
            let mut buf = Vec::with_capacity(len as usize);
            buf.set_len((len as usize) - 1); // subtract 1 to skip the trailing null character
            gl::GetShaderInfoLog(
                shader,
                len,
                ptr::null_mut(),
                buf.as_mut_ptr() as *mut GLchar,
            );
            panic!(
                "{}",
                str::from_utf8(&buf)
                    .ok()
                    .expect("ShaderInfoLog not valid utf8")
            );
        }
    }
    shader
}

fn link_program(vs: GLuint, fs: GLuint) -> GLuint {
    unsafe {
        let program = gl::CreateProgram();
        gl::AttachShader(program, vs);
        gl::AttachShader(program, fs);
        gl::LinkProgram(program);
        // Get the link status
        let mut status = gl::FALSE as GLint;
        gl::GetProgramiv(program, gl::LINK_STATUS, &mut status);

        // Fail on error
        if status != (gl::TRUE as GLint) {
            let mut len: GLint = 0;
            gl::GetProgramiv(program, gl::INFO_LOG_LENGTH, &mut len);
            let mut buf = Vec::with_capacity(len as usize);
            buf.set_len((len as usize) - 1); // subtract 1 to skip the trailing null character
            gl::GetProgramInfoLog(
                program,
                len,
                ptr::null_mut(),
                buf.as_mut_ptr() as *mut GLchar,
            );
            panic!(
                "{}",
                str::from_utf8(&buf)
                    .ok()
                    .expect("ProgramInfoLog not valid utf8")
            );
        }
        program
    }
}

struct App {
    program: GLuint,
    vao: GLuint, vbo: GLuint,
    vp_unif: GLint,
    vertex_data: Vec<Point3<GLfloat>>,
    inf: zip::read::ZipArchive<std::fs::File>,
    current_step: usize, t: GLfloat,
    simm: SimulationMetadata
}

impl App {
    fn init() -> App {
        let vs = compile_shader(VS_SRC, gl::VERTEX_SHADER);
        let fs = compile_shader(FS_SRC, gl::FRAGMENT_SHADER);
        let program = link_program(vs, fs);

        let mut vao = 0;
        let mut vbo = 0;
        let mut vp_unif = 0;

        unsafe {
            // Create Vertex Array Object
            gl::GenVertexArrays(1, &mut vao);
            gl::BindVertexArray(vao);

            // Create a Vertex Buffer Object and copy the vertex data to it
            gl::GenBuffers(1, &mut vbo);
            gl::BindBuffer(gl::ARRAY_BUFFER, vbo);
            /*gl::BufferData(
                gl::ARRAY_BUFFER,
                (VERTEX_DATA.len() * mem::size_of::<GLfloat>()) as GLsizeiptr,
                mem::transmute(&VERTEX_DATA[0]),
                gl::STATIC_DRAW,
                );*/

            // Use shader program
            gl::UseProgram(program);
//            gl::BindFragDataLocation(program, 0, CString::new("out_color").unwrap().as_ptr());

            // Specify the layout of the vertex data
            let pos_attr = gl::GetAttribLocation(program, CString::new("position").unwrap().as_ptr());
            gl::EnableVertexAttribArray(pos_attr as GLuint);
            gl::VertexAttribPointer(
                pos_attr as GLuint,
                3,
                gl::FLOAT,
                gl::FALSE as GLboolean,
                0,
                ptr::null(),
                );

            vp_unif = gl::GetUniformLocation(program, CString::new("vp").unwrap().as_ptr());
            println!("{}", vp_unif);
        }

        let mut inf = zip::read::ZipArchive::new(std::fs::File::open(std::env::args().skip(1).next().expect("file path as first argument")).expect("open file")).expect("open archive");
        let simm = bincode::deserialize_from(&mut inf.by_name("metadata").expect("open metadata file"), bincode::Infinite).expect("load sim metadata");
        println!("simm = {:?}", simm);

        let positions: Vec<Vec3> = bincode::deserialize_from(&mut inf.by_name(&format!("{}", 0)).expect("open step file"), bincode::Infinite).expect("load data");
        unsafe { gl::BufferData(
            gl::ARRAY_BUFFER,
            (positions.len() * mem::size_of::<Vec3>()) as GLsizeiptr,
            mem::transmute(positions.as_ptr()),
            gl::DYNAMIC_DRAW,
            );
        }



        App {
            program, vao, vbo, vp_unif, inf, simm, current_step: 0, t: 0.0,
            vertex_data: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(-1.0, 0.0, 0.0),
                Point3::new(0.0, -1.0, 0.0),
                Point3::new(0.0, 0.0, 1.0),
                Point3::new(0.0, 0.0, -1.0),
            ]
        }
    }

    unsafe fn draw(&mut self) {
        gl::BindVertexArray(self.vao);
        gl::BindBuffer(gl::ARRAY_BUFFER, self.vbo);


        let positions: Vec<Vec3> = bincode::deserialize_from(&mut self.inf.by_name(&format!("{}", self.current_step)).expect("open step file"), bincode::Infinite).expect("load data");
        gl::BufferData(
            gl::ARRAY_BUFFER,
            (positions.len() * mem::size_of::<Vec3>()) as GLsizeiptr,
            mem::transmute(positions.as_ptr()),
            gl::DYNAMIC_DRAW,
        );

        self.current_step+=2;
        if self.current_step > self.simm.steps-1 {
            println!("loop");
            self.current_step = 0;
        }

        gl::Enablei(gl::BLEND, 0);
        gl::BlendEquation(gl::FUNC_ADD);
        gl::BlendFunc(gl::SRC_ALPHA, gl::DST_ALPHA);
        gl::UseProgram(self.program);
        let v: Matrix4<GLfloat> = Matrix4::look_at(Point3::new((0.05*self.t).cos()*30.0, 14.0, (0.05*self.t).sin()*30.0), Point3::new(0.0, 0.0, 0.0), Vector3::unit_y());
        let p: Matrix4<GLfloat> = cgmath::perspective(Rad(std::f32::consts::FRAC_PI_4), 1.333333, 0.1, 1000.0);
        gl::UniformMatrix4fv(self.vp_unif, 1, 0, (p*v).as_ptr());
        gl::PointSize(3.0);
        gl::DrawArrays(gl::POINTS, 0, positions.len() as i32);
        self.t += 0.01;
    }
}


fn main() {
    let mut events_loop = glutin::EventsLoop::new();
    let window = glutin::WindowBuilder::new()
        .with_title("bndb view")
        .with_dimensions(640*4, 480*4);
    let context = glutin::ContextBuilder::new()
        .with_vsync(true);
    let gl_window = glutin::GlWindow::new(window, context, &events_loop).unwrap();

    unsafe {
        gl_window.make_current().unwrap();
    }

    unsafe {
        gl::load_with(|symbol| gl_window.get_proc_address(symbol) as *const _);
        gl::ClearColor(0.0, 0.0, 0.0, 0.0);
    }

    let mut app = App::init();

    let mut running = true;
    while running {
        events_loop.poll_events(|event| {
            match event {
                glutin::Event::WindowEvent{ event, .. } => match event {
                    glutin::WindowEvent::Closed => running = false,
                    glutin::WindowEvent::Resized(w, h) => gl_window.resize(w, h),
                    _ => ()
                },
                _ => ()
            }
        });

        unsafe {
            gl::Clear(gl::COLOR_BUFFER_BIT);
            app.draw();
        }

        gl_window.swap_buffers().unwrap();
    }
}
