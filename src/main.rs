use std::{f64::consts::PI, time, thread};

const SX: usize = 32; // resolution
const SY: usize = 32;
const STEP_SIZE: f64 = 0.007;

fn main() {
    let e = [0.0, 0.0, -1.3, 1.0];
    let cam = Camera::new(
        e, // camera center
        [0.0 - e[0], 0.0 - e[1], 0.0 - e[2], 1.0], // direction camera is looking at
        [0.0, 1.0, 0.0, 1.0], // up direction
        1.0, // f_close, distance to image plane
        15.0, // f_far, distance to far plane
        0.7, // wsize, camera size
        -0.35, // wl, left edge of camera in camera coordinates
        0.35 // wt, top edge of camera in camera coordinates
    );

    render(&cam);
}

fn render(cam: &Camera) {
    let mut t = 0.0;
    let surface = Torus{};

    loop {
        print!("{esc}c", esc = 27 as char);
        render_frame(&cam, &surface, t);
        thread::sleep(time::Duration::from_millis(2));
        t += 0.01;
    }
}

fn render_frame(cam: &Camera, surface: &impl Surface, t: f64) {
    let domain = surface.get_domain();
    let mut u = domain.0.0;
    let mut v = domain.1.0;
    let u_step = (domain.0.1 - domain.0.0) * STEP_SIZE;
    let v_step = (domain.1.1 - domain.1.0) * STEP_SIZE;

    let mut min_max_i = (f64::INFINITY, f64::NEG_INFINITY);
    let mut i_buff = [[0f64; SX]; SY]; // intensity buffer
    let mut d_buff = [[f64::INFINITY; SX]; SY]; // depth buffer

    while u <= domain.0.1 {
        while v <= domain.1.1 {
            let mut p = surface.phi(u, v); // get point from surface
            
            // normal of surface at u,v
            let mut n = surface.normal_field(u, v);

            rotate_x(&mut p, t/2.0);
            rotate_y(&mut p, t);
            rotate_z(&mut p, t.sin().abs());
            rotate_x(&mut n, t/2.0);
            rotate_y(&mut n, t);
            rotate_z(&mut n, t.sin().abs());

            // get point on image plane and pseudodepth
            let x = cam.perspective_transform(p);

            // continue if point not in vision
            let in_horiz_vision = cam.wl < x[0] && x[0] < cam.wl + cam.wsize;
            let in_vert_vision = cam.wt - cam.wsize < x[1] && x[1] < cam.wt;
            if !in_horiz_vision || !in_vert_vision { v += v_step; continue; }

            // convert to pixel
            let i = (SY as f64 * (x[0] - cam.wl) / (cam.wsize)) as usize;
            let j = (SX as f64 * (x[1] - cam.wt + cam.wsize) / (cam.wsize)) as usize;
            if d_buff[i][j] < x[2] { v += v_step; continue; }

            // compute brightness for pixel
            let mut d = x.clone();
            normalize(&mut d);

            let intensity = dot(&n, &d);
            if intensity <= 0.0 { v += v_step; continue; } // continue if point invisible
            if intensity < min_max_i.0 { min_max_i.0 = intensity }
            if intensity > min_max_i.1 { min_max_i.1 = intensity }
            
            i_buff[i][j] = intensity;
            d_buff[i][j] = x[2];

            v += v_step;
        }
        u += u_step;
        v = 0.0;
    }

    // now draw the i_buf
    for i in 0..SY {
        for j in 0..SX {
            // normalized intensity
            print!("{}", match (i_buff[i][j] - min_max_i.0) / min_max_i.1 {
                x if x <= 0.0 => ' ',
                x if x <= 0.1 => '.',
                x if x <= 0.2 => ',',
                x if x <= 0.3 => ':',
                x if x <= 0.4 => 'i',
                x if x <= 0.5 => 'l',
                x if x <= 0.6 => 'w',
                x if x <= 0.7 => 'W',
                x if x <= 0.8 => 'B',
                x if x <= 0.9 => '@',
                x if x <= 1.0 => '$',
                _ => '\0'
            });
        }

        println!();
    }
}

fn rotate_x(v: &mut [f64; 4], theta: f64) {
    let ct = theta.cos();
    let st = theta.sin();

    *v = [
        v[0],
        v[1] * ct - v[2] * st,
        v[1] * st + v[2] * ct,
        1.0
    ];
}

fn rotate_y(v: &mut [f64; 4], theta: f64) {
    let ct = theta.cos();
    let st = theta.sin();

    *v = [
        v[0] * ct + v[2] * st,
        v[1],
        -v[0] * st + v[2] * ct,
        1.0
    ];
}

fn rotate_z(v: &mut [f64; 4], theta: f64) {
    let ct = theta.cos();
    let st = theta.sin();

    *v = [
        v[0] * ct - v[1] * st,
        v[0] * st + v[1] * ct,
        v[2],
        1.0
    ];
}

fn dot(p1: &[f64; 4], p2: &[f64; 4]) -> f64 {
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]
}

// cross product ignoring the homogenous component
fn cross(u: &[f64; 4], v: &[f64; 4]) -> [f64; 4] {
    return [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
        1.0
    ];
}

// normalizes the vector ignoring the 4th homogeneous component
fn normalize(v: &mut[f64; 4]) {
    let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    v[0] = v[0] / norm;
    v[1] = v[1] / norm;
    v[2] = v[2] / norm;
}

fn mat_vec_mult(m: [[f64; 4]; 4], v: [f64; 4]) -> [f64; 4] {
    return [
        m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2] + m[0][3]*v[3],
        m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2] + m[1][3]*v[3],
        m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2] + m[2][3]*v[3],
        m[3][0]*v[0] + m[3][1]*v[1] + m[3][2]*v[2] + m[3][3]*v[3]
    ];
}

#[allow(dead_code)]
struct Camera {
    e: [f64; 4], // camera origin
    u: [f64; 4], // camera orthonormal basis vectors
    v: [f64; 4],
    w: [f64; 4],
    f_close: f64, // distance to image plane
    f_far: f64, // distance to far plane
    a: f64, // precomputed pseudodepth values
    b: f64,
    w2c: [[f64; 4]; 4], // matrix to convert to camera coordinates
    wsize: f64, // window size in distance units
    wl: f64, // left edge of camera in camera coordinates
    wt: f64 // top edge of camera in camera coordinates
}

impl Camera {
    fn new(e: [f64; 4], g: [f64; 4], up: [f64; 4], f_close: f64, f_far: f64, wsize: f64, wl: f64, wt: f64) -> Camera {
        let mut g = g.to_owned();
        normalize(&mut g);
        let w = [-g[0], -g[1], -g[2], 1.0];

        let mut u = cross(&up, &w);
        normalize(&mut u);

        let v = cross(&w, &u);

        let a = (-1.0 / f_close) * ((f_close + f_far)/(f_close - f_far));
        let b = (2.0 * f_far) / (f_close - f_far);

        return Camera {
            e, u, v, w, f_close, f_far, a, b,
            wsize, wl, wt,
            w2c: [
                [u[0], u[1], u[2], -dot(&u, &e)],
                [v[0], v[1], v[2], -dot(&v, &e)],
                [w[0], w[1], w[2], -dot(&w, &e)],
                [0.0, 0.0, 0.0, 1.0]
            ]
        }
    }

    // transforms a point from world coordinates to the image plane with pseudodepth
    fn perspective_transform(&self, p: [f64; 4]) -> [f64; 4] {
        let pc = mat_vec_mult(self.w2c, p); // camera coordinates
        return [
            self.f_close * pc[0]/pc[2],
            self.f_close * pc[1]/pc[2],
            self.f_close * (self.a*pc[2] + self.b) / pc[2],
            1.0
        ];
    }
}

trait Surface {
    // get domain of parametrization: ((u_min, u_max), (v_min, v_max))
    fn get_domain(&self) -> ((f64, f64), (f64, f64));

    // parametrization of surface
    fn phi(&self, u: f64, v: f64) -> [f64; 4];

    // normal field for this surface
    fn normal_field(&self, u: f64, v: f64) -> [f64; 4];
}

struct Sphere {}
impl Surface for Sphere {
    fn get_domain(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, PI), (0.0, 2.0*PI))
    }

    fn phi(&self, u: f64, v: f64) -> [f64; 4] {
        // generate point from sphere parametric equation
        let r = 0.3;
        return [
            r * u.sin() * v.cos(),
            r * u.sin() * v.sin(),
            r * u.cos() + 0.14,
            1.0
        ];
    }

    // normal vector field for sphere
    fn normal_field(&self, u: f64, v: f64) -> [f64; 4] {
        let mut n = self.phi(u, v);
        normalize(&mut n);
        return n;
    }
}

struct Torus {}
impl Surface for Torus {
    fn get_domain(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, 2.0*PI), (0.0, 2.0*PI))
    }

    fn normal_field(&self, u: f64, v: f64) -> [f64; 4] {
        let cu = u.cos();
        let su = u.sin();
        let cv = v.cos();
        let sv = v.sin();

        let mut n = [
            -cu*cv,
            -cu*sv,
            su * sv*sv - su*cv*cv,
            1.0
        ];

        normalize(&mut n);
        return n;
    }

    fn phi(&self, u: f64, v: f64) -> [f64; 4] {
        // generate point from sphere parametric equation
        let r1 = 0.3;
        let r2 = 0.1;
        return [
            (r1 + r2*v.cos())*u.cos(),
            (r1+r2*v.cos())*u.sin(),
            r2*v.sin(),
            1.0
        ];
    }
}

