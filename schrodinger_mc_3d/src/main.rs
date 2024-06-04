use file::save_columns_to_file;
use rand::{self, Rng};
use rand::distributions::{Distribution, Uniform};

mod file;
mod plot;
mod draw;

const H2M0: f64 = 0.0762;
const ME: f64 = 0.067;
const A1: f64 = 20.0;
const A2: f64 = 8.0;
const B1: f64 = 20.0;
const B2: f64 = 8.0;
const QR_HEIGHT: f64 = 5.0;
const QW_HEIGHT: f64 = 3.0;
const LX: f64 = A1*5.0;
const LZ: f64 = 2.0*(QW_HEIGHT + QR_HEIGHT);
const X1: f64 = LX/2.0;
const X2: f64 = LX/2.0;
const Y1: f64 = LX/2.0;
const Y2: f64 = LX/2.0;
const U0: f64 = 0.25;
const DX: f64 = 0.2;
const DX2: f64 = DX*DX;
const NX: usize = (LX/DX) as usize;
const NY: usize = NX;
const NZ: usize = (LZ/DX) as usize;
//const NZ0: usize = (2.0*QW_HEIGHT/DX) as usize;
const D: f64 = H2M0/ME/2.0;
const DTAU: f64 = DX2/4.0/D; // units 1/eV
const FINAL_TIME: f64 = 150.0; // units 1/eV
const N_PARTICLES: usize = 50000000;
const N_PARTICLES_F: f64 = N_PARTICLES as f64;
const N_TIME_STEPS: usize = (FINAL_TIME/DTAU) as usize;
const N_VALUES_CHECK: usize = 100;
const CHECK_INTERVAL: usize = N_TIME_STEPS/N_VALUES_CHECK;
const N_ANGLES_ROTATION_PLOT: usize = 20;
const D_ANGLE: f64 = 2.0*3.14159/(N_ANGLES_ROTATION_PLOT as f64);

fn quantum_ring_shape(x: f64, y: f64) -> f64 {
    let s1 = (x-X1).powi(2)/A1.powi(2)+(y-Y1).powi(2)/B1.powi(2);
    let s2 = (x-X2).powi(2)/A2.powi(2)+(y-Y2).powi(2)/B2.powi(2);
    (1.0 + s1)*(-s1).exp() - (1.0 + s2)*(-s2).exp()
}

fn quantum_ring_z(x: f64, y: f64, a: f64) -> f64 {
    2.0*QW_HEIGHT + QR_HEIGHT*quantum_ring_shape(x, y)/a
}

fn potential(x: f64, y: f64, z: f64, a: f64) -> f64 {
    if z < quantum_ring_z(x, y, a) && z > QW_HEIGHT {
        0.0
    } else {
        U0
    }
}

fn index_xy(jx: usize, jy: usize) -> usize {
    NY*jx+jy
}

fn index_xz(jx: usize, jz: usize) -> usize {
    NZ*jx+jz
}

fn index_3d(jx: usize, jy: usize, jz: usize) -> usize {
    NZ*(NY*jx+jy)+jz
}

fn mean_last_half(v: &Vec<f64>) -> f64 {
    let l2 = v.len()/2;
    let l2f = l2 as f64;
    let mut a = 0.0;
    for e in 0..l2 {
        a += v[e+l2]/l2f
    }
    a
}

fn main() {
    let radius = 2.0*A1*A2*((A1.ln()-A2.ln())/(A1*A1-A2*A2)).sqrt();
    println!("QR radius = {} nm", radius);
    let amplitude = quantum_ring_shape(radius+X1, Y1);
    println!("QR amplitude = {} U0", amplitude);
    let mut rng = rand::thread_rng();
    let x_dist = Uniform::from(0..NX);
    let y_dist = Uniform::from(0..NY);
    let z_dist = Uniform::from(0..NZ);
    let p_dist = Uniform::from(0..N_PARTICLES);
    let d_dist = Uniform::from(0..6);
    let xn: Vec<f64> = (0..NX).map(|n| n as f64 * DX).collect();
    let yn: Vec<f64> = (0..NY).map(|n| n as f64 * DX).collect();
    let zn: Vec<f64> = (0..NZ).map(|n| n as f64 * DX).collect();
    let phin: Vec<f64> = (0..N_ANGLES_ROTATION_PLOT).map(|n| n as f64 * D_ANGLE).collect();
    let cosphin: Vec<f64> = phin.clone().into_iter().map(|x| x.cos()).collect();
    let sinphin: Vec<f64> = phin.clone().into_iter().map(|x| x.sin()).collect();
    let un1d: Vec<f64> = xn.clone().into_iter().map(|x| potential(x, Y1, QW_HEIGHT+0.5*QR_HEIGHT, amplitude)).collect();

    let mut zn2d: Vec<f64> = vec![0.0; NX*NY];
    for jx in 0..NX {
        for jy in 0..NY {
            zn2d[index_xy(jx,jy)] = quantum_ring_z(xn[jx], yn[jy], amplitude) - 2.0*QW_HEIGHT
        }
    }
    draw::plot_2d(&zn2d, NX, NY, 2000, 2000, "QR_shape_xy.png", "Turbo");

    let mut un2d: Vec<f64> = vec![0.0; NX*NZ];
    for jx in 0..NX {
        for jz in 0..NZ {
            un2d[index_xz(jx,jz)] = potential(xn[jx], Y1, LZ-zn[jz], amplitude)
        }
    }
    draw::plot_2d(&un2d, NX, NZ, 2000, (2000*NZ)/NX, "Potential_shape_xz.png", "BnW");

    let mut un3d: Vec<f64> = vec![0.0; NX*NY*NZ];
    for jx in 0..NX {
        for jy in 0..NY {
            for jz in 0..NZ {
                un3d[index_3d(jx,jy,jz)] = potential(xn[jx], yn[jy], zn[jz], amplitude)
            }
        }
    }
    let p_evap: Vec<f64> = un3d.clone().into_iter().map(|u| u*DTAU).collect();

    let mut particle_distribution_xyz: Vec<f64> = vec![0.0; NX*NY*NZ];
    let mut particle_distribution_xy: Vec<f64> = vec![0.0; NX*NY];
    let mut particle_distribution_xz: Vec<f64> = vec![0.0; NX*NZ];
    let mut particle_distribution_x: Vec<f64> = vec![0.0; NX];
    
    let mut energy: f64 = 0.0;
    let mut energy_values: Vec<f64> = Vec::new();
    let mut tau_values: Vec<f64> = Vec::new();
    let mut evaporated_particles: f64 = 0.0;
    let mut tau_current: f64 = 0.0;
    let mut particle_array_x = vec![1; N_PARTICLES];
    let mut particle_array_y = vec![1; N_PARTICLES];
    let mut particle_array_z = vec![1; N_PARTICLES];
    let mut particle_distances = vec![1; N_PARTICLES];
    let mut r1: f64;
    let mut p1: usize;

    for jp in 0..N_PARTICLES {
        particle_array_x[jp] = x_dist.sample(&mut rng);
        particle_array_y[jp] = y_dist.sample(&mut rng);
        particle_array_z[jp] = z_dist.sample(&mut rng);
    }

    for jt in 0..N_TIME_STEPS {
        for jp in 0..N_PARTICLES {
            p1 = d_dist.sample(&mut rng);
            match p1 {
                0 => {particle_array_x[jp] -= 1},
                1 => {particle_array_x[jp] += 1},
                2 => {particle_array_y[jp] -= 1},
                3 => {particle_array_y[jp] += 1},
                4 => {particle_array_z[jp] -= 1},
                5 => {particle_array_z[jp] += 1},
                _ => {},
            }
            if particle_array_x[jp] < NX && particle_array_y[jp] < NY && particle_array_z[jp] < NZ {
                r1 = rng.gen();
                if r1 < p_evap[index_3d(particle_array_x[jp],particle_array_y[jp], particle_array_z[jp])] {
                    evaporated_particles += 1.0;
                    let p = p_dist.sample(&mut rng);
                    particle_array_x[jp] = particle_array_x[p];
                    particle_array_y[jp] = particle_array_y[p];
                    particle_array_z[jp] = particle_array_z[p];
                }
            } else {
                let p = p_dist.sample(&mut rng);
                particle_array_x[jp] = particle_array_x[p];
                particle_array_y[jp] = particle_array_y[p];
                particle_array_z[jp] = particle_array_z[p];
            }
        }

        if jt % CHECK_INTERVAL == 0 {
            let tau = jt as f64 * DTAU;
            tau_values.push(tau);
            if N_PARTICLES_F > evaporated_particles {
                energy = (N_PARTICLES_F/(N_PARTICLES_F - evaporated_particles)).ln()/(tau - tau_current)
            }
            tau_current = tau;
            evaporated_particles = 0.0;
            energy_values.push(energy)
        }
    }

    println!("Main loop done");

    let energy_correct = mean_last_half(&energy_values);
    let n_values = energy_values.len();

    let flnm = format!("schr_MC_3D_energy-{}-{}-{}-{}", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX);
    let title = format!("{} particles, {} time steps, {}x{}x{} grid points.", N_PARTICLES, N_TIME_STEPS, NX, NY, NZ);
    let plot_par = plot::PlotPar::new(
        "Time, 1/eV", 
        "Energy, eV", 
        &title, 
        &flnm,
        vec![format!("Energy average"), format!("Energy vs time")],
    );
    plot::line_and_scatter_plot(
        &tau_values, 
        &vec![energy_correct; n_values], 
        &tau_values, 
        &energy_values, 
        &plot_par
    );

    save_columns_to_file(
        &vec![
            tau_values, energy_values, vec![energy_correct; n_values],
        ], "results", 
        &format!("schr_MC_1D_energy-{}-{}-{}-{}.dat", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX),
    );

    // Full denstiy distribution

    for jp in 0..N_PARTICLES {
        if particle_array_x[jp] < NX && particle_array_y[jp] < NY && particle_array_z[jp] < NZ {
            particle_distances[jp] = ((particle_array_x[jp] as f64).powi(2) + (particle_array_y[jp] as f64).powi(2) + (particle_array_z[jp] as f64).powi(2)).sqrt().floor() as usize;
            particle_distribution_xyz[index_3d(particle_array_x[jp],particle_array_y[jp],particle_array_z[jp])] += 1.0
        }
    }

    // 1D distribution by summation for all z and rotation around z

    let mut sum = 0.0;
    let mut square;
    for jr in 0..NX {
        let r = jr as f64 * DX - LX/2.0;
        for jphi in 0..N_ANGLES_ROTATION_PLOT {
            let jx = ((r*cosphin[jphi] + LX/2.0)/DX) as usize;
            let jy = ((r*sinphin[jphi] + LX/2.0)/DX) as usize;
            if jx < NX && jy < NY {
                for jz in 0..NZ {
                    square = particle_distribution_xyz[index_3d(jx, jy, jz)]*particle_distribution_xyz[index_3d(jx, jy, jz)];
                    particle_distribution_x[jr] += square;
                    sum += square;
                }
            }
        }
    }
    for jx in 0..NX {
        particle_distribution_x[jx] = particle_distribution_x[jx]/sum/DX + energy_correct
    }

    // 2D distribution xy by summation for all z

    for jx in 0..NX {
        for jy in 0..NY {
            for jz in 0..NZ {
                square = particle_distribution_xyz[index_3d(jx, jy, jz)]*particle_distribution_xyz[index_3d(jx, jy, jz)];
                particle_distribution_xy[index_xy(jx, jy)] += square;
            }
        }
        for jz in 0..NZ {
            for jr in 0..NX {
                let r = jr as f64 * DX - LX/2.0;
                for jphi in 0..N_ANGLES_ROTATION_PLOT {
                    let jx = ((r*cosphin[jphi] + LX/2.0)/DX) as usize;
                    let jy = ((r*sinphin[jphi] + LX/2.0)/DX) as usize;
                    if jx < NX && jy < NY {
                        square = particle_distribution_xyz[index_3d(jx, jy, jz)]*particle_distribution_xyz[index_3d(jx, jy, jz)];
                        particle_distribution_xz[index_xz(jr, NZ-jz-1)] += square;
                    }
                }
            }
        }
    }

    let flnm = format!("schr_MC_3D_density1D-{}-{}-{}-{}", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX);
    let title = format!("{} particles, {} time steps, {}x{}x{} grid points.", N_PARTICLES, N_TIME_STEPS, NX, NY, NZ);
    let plot_par = plot::PlotPar::new(
        "x, nm", 
        "Energy, eV", 
        &title, 
        &flnm,
        vec![format!("U(x)"), format!("Probability density")],
    );
    plot::line_and_scatter_plot(
        &xn, 
        &un1d, 
        &xn, 
        &particle_distribution_x, 
        &plot_par
    );
    save_columns_to_file(
        &vec![
            xn, un1d, particle_distribution_x,
        ], "results", 
        &format!("schr_MC_1D_density-{}-{}-{}-{}.dat", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX),
    );

    draw::plot_2d(&particle_distribution_xy, NX, NY, 2000, 2000, &format!("schr_MC_3D_density_xy-{}-{}-{}-{}.png", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX), "Turbo");
    draw::plot_2d(&particle_distribution_xz, NX, NZ, 2000, (2000 * NZ) / NX, &format!("schr_MC_3D_density_xz-{}-{}-{}-{}.png", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX), "Turbo");
}
