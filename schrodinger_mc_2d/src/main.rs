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
const L: f64 = 80.0;
const X1: f64 = L/2.0;
const X2: f64 = L/2.0;
const Y1: f64 = L/2.0;
const Y2: f64 = L/2.0;
const U0: f64 = 0.25;
const DX: f64 = 0.5;
const DX2: f64 = DX*DX;
const NX: usize = (L/DX) as usize;
const D: f64 = H2M0/ME/2.0;
const DTAU: f64 = DX2/4.0/D; // units 1/eV
const FINAL_TIME: f64 = 100.0; // units 1/eV
const N_PARTICLES: usize = 500000;
const N_PARTICLES_F: f64 = N_PARTICLES as f64;
const N_TIME_STEPS: usize = (FINAL_TIME/DTAU) as usize;
const N_VALUES_CHECK: usize = 100;
const CHECK_INTERVAL: usize = N_TIME_STEPS/N_VALUES_CHECK;

fn potential_shape(x: f64, y: f64) -> f64 {
    let s1 = (x-X1).powi(2)/A1.powi(2)+(y-Y1).powi(2)/B1.powi(2);
    let s2 = (x-X2).powi(2)/A2.powi(2)+(y-Y2).powi(2)/B2.powi(2);
    (1.0 + s1)*(-s1).exp() - (1.0 + s2)*(-s2).exp()
}

fn potential(x: f64, y: f64, a: f64) -> f64 {
    U0*(1.0-potential_shape(x, y)/a)
}

fn index(jx: usize, jy: usize) -> usize {
    NX*jx+jy
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
    let amplitude = potential_shape(radius+X1, Y1);
    println!("QR amplitude = {} U0", amplitude);
    let mut rng = rand::thread_rng();
    let x_dist = Uniform::from(0..NX);
    let p_dist = Uniform::from(0..N_PARTICLES);
    let d_dist = Uniform::from(0..4);
    let xn: Vec<f64> = (0..NX).map(|n| n as f64 * DX).collect();
    let yn: Vec<f64> = (0..NX).map(|n| n as f64 * DX).collect();
    let un1d: Vec<f64> = xn.clone().into_iter().map(|x| potential(x, Y1, amplitude)).collect();
    let mut un2d: Vec<f64> = vec![0.0; NX*NX];
    for jx in 0..NX {
        for jy in 0..NX {
            un2d[index(jx,jy)] = potential(xn[jx], yn[jy], amplitude)
        }
    }
    let p_evap: Vec<f64> = un2d.clone().into_iter().map(|u| u*DTAU).collect();
    let mut particle_distribution: Vec<f64> = vec![0.0; NX*NX];
    let mut particle_distribution_1d: Vec<f64> = vec![0.0; NX];
    
    let mut energy: f64 = 0.0;
    let mut energy_values: Vec<f64> = Vec::new();
    let mut tau_values: Vec<f64> = Vec::new();
    let mut evaporated_particles: f64 = 0.0;
    let mut tau_current: f64 = 0.0;
    let mut particle_array_x = vec![1; N_PARTICLES];
    let mut particle_array_y = vec![1; N_PARTICLES];
    let mut particle_distances = vec![1; N_PARTICLES];
    let mut r1: f64;
    let mut p1: usize;

    for jp in 0..N_PARTICLES {
        particle_array_x[jp] = x_dist.sample(&mut rng);
        particle_array_y[jp] = x_dist.sample(&mut rng);
    }

    for jt in 0..N_TIME_STEPS {
        for jp in 0..N_PARTICLES {
            p1 = d_dist.sample(&mut rng);
            match p1 {
                0 => {particle_array_x[jp] -= 1},
                1 => {particle_array_x[jp] += 1},
                2 => {particle_array_y[jp] -= 1},
                3 => {particle_array_y[jp] += 1},
                _ => {},
            }
            if particle_array_x[jp] < NX && particle_array_y[jp] < NX {
                r1 = rng.gen();
                if r1 < p_evap[index(particle_array_x[jp],particle_array_y[jp])] {
                    evaporated_particles += 1.0;
                    let p = p_dist.sample(&mut rng);
                    particle_array_x[jp] = particle_array_x[p];
                    particle_array_y[jp] = particle_array_y[p];
                }
            } else {
                let p = p_dist.sample(&mut rng);
                particle_array_x[jp] = particle_array_x[p];
                particle_array_y[jp] = particle_array_y[p];
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

    let energy_correct = mean_last_half(&energy_values);
    let n_values = energy_values.len();

    for jp in 0..N_PARTICLES {
        if particle_array_x[jp] < NX && particle_array_y[jp] < NX {
            particle_distances[jp] = ((particle_array_x[jp] as f64).powi(2) + (particle_array_y[jp] as f64).powi(2)).sqrt().floor() as usize;
            particle_distribution[index(particle_array_x[jp],particle_array_y[jp])] += 1.0
        }
    }
    let mut sum = 0.0;
    let mut square;
    for jx in 0..NX {
        square = particle_distribution[index(jx, NX/2)]*particle_distribution[index(jx, NX/2)];
        particle_distribution_1d[jx] = square;
        sum += square;
        for jy in 0..NX {
            square = particle_distribution[index(jx, jy)]*particle_distribution[index(jx, jy)];
            particle_distribution[index(jx, jy)] = square;
        }
    }
    for jx in 0..NX {
        particle_distribution_1d[jx] = particle_distribution_1d[jx]/sum/DX + energy_correct
    }


    let flnm = format!("schr_MC_2D_energy-{}-{}-{}-{}", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX);
    let title = format!("{} particles, {} time steps, {} values saved, {} grid points.", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX);
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

    let flnm = format!("schr_MC_2D_density1D-{}-{}-{}-{}", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX);
    let title = format!("{} particles, {} time steps, {} values saved, {} grid points.", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX);
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
        &particle_distribution_1d, 
        &plot_par
    );
    save_columns_to_file(
        &vec![
            xn, un1d, particle_distribution_1d,
        ], "results", 
        &format!("schr_MC_1D_density-{}-{}-{}-{}.dat", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX),
    );
    save_columns_to_file(
        &vec![
            tau_values, energy_values, vec![energy_correct; n_values],
        ], "results", 
        &format!("schr_MC_1D_energy-{}-{}-{}-{}.dat", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX),
    );

    draw::plot_2d(&un2d, NX, NX, 2000, 2000, "potential_plot_2d.png");

    draw::plot_2d(&particle_distribution, NX, NX, 2000, 2000, &format!("schr_MC_2D_density-{}-{}-{}-{}.png", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX));
}
