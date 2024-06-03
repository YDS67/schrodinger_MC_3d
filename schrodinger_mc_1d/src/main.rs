use file::save_columns_to_file;
use rand::{self, Rng, rngs::ThreadRng};

mod file;
mod plot;

const H2M0: f64 = 0.0762;
const ME: f64 = 0.067;
const A: f64 = 4.0;
const L: f64 = 60.0;
const X1: f64 = L/3.0;
const X2: f64 = 2.0*L/3.0;
const U0: f64 = 0.25;
const DX: f64 = 0.1;
const DX2: f64 = DX*DX;
const NX: usize = (L/DX) as usize;
const D: f64 = H2M0/ME/2.0;
const DTAU: f64 = DX2/2.0/D; // units 1/eV
const FINAL_TIME: f64 = 50.0; // units 1/eV
const N_PARTICLES: usize = 100000;
const N_PARTICLES_F: f64 = N_PARTICLES as f64;
const N_TIME_STEPS: usize = (FINAL_TIME/DTAU) as usize;
const N_VALUES_CHECK: usize = 200;
const CHECK_INTERVAL: usize = N_TIME_STEPS/N_VALUES_CHECK;

fn potential(x: f64) -> f64 {
    U0*(1.0 - (-(x-X1).powi(2)/A.powi(2)).exp() - (-(x-X2).powi(2)/A.powi(2)).exp())
}

fn sample(v: &Vec<usize>, rng: &mut ThreadRng) -> usize {
    let l = v.len() as f64;
    let r1: f64 = rng.gen();
    let e = (r1*(l-1.0)).floor() as usize;
    v[e]
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
    let mut rng = rand::thread_rng();
    let indices: Vec<usize> = (0..NX).collect();
    let xn: Vec<f64> = (0..NX).map(|n| n as f64 * DX).collect();
    let un: Vec<f64> = xn.clone().into_iter().map(|x| potential(x)).collect();
    let p_evap: Vec<f64> = un.clone().into_iter().map(|u| u*DTAU).collect();
    let mut energy: f64 = 0.0;
    let mut energy_values: Vec<f64> = Vec::new();
    let mut tau_values: Vec<f64> = Vec::new();
    let mut evaporated_particles: f64 = 0.0;
    let mut tau_current: f64 = 0.0;
    let mut particle_array = vec![1; N_PARTICLES];
    let mut r1: f64;

    for jp in 0..N_PARTICLES {
        particle_array[jp] = sample(&indices, &mut rng);
    }

    for jt in 0..N_TIME_STEPS {
        for jp in 0..N_PARTICLES {
            r1 = rng.gen();
            if r1 < 0.5 {
                particle_array[jp] = particle_array[jp] - 1
            } else {
                particle_array[jp] = particle_array[jp] + 1
            }
            if particle_array[jp] < NX {
                r1 = rng.gen();
                if r1 < p_evap[particle_array[jp]] {
                    evaporated_particles += 1.0;
                    particle_array[jp] = sample(&particle_array, &mut rng)
                }
            } else {
                particle_array[jp] = sample(&particle_array, &mut rng)
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

    let mut particle_distribution: Vec<f64> = vec![0.0; NX];

    for jp in 0..N_PARTICLES {
        if particle_array[jp] < NX {
            particle_distribution[particle_array[jp]] += 1.0
        }
    }
    let mut sum = 0.0;
    let mut square;
    for jx in 0..NX {
        square = particle_distribution[jx]*particle_distribution[jx];
        sum += square;
        particle_distribution[jx] = square;
    }
    for jx in 0..NX {
        particle_distribution[jx] = particle_distribution[jx]/sum/DX + energy_correct
    }


    let flnm = format!("schr_FD_MC_1D_energy-{}-{}-{}-{}", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX);
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

    let flnm = format!("schr_FD_MC_1D_density-{}-{}-{}-{}", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX);
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
        &un, 
        &xn, 
        &particle_distribution, 
        &plot_par
    );
    save_columns_to_file(
        &vec![
            xn, un, particle_distribution,
        ], "results", 
        &format!("schr_FD_MC_1D_density-{}-{}-{}-{}.dat", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX),
    );
    save_columns_to_file(
        &vec![
            tau_values, energy_values, vec![energy_correct; n_values],
        ], "results", 
        &format!("schr_FD_MC_1D_energy-{}-{}-{}-{}.dat", N_PARTICLES, N_TIME_STEPS, N_VALUES_CHECK, NX),
    );
}
