use crate::file;
use image::{DynamicImage, ImageBuffer, RgbImage};

pub fn plot_2d(map: &Vec<f64>, nx: usize, ny: usize, width: usize, height: usize, flnm: &str) {
    let mut img: RgbImage = ImageBuffer::new(nx as u32, ny as u32);
    let mut idx: usize;
    let maxmin = max_min(&map);

    let palette: Vec<Vec<u8>> = file::read_palette("palettes", "Turbo.dat");

    for jx in 0..nx {
        for jy in 0..ny {
            idx = jx * ny + jy;
            let u = map[idx];
            let s = scale(u, maxmin[1], maxmin[0]);
            let col = &palette[s];
            img.put_pixel(jx as u32, jy as u32, image::Rgb([col[0],col[1],col[2]]));
        }
    }

    let mut img2 = DynamicImage::ImageRgb8(img);
    img2 = img2.resize(
        width as u32,
        height as u32,
        image::imageops::FilterType::Nearest,
    );

    img2.save(&flnm).unwrap();
}

fn scale(f: f64, f_min: f64, f_max: f64) -> usize {
    let z = ((f - f_min) / (f_max - f_min) * 255.0) as usize;
    z
}

fn max_min(v: &Vec<f64>) -> [f64; 2] {
    let mut umax = 0.0;
    let mut umin = 0.0;
    let mut u;
    for j in 0..v.len() {
        u = v[j];
        if umax < u {umax = u}
        if umin > u {umin = u}
    }
    [umax, umin]
}