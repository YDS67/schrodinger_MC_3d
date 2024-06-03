// Saving data to files
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::io::BufRead;
use std::io::BufReader;

pub fn save_columns_to_file(vecs: &Vec<Vec<f64>>, dir_name: &str, file_name: &str) {
    let try_create_dir = std::fs::create_dir(dir_name);
    match try_create_dir {
        Ok(create_dir) => {create_dir},
        Err(_) => {},
    };
    let file_path: PathBuf = [dir_name, file_name].iter().collect();
    let mut my_file = File::create(file_path).expect("Error creating file");
    let vec_num = vecs.len();
    let mut el_nums = Vec::new();
    for i in 0..vec_num {
        el_nums.push(vecs[i].len())
    }
    let num = el_nums.iter().max().unwrap().to_owned();
    for j in 0..num {
        for i in 0..vec_num {
            if el_nums[i] > j {
                write!(my_file, "{:.6} ", vecs[i][j]).expect(&format!("Can't write line {} to file", j));
            } else {
                write!(my_file, "{} ", "NAN").expect(&format!("Can't write line {} to file", j));
            }
        }
        write!(my_file, "{}", "\n").expect(&format!("Can't write line {} to file", j));
    }
}

pub fn read_palette(dir_name: &str, fl_name: &str) -> Vec<Vec<u8>> {

    let file_path: PathBuf = [dir_name, fl_name].iter().collect();
    let file = File::open(file_path).expect("Error opening file");
    let reader = BufReader::new(file);
    let mut file_contents: Vec<String> = Vec::new();

    for line in reader.lines() {
        file_contents.push(line.unwrap());
    }

    let mut data_rows: Vec<Vec<u8>> = Vec::new();

    for line in file_contents {
        let row: Vec<u8> = line
            .split(|c| c == ' ' || c == '\t')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().unwrap())
            .collect();

        data_rows.push(row)
    }

    data_rows
}
