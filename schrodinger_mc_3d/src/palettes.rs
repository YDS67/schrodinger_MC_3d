const BNW: &[u8] = include_bytes!("../palettes/BnW.dat");
const CIVIDIS: &[u8] = include_bytes!("../palettes/Cividis.dat");
const INFERNO: &[u8] = include_bytes!("../palettes/Inferno.dat");
const LAJOLLA: &[u8] = include_bytes!("../palettes/Lajolla.dat");
const RDBU: &[u8] = include_bytes!("../palettes/RdBu.dat");
const TURBO: &[u8] = include_bytes!("../palettes/Turbo.dat");
const VIRIDIS: &[u8] = include_bytes!("../palettes/Viridis.dat");
const YLORRD: &[u8] = include_bytes!("../palettes/YlOrRd.dat");

pub fn load_palettes() -> Vec<Vec<Vec<u8>>> {

    let mut palettes = Vec::new();
    palettes.push(read_palette(BNW));
    palettes.push(read_palette(CIVIDIS));
    palettes.push(read_palette(INFERNO));
    palettes.push(read_palette(LAJOLLA));
    palettes.push(read_palette(RDBU));
    palettes.push(read_palette(TURBO));
    palettes.push(read_palette(VIRIDIS));
    palettes.push(read_palette(YLORRD));

    palettes
}

pub fn read_palette(file_content: &[u8]) -> Vec<Vec<u8>> {
    use std::io::BufRead;
    use std::io::BufReader;

    let reader = BufReader::new(file_content);
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