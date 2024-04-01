use std::error::Error;
use std::env;
use std::collections::VecDeque;
use std::io::Write;
use std::fs::File;
use std::path::Path;

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::*;

    //#[test]
    //fn test_build_sparse_mat() {

    //    let pos_num = 15;
    //    let mut ctg_ends = VecDeque<i64>::new();
    //    ctg_ends.push_back(11);
    //    ctg_ends.push_back(14);
    //    let mut values: Vec<(Vec<f64>, Vec<i64>)> = Vec::with_capacity(pos_num);
    //    for i in 0..pos_num as usize {
    //        values.push((Vec::with_capacity(500), Vec::with_capacity(500)));
    //    }

    //    let mut ctg_end = ctg_ends.pop_front().unwrap();
    //    let mut ctg_start = 0;
    //    let mut genome_top: i64 = 0;
    //    let mut genome_bottom: i64 = 1000000;

    //    for i in 0..sub_L {
    //        let genome_center = i*sub_dist + sub_dist/2;
    //        genome_top = genome_center - half_kern_size;
    //        genome_bottom = genome_center + half_kern_size;
    //        let mut kern_top = 0;
    //        let mut kern_bottom = kern_size-1;

    //    }

    //    assert_eq!(0.0, 1.0);
    //}

    #[test]
    fn test_exp_kern() {
        let kern = set_up_exponential_kernel(6.0);
        let targets = vec![0.004374576245646989, 0.006105212962911393, 0.0085205110688358, 0.0118913311, 0.0165956895, 0.0231611504, 0.0323239893, 0.0451117611, 0.0629585343, 0.0878657127, 0.1226264804, 0.1711390397, 0.2388437702, 0.3333333333, 0.2388437702, 0.1711390397, 0.1226264804, 0.0878657127, 0.0629585343, 0.0451117611, 0.0323239893, 0.0231611504, 0.0165956895, 0.0118913311, 0.0085205110688358, 0.006105212962911393, 0.004374576245646989];
        assert_eq!(kern.len(), targets.len());
        for (i,target) in targets.iter().enumerate() {
            let val = kern[i];
            assert_approx_eq!(&val, target, 0.0000000001);
        }
    }
}

fn parse_args(args: &Vec<String>) -> Result<(i64,i64,i64,f64,VecDeque<i64>), Box<dyn Error>> {
    // let L=number of positions
    let L: i64 = args[1].parse::<i64>()?;
    let sub_L: i64 = args[2].parse::<i64>()?;
    let C: i64 = args[3].parse::<i64>()?;
    // K is the kernel width in at the resolution of the original analysis. That's set in the R
    // code.
    let K: f64 = args[4].parse::<f64>()?;
    let ctg_ends: VecDeque<i64> = args[5]
        .split(',')
        .map(|x| x.parse().unwrap())
        .collect();
    Ok((L,sub_L,C,K,ctg_ends))
}

fn normal_pdf(sd: f64, x: f64) -> f64 {
    let fac1 = 1.0 / (2.0*sd.powi(2) * std::f64::consts::PI).sqrt();
    let fac2 = ((-1.0 * (x).powi(2)) / (2.0 * sd.powi(2))).exp();
    fac1 * fac2
}

fn exponential_pdf(lambda: f64, x: f64) -> f64 {
    lambda * (-lambda * x).exp()
}

fn set_up_gaussian_kernel(kern_width: f64) -> Vec<f64> {
    let sd: f64 = kern_width / 2.0;
    let kern: Vec<f64> = (-5000..5000).filter_map(|x| {
            let y = normal_pdf(sd, x as f64);
            if y >= 0.01 { 
                Some(y)
            } else {
                None
            }
        }).collect();
    kern
}

fn set_up_exponential_kernel(kern_width: f64) -> Vec<f64> {
    let rate: f64 = 1.0 / (kern_width / 2.0);
    let mut kern: Vec<f64> = (0..5000).filter_map(|x| {
            let y = exponential_pdf(rate, x as f64);
            if y >= rate * 0.01 { 
                Some(y)
            } else {
                None
            }
        }).collect();
    kern.reverse();
    let mut right_side: Vec<f64> = (1..5000).filter_map(|x| {
            let y = exponential_pdf(rate, x as f64);
            if y >= rate * 0.01 { 
                Some(y)
            } else {
                None
            }
        }).collect();
    kern.append(&mut right_side);
    kern
}

fn place_value_in_vals(
        vals: &mut Vec<(Vec<f64>, Vec<i64>)>,
        row: usize,
        col: i64,
        value: f64,
) {
    vals[row].0.push(value);
    vals[row].1.push(col+1);
}

fn write_sparse_vals(
        vals: &mut Vec<(Vec<f64>, Vec<i64>)>,
        w_file: &mut File,
        v_file: &mut File,
        u_file: &mut File,
) -> Result<(), Box<dyn Error>> {

    let mut length: u64 = 1;
    write!(u_file, "1\n")?;

    for i in 0..vals.len() {
        let mut l_info = &mut vals[i];
        length += l_info.0.len() as u64;
        for val in &l_info.0 {
            write!(w_file, "{}\n", val)?;
        }
        for val in &l_info.1 {
            write!(v_file, "{}\n", val)?;
        }
        write!(u_file, "{}\n", length)?;
    }
    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let (pos_num, sub_L, sub_dist, kern_width, mut ctg_ends) = parse_args(&args).unwrap();

    let w_outname = Path::new(&args[6]);
    let v_outname = Path::new(&args[7]);
    let u_outname = Path::new(&args[8]);

    let kernel: Vec<f64> = set_up_exponential_kernel(kern_width);
    let kern_size = kernel.len() as i64;
    let half_kern_size = (kern_size / 2) as i64;
    let mut values: Vec<(Vec<f64>, Vec<i64>)> = Vec::with_capacity(
        pos_num.try_into().unwrap()
    );
    for i in 0..pos_num as usize {
        values.push((Vec::with_capacity(500), Vec::with_capacity(500)));
    }
    //println!("{}", values.len());
    //println!("{}", values.capacity());
    let mut ctg_end = ctg_ends.pop_front().unwrap();
    let mut ctg_start = 0;
    let mut genome_top: i64 = 0;
    let mut genome_bottom: i64 = 1000000;
    for i in 0..sub_L {
        let genome_center = i*sub_dist + sub_dist/2;
        genome_top = genome_center - half_kern_size;
        genome_bottom = genome_center + half_kern_size;
        let mut kern_top = 0;
        let mut kern_bottom = kern_size-1;

        if genome_center > ctg_end {
            ctg_start = ctg_end + 1;
            let opt = ctg_ends.pop_front();
            if opt.is_none() {
                println!("None in ctg_ends: {:?}", ctg_ends);
                println!("genome_center: {}", genome_center);
            }
            ctg_end = opt.unwrap()
        }

        if genome_top < ctg_start {
            kern_top = ctg_start - genome_top;
            genome_top = ctg_start;
        }

        if genome_bottom > ctg_end {
            let overreach = genome_bottom - ctg_end;
            genome_bottom = ctg_end;
            kern_bottom = kern_size-1-overreach;
        }

        //println!("Genome top: {}", genome_top);
        //println!("Genome bottom: {}", genome_bottom);
        //println!("Kernel top: {}", kern_top);
        //println!("Kernel bottom: {}", kern_bottom);
        let rows: Vec<usize> = (genome_top as usize..genome_bottom as usize).collect();
        //println!("{:?}", rows);
        //println!("{:?}", rows.len());
        let mut kern_values = &kernel[kern_top as usize..kern_bottom as usize];
        //println!("{:?}", kern_values);
        //println!("{:?}", kern_values.len());
        if kern_values.len() > rows.len() {
            kern_values = &kern_values[0..rows.len()-1];
        }

        for j in 0..kern_values.len() {
            place_value_in_vals(
                &mut values,
                rows[j],
                i,
                kern_values[j],
            );
        }
    }

    // fill in final spots with 1.0
    let remaining = pos_num - genome_bottom;
    if remaining > 0 {
        for j in (genome_bottom + 1) as usize..pos_num as usize {
            place_value_in_vals(
                &mut values,
                j,
                sub_L,
                1.0,
            );
        }
    }

    //println!("{:?}", kernel);
    //println!("Kernel length: {}", kern_size);
    //println!("Contig ends: {:?}", ctg_ends);
    //println!("Values 0: {:?}", &values[0].0);
    for l in 0..pos_num as usize {
        let these_vals = &values[l].0;
        let summed: f64 = these_vals.iter().sum();
        values[l].0 = these_vals.iter().map(|x| x / summed).collect();
    }

    let mut w_f = File::create(w_outname).unwrap();
    let mut v_f = File::create(v_outname).unwrap();
    let mut u_f = File::create(u_outname).unwrap();
    let _ = write_sparse_vals(&mut values, &mut w_f, &mut v_f, &mut u_f);
}

