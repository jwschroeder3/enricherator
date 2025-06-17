use nlopt::{Nlopt, Algorithm::Lbfgs, Target::Maximize, FailState};
//use statrs::function::gamma;
use parse_csv_to_blob::BoxedResult;
use std::error::Error;
use plotters::prelude::*;

/// Function to plot the fitted relationship
pub fn plot_relationships(
        path: &str,
        geno_name: &str,
        non_sub: Vec<f64>,
        beta: Vec<f64>,
        gamma: Vec<f64>,
        a: f64,
        b: f64,
        c: f64,
) -> BoxedResult<()> {

    let y_min = non_sub.iter().cloned().reduce(f64::min).unwrap() - 0.1;
    let x_min = beta.iter().cloned().reduce(f64::min).unwrap() - 0.1;
    let y_max = non_sub.iter().cloned().reduce(f64::max).unwrap() + 0.1;
    let x_max = beta.iter().cloned().reduce(f64::max).unwrap() + 0.1;

    // get 100 values between x_min and x_max
    let step = x_max / (100 - 1) as f64;
    let line_x: Vec<f64> = (0..100).map(|i| 0.0 + step * i as f64).collect();
    let line_y: Vec<f64> = line_x.iter().map(|x| cubic_interpolation(x, a, b, c)).collect();
    let mut line: Vec<(f64, f64)> = vec![(x_min, 0.0), (0.0,0.0)];
    let mut curve: Vec<(f64, f64)> = line_x.iter().zip(line_y).map(|(x,y)| (*x as f64,y as f64)).collect();
    line.append(&mut curve);
    let points: Vec<(f64,f64)> = beta.iter().zip(non_sub).map(|(x,y)| (*x,y)).collect();
    let gamma_points: Vec<(f64,f64)> = beta.iter().zip(gamma).map(|(x,y)| (*x,y)).collect();

    let fname = format!("{}/{}_non_sub_vs_beta.png", path, geno_name);
    let root = BitMapBackend::new(&fname, (640, 480)).into_drawing_area();
    let _ = root.fill(&WHITE);
    let root = root.margin(10, 10, 10, 10);
    // After this point, we should be able to construct a chart context
    let mut chart = ChartBuilder::on(&root)
        // Set the caption of the chart
        .caption(format!("{} extract vs control relationship", geno_name), ("sans-serif", 40).into_font())
        // Set the size of the label region
        .x_label_area_size(20)
        .y_label_area_size(40)
        // Finally attach a coordinate on the drawing area and make a chart context
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    // Then we can draw a mesh
    chart
        .configure_mesh()
        // We can customize the maximum number of labels allowed for each axis
        .x_labels(5)
        .y_labels(5)
        // We can also change the format of the label text
        .y_label_formatter(&|y| format!("{:.1}", y))
        .draw()?;

    // We can draw point series
    chart.draw_series(PointSeries::of_element(
        points,
        5,
        &BLACK,
        &|c, s, st| {
            return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
            + Circle::new((0,0),s,st.filled()); // At this point, the new pixel coordinate is established
            //+ Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
        },
    ))?;

    // And we can draw something in the drawing area
    chart.draw_series(LineSeries::new(line, &RED,))?;
    root.present()?;
    println!("Plot {fname} written to disk");

    let fname = format!("{}/{}_gamma_vs_beta.png", path, geno_name);
    let gamma_root = BitMapBackend::new(&fname, (640, 480)).into_drawing_area();
    let _ = gamma_root.fill(&WHITE);
    let gamma_root = gamma_root.margin(10, 10, 10, 10);
    // After this point, we should be able to construct a chart context
    let mut chart = ChartBuilder::on(&gamma_root)
        // Set the caption of the chart
        .caption(format!("{} subtracted (Gamma) vs control relationship", geno_name), ("sans-serif", 40).into_font())
        // Set the size of the label region
        .x_label_area_size(20)
        .y_label_area_size(40)
        // Finally attach a coordinate on the drawing area and make a chart context
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    // Then we can draw a mesh
    chart
        .configure_mesh()
        // We can customize the maximum number of labels allowed for each axis
        .x_labels(5)
        .y_labels(5)
        // We can also change the format of the label text
        .y_label_formatter(&|y| format!("{:.1}", y))
        .draw()?;

    // We can draw point series
    chart.draw_series(PointSeries::of_element(
        gamma_points,
        5,
        &BLACK,
        &|c, s, st| {
            return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
            + Circle::new((0,0),s,st.filled()); // At this point, the new pixel coordinate is established
            //+ Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
        },
    ))?;

    gamma_root.present()?;
    println!("Plot {fname} written to disk");

    Ok(())
}

/// Predicts the y value based on x and the three cubic parameters
fn cubic_interpolation(x: &f64, a: f64, b: f64, c: f64) -> f64 {
    a*x + b*x.powf(2.0) + c*x.powf(3.0)
}

/// Simplified student t ln(pdf) for nu = 2.0, mu = 0.0
fn student_t_ln_pdf(x: &f64, sd: &f64) -> f64 {
    // When df set to 2, and mean to 0, the simple ln(p(x|sd)) is,
    let fact1 = 2.0 * sd * 2.0_f64.sqrt();
    let fact2 = 1.0 + 0.5 * (x / sd).powf(2.0);
    let ln_pdf = (1.0 / (fact1 * fact2.powf(1.5))).ln();
    ln_pdf
}

/// Objective function to be minimized by nlopt
fn objective(x: &[f64], _grad: Option<&mut [f64]>, user_data: &mut Vec<f64>) -> f64 {
    let n = user_data[0] as usize;
    let independ = &user_data[1..n+1];
    let depend = &user_data[n+1..n*2+1];
    let data: Vec<(&f64, &f64)> = independ.iter().zip(depend).collect();
    // get the residuals
    let residuals: Vec<f64> = data.iter().map(|(indep,dep)| {
        let pred = cubic_interpolation(indep, x[0], x[1], x[2]);
        let residual = pred - *dep;
        residual
    }).collect();
    // calculate the log-likelihood
    let ll: f64 = residuals.iter().map(|r| student_t_ln_pdf(r, &x[3])).sum();
    ll
}

/// Performs the fit of a cubic relationship between y and x.
/// Uses L-BFGS to do MLE, where likelihood is student-t with 2 df.
pub fn cubic_fit(
        x: &mut Vec<f64>,
        y: &mut Vec<f64>,
) -> Result<[f64; 4], Box<FailState>> {
    let n = x.len() as f64;
    let mut data: Vec<f64> = Vec::new();
    data.push(n);
    data.append(x);
    data.append(y);
    // 4 dims to optimize, one for each of the three cubic params, and one
    // for the sd of the student t distribution. We're fixing df to 2 and mean
    // to 0.
    let mut optimizer = Nlopt::new(
        Lbfgs,
        4,
        objective,
        Maximize,
        data,
    );
    // contrain linear component > 0 and sd > 0
    let _ = optimizer.set_lower_bounds(&[0.0,f64::NEG_INFINITY,f64::NEG_INFINITY,0.0]);
    // instantiate start values for optimizer
    let mut opt_vals = [1.0, 0.0, 0.0, 1.0];
    // actually do the optimization
    let opt_res = optimizer.optimize(&mut opt_vals);
    // handle the result
    match opt_res {
        Ok((succ_state, succ_val)) => {
            println!("succ_state: {:?}", &succ_state);
            println!("succ_val: {:?}", &succ_val);
            return Ok(opt_vals);
        },
        Err((fail_state, fail_val)) => {
            println!("fail_state: {:?}", &fail_state);
            println!("fail_val: {:?}", &fail_val);
            return Err(Box::new(fail_state));
        },
    };
}
