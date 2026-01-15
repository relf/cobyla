use cobyla::{Func, RhoBeg, StopTols, minimize};

/// Problem cost function
fn paraboloid(x: &[f64], _data: &mut ()) -> f64 {
    10. * (x[0] + 1.).powf(2.) + x[1].powf(2.)
}

fn main() {
    println!("*** Solve paraboloid problem using nlopt_cobyla");

    // Initial guess
    let xinit = vec![1., 1.];

    // Define a constraint: x0 > 0
    let mut cons: Vec<&dyn Func<()>> = vec![];
    let cstr1 = |x: &[f64], _u: &mut ()| x[0];
    cons.push(&cstr1);

    // Define a stopping criterion on objective function change
    let stop_tol = StopTols {
        ftol_rel: 1e-4,
        ..StopTols::default()
    };

    match minimize(
        paraboloid,
        &xinit,
        &[(-10., 10.), (-10., 10.)],
        &cons,
        (),
        200,
        RhoBeg::All(0.5),
        Some(stop_tol),
    ) {
        Ok((status, x_opt, y_opt)) => {
            println!("status = {:?}", status);
            println!("x_opt = {:?}", x_opt);
            println!("y_opt = {}", y_opt);
        }
        Err((e, _, _)) => println!("Optim error: {:?}", e),
    }
}
