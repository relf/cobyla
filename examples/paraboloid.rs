use argmin::core::observers::{ObserverMode, SlogLogger};
use argmin::core::{CostFunction, Error, Executor};
use cobyla::{minimize, CobylaSolver, Func};

/// Problem cost function
fn paraboloid(x: &[f64], _data: &mut ()) -> f64 {
    10. * (x[0] + 1.).powf(2.) + x[1].powf(2.)
}

fn nlopt_paraboloid(x: &[f64], _g: Option<&mut [f64]>, _user_data: &mut ()) -> f64 {
    paraboloid(x, &mut ())
}

/// Problem Definition: minimize paraboloid(x) subject to x0 >= 0
struct ParaboloidProblem;
impl CostFunction for ParaboloidProblem {
    type Param = Vec<f64>;
    type Output = Vec<f64>;

    fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> {
        let fx = paraboloid(x, &mut ());
        Ok(vec![fx, x[0]])
    }
}

fn main() {
    let mut x = vec![1., 1.];

    println!("*** Solve paraboloid problem using nlopt_cobyla");
    let mut cons: Vec<&dyn Func<()>> = vec![];
    let cstr1 = |x: &[f64], _g: Option<&mut [f64]>, _u: &mut ()| x[0];
    cons.push(&cstr1);

    let (status, x_opt) = minimize(
        nlopt_paraboloid,
        &mut x,
        &cons,
        (),
        0.5,
        1e-4,
        200,
        0,
        (-10., 10.),
    );

    // For status meaning see cobyla/nlopt/nlopt.h
    println!("status = {}", status);
    println!("x = {:?}\n\n", x_opt);

    println!(
        "*** Solve paraboloid problem using Cobyla argmin solver implemented on top of fmin_cobyla impl"
    );
    let problem = ParaboloidProblem;
    let solver = CobylaSolver::new(vec![1., 1.]);

    let res = Executor::new(problem, solver)
        .configure(|state| state.max_iters(100).iprint(0))
        .add_observer(SlogLogger::term(), ObserverMode::Always)
        .run()
        .unwrap();

    // Wait a second (lets the logger flush everything before printing again)
    std::thread::sleep(std::time::Duration::from_secs(1));
    println!("*** Result argmin solver impl ***");
    println!("Result:\n{}", res);
}
