use cobyla::{Func, RhoBeg, StopTols, minimize, minimize_with_nevals};

/// Optional COBYLA implementation as an argmin solver: CobylaSolver
#[cfg(feature = "argmin")]
use argmin::core::{CostFunction, Error, Executor, observers::ObserverMode};
#[cfg(feature = "argmin")]
use argmin_observer_slog::SlogLogger;
#[cfg(feature = "argmin")]
use cobyla::CobylaSolver;

/// Problem cost function
fn paraboloid(x: &[f64], _data: &mut ()) -> f64 {
    10. * (x[0] + 1.).powf(2.) + x[1].powf(2.)
}

/// Problem Definition for CobylaSolver : minimize paraboloid(x) subject to x0 >= 0
#[cfg(feature = "argmin")]
struct ParaboloidProblem;
#[cfg(feature = "argmin")]
impl CostFunction for ParaboloidProblem {
    type Param = Vec<f64>;
    type Output = Vec<f64>;

    fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> {
        let fx = paraboloid(x, &mut ());
        Ok(vec![fx, x[0]])
    }
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

    // Example 1: Using minimize without function evaluation count
    match minimize(
        paraboloid,
        &xinit,
        &[(-10., 10.), (-10., 10.)],
        &cons,
        (),
        200,
        RhoBeg::All(0.5),
        Some(stop_tol.clone()),
    ) {
        Ok((status, x_opt, y_opt)) => {
            println!("status = {:?}", status);
            println!("x_opt = {:?}", x_opt);
            println!("y_opt = {}", y_opt);
        }
        Err((e, _, _)) => println!("Optim error: {:?}", e),
    }

    println!("\n*** With function evaluation count");
    
    // Example 2: Using minimize_with_nevals to get function evaluation count
    let (result, nfeval) = minimize_with_nevals(
        paraboloid,
        &xinit,
        &[(-10., 10.), (-10., 10.)],
        &cons,
        (),
        200,
        RhoBeg::All(0.5),
        Some(stop_tol),
    );
    match result {
        Ok((status, x_opt, y_opt)) => {
            println!("status = {:?}", status);
            println!("x_opt = {:?}", x_opt);
            println!("y_opt = {}", y_opt);
            println!("function evaluations = {}", nfeval);
        }
        Err((e, _, _)) => println!("Optim error: {:?}", e),
    }

    #[cfg(feature = "argmin")]
    {
        println!(
            "*** Solve paraboloid problem using Cobyla argmin solver implemented on top of fmin_cobyla impl"
        );
        let problem = ParaboloidProblem;
        let solver = CobylaSolver::new(vec![1., 1.]);

        let res = Executor::new(problem, solver)
            .timer(true)
            .configure(|state| state.max_iters(100).iprint(0))
            .add_observer(SlogLogger::term(), ObserverMode::Always)
            .run()
            .unwrap();

        // Wait a second (lets the logger flush everything before printing again)
        std::thread::sleep(std::time::Duration::from_secs(1));
        println!("*** Result argmin solver impl ***");
        println!("Result:\n{}", res);
    }
}
