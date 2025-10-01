#![doc = include_str!("../README.md")]

mod cobyla;
use nlopt_cobyla::nlopt_constraint;

mod nlopt_cobyla;
pub use crate::nlopt_cobyla::Func;

use crate::nlopt_cobyla::{
    cobyla_minimize,
    nlopt_constraint_raw_callback, // nlopt_eval_constraint,
    nlopt_function_raw_callback,
    nlopt_stopping,
    NLoptConstraintCfg,
    NLoptFunctionCfg,
};

mod cobyla_solver;
mod cobyla_state;
pub use crate::cobyla_solver::*;
pub use crate::cobyla_state::*;

use std::os::raw::c_void;

/// Failed termination status of the optimization process
#[derive(Debug, Clone, Copy)]
pub enum FailStatus {
    Failure,
    InvalidArgs,
    OutOfMemory,
    RoundoffLimited,
    ForcedStop,
    UnexpectedError,
}

/// Successful termination status of the optimization process
#[derive(Debug, Clone, Copy)]
pub enum SuccessStatus {
    Success,
    StopValReached,
    FtolReached,
    XtolReached,
    MaxEvalReached,
    MaxTimeReached,
}

/// Outcome when optimization process fails
type FailOutcome = (FailStatus, Vec<f64>, f64);
/// Outcome when optimization process succeeds
type SuccessOutcome = (SuccessStatus, Vec<f64>, f64);

/// Tolerances used as termination criteria.
/// For all, condition is disabled if value is not strictly positive.
/// ```rust
/// # use crate::cobyla::StopTols;
/// let stop_tol = StopTols {
///     ftol_rel: 1e-4,
///     xtol_abs: vec![1e-3; 3],   // size should be equal to x dim
///     ..StopTols::default()      // default stop conditions are disabled
/// };  
/// ```
#[derive(Debug, Clone, Default)]
pub struct StopTols {
    /// Relative tolerance on function value, algorithm stops when `func(x)` changes by less than `ftol_rel * func(x)`
    pub ftol_rel: f64,
    /// Absolute tolerance on function value, algorithm stops when `func(x)` change is less than `ftol_rel`
    pub ftol_abs: f64,
    /// Relative tolerance on optimization parameters, algorithm stops when all `x[i]` changes by less than `xtol_rel * x[i]`
    pub xtol_rel: f64,
    /// Relative tolerance on optimization parameters, algorithm stops when `x[i]` changes by less than `xtol_abs[i]`
    pub xtol_abs: Vec<f64>,
}

/// An enum for specifying the initial change of x which correspond to the `rhobeg`
/// argument of the original Powell's algorithm (hence the name)
pub enum RhoBeg {
    /// Used when all x components changes are specified with a single given value
    All(f64),
    /// Used to set the components with the given x-dim-sized vector
    Set(Vec<f64>),
}

/// Minimizes a function using the Constrained Optimization By Linear Approximation (COBYLA) method.
///
/// ## Arguments
///
/// * `func` - the function to minimize
/// * `xinit` - n-vector the initial guess
/// * `bounds` - x domain specified as a n-vector of tuple `(lower bound, upper bound)`  
/// * `cons` - slice of constraint function intended to be negative at the end
/// * `args` - user data pass to objective and constraint functions
/// * `maxeval` - maximum number of objective function evaluation
/// * `rhobeg`- initial changes to the x component
///     
/// ## Returns
///
/// The status of the optimization process, the argmin value and the objective function value
///
/// ## Panics
///
/// When some vector arguments like `bounds`, `xtol_abs` do not have the same size as `xinit`
///
/// ## Implementation note:
///
/// This implementation is a translation of [NLopt](https://github.com/stevengj/nlopt) 2.7.1
/// See also [NLopt SLSQP](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#slsqp) documentation.
///
/// ## Example
/// ```
/// # use approx::assert_abs_diff_eq;
/// use cobyla::{minimize, Func, RhoBeg};
///
/// fn paraboloid(x: &[f64], _data: &mut ()) -> f64 {
///     10. * (x[0] + 1.).powf(2.) + x[1].powf(2.)
/// }
///
/// let mut x = vec![1., 1.];
///
/// // Constraints definition to be positive eventually: here `x_0 > 0`
/// let cstr1 = |x: &[f64], _user_data: &mut ()| x[0];
/// let cons: Vec<&dyn Func<()>> = vec![&cstr1];
///
/// match minimize(
///     paraboloid,
///     &mut x,
///     &[(-10., 10.), (-10., 10.)],
///     &cons,
///     (),
///     200,
///     RhoBeg::All(0.5),
///     None
/// ) {
///     Ok((status, x_opt, y_opt)) => {
///         println!("status = {:?}", status);
///         println!("x_opt = {:?}", x_opt);
///         println!("y_opt = {}", y_opt);
/// #        assert_abs_diff_eq!(y_opt, 10.0);
///     }
///     Err((e, _, _)) => println!("Optim error: {:?}", e),
/// }
/// ```
///
/// ## Algorithm description:
///
/// COBYLA minimizes an objective function F(X) subject to M inequality
/// constraints on X, where X is a vector of variables that has N components.
///
/// The algorithm employs linear approximations to the objective and constraint
/// functions, the approximations being formed by linear interpolation at N+1
/// points in the space of the variables. We regard these interpolation points
/// as vertices of a simplex.
///
/// The parameter RHO controls the size of the simplex and it is reduced
/// automatically from RHOBEG to RHOEND. For each RHO the subroutine tries
/// to achieve a good vector of variables for the current size, and then RHO
/// is reduced until the value RHOEND is reached.
///
/// Therefore RHOBEG and RHOEND should be set to reasonable initial changes to and the
/// required accuracy in the variables respectively, but this accuracy should be
/// viewed as a subject for experimentation because it is not guaranteed.
///  
/// The subroutine has an advantage over many of its competitors, however, which is
/// that it treats each constraint individually when calculating a change to the
/// variables, instead of lumping the constraints together into a single penalty
/// function.  
///
/// The name of the algorithm is derived from the phrase Constrained
/// Optimization BY Linear Approximations.
///
/// The user can set the values of RHOBEG and RHOEND, and must provide an
/// initial vector of variables in X. Further, the value of IPRINT should be
/// set to 0, 1, 2 or 3, which controls the amount of printing during the
/// calculation. Specifically, there is no output if IPRINT=0 and there is
/// output only at the end of the calculation if IPRINT=1.  
/// Otherwise each new value of RHO and SIGMA is printed.  
///
/// Further, the vector of variables and some function information are
/// given either when RHO is reduced or when each
/// new value of F(X) is computed in the cases IPRINT=2 or IPRINT=3
/// respectively.  Here SIGMA is a penalty parameter, it being assumed that a
/// change to X is an improvement if it reduces the merit function:
///
/// F(X)+SIGMA*MAX(0.0,-C1(X),-C2(X),...,-CM(X)),
///
/// where C1, C2, ..., CM denote the constraint functions that should become
/// nonnegative eventually, at least to the precision of RHOEND.  In the
/// printed output the displayed term that is multiplied by SIGMA is
/// called MAXCV, which stands for 'MAXimum Constraint Violation'.
///
/// This implementation is a translation/adaptation of [NLopt](https://github.com/stevengj/nlopt) 2.7.1
/// See [NLopt COBYLA](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#cobyla-constrained-optimization-by-linear-approximations) documentation.
#[allow(clippy::useless_conversion)]
#[allow(clippy::too_many_arguments)]
pub fn minimize<F: Func<U>, G: Func<U>, U: Clone>(
    func: F,
    xinit: &[f64],
    bounds: &[(f64, f64)],
    cons: &[G],
    args: U,
    maxeval: usize,
    rhobeg: RhoBeg,
    stop_tol: Option<StopTols>,
) -> Result<SuccessOutcome, FailOutcome> {
    let fn_cfg = Box::new(NLoptFunctionCfg {
        objective_fn: func,
        user_data: args.clone(),
    });
    let fn_cfg_ptr = Box::into_raw(fn_cfg) as *mut c_void;
    let mut cstr_tol = 0.0; // no cstr tolerance

    let mut cstr_cfg = cons
        .iter()
        .map(|c| {
            let c_cfg = Box::new(NLoptConstraintCfg {
                constraint_fn: c as &dyn Func<U>,
                user_data: args.clone(),
            });
            let c_cfg_ptr = Box::into_raw(c_cfg) as *mut c_void;

            nlopt_constraint {
                m: 1,
                f: Some(nlopt_constraint_raw_callback::<F, U>),
                pre: None,
                mf: None,
                f_data: c_cfg_ptr,
                tol: &mut cstr_tol,
            }
        })
        .collect::<Vec<_>>();

    let mut x = vec![0.; xinit.len()];
    x.copy_from_slice(xinit);
    let n = x.len() as u32;
    let m = cons.len() as u32;

    if bounds.len() != x.len() {
        panic!(
            "{}",
            format!(
                "Minimize Error: Bounds and x should have same size! Got {} for bounds and {} for x.",
                bounds.len(),
                x.len()
            )
        )
    }
    let lbs: Vec<f64> = bounds.iter().map(|b| b.0).collect();
    let ubs: Vec<f64> = bounds.iter().map(|b| b.1).collect();
    let x_weights = vec![0.; n as usize];
    let mut dx = match rhobeg {
        RhoBeg::All(val) => vec![val; n as usize],
        RhoBeg::Set(val) => {
            if val.len() != n as usize {
                panic!(
                    "{}",
                    format!(
                        "Minimize Error: xtol_abs should have x dim size ({}), got {}",
                        n,
                        val.len()
                    )
                )
            } else {
                val
            }
        }
    };
    let mut minf = f64::INFINITY;
    let mut nevals_p = 0;
    let mut force_stop = 0;

    let stop_tol = stop_tol.unwrap_or_default();
    let xtol_abs = if stop_tol.xtol_abs.is_empty() {
        std::ptr::null()
    } else if stop_tol.xtol_abs.len() != n as usize {
        panic!(
            "{}",
            format!(
                "Minimize Error: xtol_abs should have x dim size ({}), got {}",
                n,
                stop_tol.xtol_abs.len()
            )
        );
    } else {
        stop_tol.xtol_abs.as_ptr()
    };
    let mut stop = nlopt_stopping {
        n,
        minf_max: -f64::INFINITY,
        ftol_rel: stop_tol.ftol_rel,
        ftol_abs: stop_tol.ftol_abs,
        xtol_rel: stop_tol.xtol_rel,
        xtol_abs,
        x_weights: x_weights.as_ptr(),
        nevals_p: &mut nevals_p,
        maxeval: maxeval as i32,
        maxtime: 0.0,
        start: 0.0,
        force_stop: &mut force_stop,
        stop_msg: "".to_string(),
    };

    let status = unsafe {
        cobyla_minimize::<U>(
            n.into(),
            Some(nlopt_function_raw_callback::<F, U>),
            fn_cfg_ptr,
            m.into(),
            cstr_cfg.as_mut_ptr(),
            0,
            std::ptr::null_mut(),
            lbs.as_ptr(),
            ubs.as_ptr(),
            x.as_mut_ptr(),
            &mut minf,
            &mut stop,
            dx.as_mut_ptr(),
        )
    };

    // Convert the raw pointer back into a Box with the B::from_raw function,
    // allowing the Box destructor to perform the cleanup.
    unsafe {
        let _ = Box::from_raw(fn_cfg_ptr as *mut NLoptFunctionCfg<F, U>);
    };

    match status {
        -1 => Err((FailStatus::Failure, x, minf)),
        -2 => Err((FailStatus::InvalidArgs, x, minf)),
        -3 => Err((FailStatus::OutOfMemory, x, minf)),
        -4 => Err((FailStatus::RoundoffLimited, x, minf)),
        -5 => Err((FailStatus::ForcedStop, x, minf)),
        1 => Ok((SuccessStatus::Success, x, minf)),
        2 => Ok((SuccessStatus::StopValReached, x, minf)),
        3 => Ok((SuccessStatus::FtolReached, x, minf)),
        4 => Ok((SuccessStatus::XtolReached, x, minf)),
        5 => Ok((SuccessStatus::MaxEvalReached, x, minf)),
        6 => Ok((SuccessStatus::MaxTimeReached, x, minf)),
        _ => Err((FailStatus::UnexpectedError, x, minf)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    use crate::nlopt_cobyla::cobyla_minimize;
    use crate::nlopt_cobyla::nlopt_stopping;

    /////////////////////////////////////////////////////////////////////////
    // Second problem (see cobyla.c case 6)

    fn raw_paraboloid(
        _n: libc::c_uint,
        x: *const libc::c_double,
        _gradient: *mut libc::c_double,
        _func_data: *mut libc::c_void,
    ) -> libc::c_double {
        unsafe {
            let r1 = *x.offset(0) + 1.0;
            let r2 = *x.offset(1);
            10.0 * (r1 * r1) + (r2 * r2) as libc::c_double
        }
    }

    #[test]
    fn test_cobyla_minimize() {
        let mut x = vec![1., 1.];
        let mut lb = vec![-10.0, -10.0];
        let mut ub = vec![10.0, 10.0];
        let mut dx = vec![0.5, 0.5];
        let mut minf = f64::INFINITY;
        let mut nevals_p = 0;
        let mut force_stop = 0;

        let mut stop = nlopt_stopping {
            n: 2,
            minf_max: -f64::INFINITY,
            ftol_rel: 0.0,
            ftol_abs: 0.0,
            xtol_rel: 0.0,
            xtol_abs: std::ptr::null(),
            x_weights: std::ptr::null(),
            nevals_p: &mut nevals_p,
            maxeval: 1000,
            maxtime: 0.0,
            start: 0.0,
            force_stop: &mut force_stop,
            stop_msg: "".to_string(),
        };

        let res = unsafe {
            cobyla_minimize::<()>(
                2,
                Some(raw_paraboloid),
                std::ptr::null_mut(),
                0,
                std::ptr::null_mut(),
                0,
                std::ptr::null_mut(),
                lb.as_mut_ptr(),
                ub.as_mut_ptr(),
                x.as_mut_ptr(),
                &mut minf,
                &mut stop,
                dx.as_mut_ptr(),
            )
        };

        println!("status = {:?}", res);
        println!("x = {:?}", x);

        assert_abs_diff_eq!(x.as_slice(), [-1., 0.].as_slice(), epsilon = 1e-4);
    }

    fn paraboloid(x: &[f64], _data: &mut ()) -> f64 {
        10. * (x[0] + 1.).powf(2.) + x[1].powf(2.)
    }

    #[test]
    fn test_paraboloid() {
        let xinit = vec![1., 1.];

        // let mut cons: Vec<&dyn Func<()>> = vec![];
        let mut cons: Vec<&dyn Func<()>> = vec![];
        let cstr1 = |x: &[f64], _user_data: &mut ()| x[0];
        cons.push(&cstr1 as &dyn Func<()>);

        // x_opt = [0, 0]
        match minimize(
            paraboloid,
            &xinit,
            &[(-10., 10.), (-10., 10.)],
            &cons,
            (),
            200,
            RhoBeg::All(0.5),
            None,
        ) {
            Ok((_, x, _)) => {
                let exp = [0., 0.];
                for (act, exp) in x.iter().zip(exp.iter()) {
                    assert_abs_diff_eq!(act, exp, epsilon = 1e-3);
                }
            }
            Err((status, _, _)) => {
                panic!("{}", format!("Error status : {:?}", status));
            }
        }
    }

    fn fletcher9115(x: &[f64], _user_data: &mut ()) -> f64 {
        -x[0] - x[1]
    }

    fn cstr1(x: &[f64], _user_data: &mut ()) -> f64 {
        x[1] - x[0] * x[0]
    }
    fn cstr2(x: &[f64], _user_data: &mut ()) -> f64 {
        1. - x[0] * x[0] - x[1] * x[1]
    }

    #[test]
    fn test_fletcher9115() {
        let xinit = vec![1., 1.];

        let cons = vec![&cstr1 as &dyn Func<()>, &cstr2 as &dyn Func<()>];

        let stop_tol = StopTols {
            ftol_rel: 1e-4,
            xtol_rel: 1e-4,
            ..StopTols::default()
        };

        match minimize(
            fletcher9115,
            &xinit,
            &[(-10., 10.), (-10., 10.)],
            &cons,
            (),
            200,
            RhoBeg::All(0.5),
            Some(stop_tol),
        ) {
            Ok((_, x, _)) => {
                let sqrt_0_5: f64 = 0.5_f64.sqrt();
                let exp = [sqrt_0_5, sqrt_0_5];
                for (act, exp) in x.iter().zip(exp.iter()) {
                    assert_abs_diff_eq!(act, exp, epsilon = 1e-3);
                }
            }
            Err((status, _, _)) => {
                println!("Error status : {:?}", status);
                panic!("Test fail");
            }
        }
    }

    fn xsinx(x: &[f64], _user_data: &mut ()) -> f64 {
        //(x - 3.5) * ((x - 3.5) / std::f64::consts::PI).mapv(|v| v.sin())
        (x[0] - 3.5) * f64::sin((x[0] - 3.5) / std::f64::consts::PI)
    }

    #[test]
    fn test_xsinx() {
        let xinit = vec![10.];

        // let mut cons: Vec<&dyn Func<()>> = vec![];
        let mut cons: Vec<&dyn Func<()>> = vec![];
        let cstr1 = |x: &[f64], _user_data: &mut ()| 17. - x[0];
        cons.push(&cstr1 as &dyn Func<()>);

        // x_opt = [0, 0]
        match minimize(
            xsinx,
            &xinit,
            &[(0., 25.)],
            &cons,
            (),
            200,
            RhoBeg::All(0.5),
            None,
        ) {
            Ok((_, x, _)) => {
                let exp = [18.935];
                let exp = [17.];
                for (act, exp) in x.iter().zip(exp.iter()) {
                    assert_abs_diff_eq!(act, exp, epsilon = 1e-2);
                }
            }
            Err((status, _, _)) => {
                panic!("{}", format!("Error status : {:?}", status));
            }
        }
    }
}
