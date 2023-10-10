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

/// Minimizes a function using the Constrained Optimization By Linear Approximation (COBYLA) method.
///
/// # Example
/// ```
/// use cobyla::{minimize, Func};
///
/// fn paraboloid(x: &[f64], _data: &mut ()) -> f64 {
///     10. * (x[0] + 1.).powf(2.) + x[1].powf(2.)
/// }
///
/// let mut x = vec![1., 1.];
///
/// // Constraints definition to be positive eventually
/// let mut cons: Vec<&dyn Func<()>> = vec![];
/// cons.push(&|x: &[f64], _data: &mut ()| x[1] - x[0] * x[0]);
/// cons.push(&|x: &[f64], _data: &mut ()| 1. - x[0] * x[0] - x[1] * x[1]);
///
/// let (status, x_opt) = minimize(paraboloid, &mut x, &cons, (), 0.5, 1e-4, 200, 0, (-10., 10.));
/// println!("status = {}", status);
/// println!("x = {:?}", x_opt);
/// ```
///
/// # Algorithm description:
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
/// This implementation is a translation of [NLopt](https://github.com/stevengj/nlopt) 2.7.1
/// See [NLopt COBYLA](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#cobyla-constrained-optimization-by-linear-approximations) documentation.
#[allow(clippy::useless_conversion)]
#[allow(clippy::too_many_arguments)]
pub fn minimize<'a, F: Func<U>, G: Func<U>, U: Clone>(
    func: F,
    x0: &'a mut [f64],
    cons: &[G],
    args: U,
    rhobeg: f64,
    rhoend: f64,
    maxfun: i32,
    _iprint: i32,
    bounds: (f64, f64),
) -> (i32, &'a [f64]) {
    let fn_cfg = Box::new(NLoptFunctionCfg {
        objective_fn: func,
        user_data: args.clone(), // move user_data into FunctionCfg
    });
    let fn_cfg_ptr = Box::into_raw(fn_cfg) as *mut c_void;
    let mut cstr_tol = 2e-4;

    let mut cstr_cfg = cons
        .iter()
        .map(|c| {
            let c_cfg = Box::new(NLoptConstraintCfg {
                constraint_fn: c as &dyn Func<U>,
                user_data: args.clone(), // move user_data into FunctionCfg
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

    let x = x0;
    let n = x.len() as u32;
    let m = cons.len() as u32;

    let lb = vec![bounds.0; n as usize];
    let ub = vec![bounds.1; n as usize];
    let xtol_abs = vec![0.; n as usize];
    let x_weights = vec![0.; n as usize];
    let mut dx = vec![rhobeg; n as usize];
    let mut minf = f64::INFINITY;
    let mut nevals_p = 0;
    let mut force_stop = 0;
    let mut stop = nlopt_stopping {
        n,
        minf_max: -f64::INFINITY,
        ftol_rel: 1e-4,
        ftol_abs: 0.0,
        xtol_rel: rhoend / rhobeg,
        xtol_abs: xtol_abs.as_ptr(),
        x_weights: x_weights.as_ptr(),
        nevals_p: &mut nevals_p,
        maxeval: maxfun.into(),
        maxtime: 0.0,
        start: 0.0,
        force_stop: &mut force_stop,
        stop_msg: "".to_string(),
    };

    // XXX: Weird bug. Can not pass nlopt_constraint_raw_callback
    // Work around is to patch nlopt_eval_constraint to use nlopt_constraint_raw_callback directly
    // if !cons.is_empty() {
    //     let mut xtest = vec![1., 1.];
    //     unsafe {
    //         let mut result = -666.;
    //         let nc = cstr_cfg[0];
    //         let fc = &nc as *const nlopt_constraint;

    //         // It works: cstr1 is called
    //         let _res = nlopt_constraint_raw_callback::<&dyn Func<()>, ()>(
    //             2,
    //             xtest.as_mut_ptr(),
    //             std::ptr::null_mut::<libc::c_double>(),
    //             (nc).f_data,
    //         );
    //         // println!(
    //         //     "###################################### JUST direct nlopt_constraint_raw_callback is OK = {}",
    //         //     res,
    //         // );

    //         // XXX: Weird bug!
    //         // If fails : (*fc).f does call nlopt_constraint_raw_callback but
    //         // when unpacking (*fc).f_data with unsafe f = { &mut *(params as *mut NLoptConstraintCfg<F, T>) };
    //         // (*f).constraint_fn(...) calls the objective function instead of the constraint function!!!
    //         let _res = ((*fc).f.expect("func"))(
    //             2,
    //             xtest.as_mut_ptr(),
    //             std::ptr::null_mut::<libc::c_double>(),
    //             (*fc).f_data,
    //         );
    //         // println!(
    //         //     "###################################### JUST stored nlopt_constraint_raw_callback is NOT OK = {}",
    //         //     res,
    //         // );

    //         // It works: cstr1 is called
    //         // we use directly a copy of specialized nlopt_constraint_raw_callback
    //         nlopt_eval_constraint::<()>(
    //             &mut result,
    //             std::ptr::null_mut::<libc::c_double>(),
    //             fc,
    //             n,
    //             x.as_mut_ptr(),
    //         );
    //         // println!(
    //         //     "############################### TEST nlopt_eval_constraint (OK if opposite previous OK result) = {}",
    //         //     result
    //         // );
    //     }
    // }

    let status = unsafe {
        cobyla_minimize::<U>(
            n.into(),
            Some(nlopt_function_raw_callback::<F, U>),
            fn_cfg_ptr,
            m.into(),
            cstr_cfg.as_mut_ptr(),
            0,
            std::ptr::null_mut(),
            lb.as_ptr(),
            ub.as_ptr(),
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
    (status, x)
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
    fn test_nlopt_cobyla_minimize() {
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
        let mut x = vec![1., 1.];

        // #[allow(bare_trait_objects)]
        let mut cons: Vec<&dyn Func<()>> = vec![];
        let cstr1 = |x: &[f64], _user_data: &mut ()| x[0];
        cons.push(&cstr1 as &dyn Func<()>);

        // x_opt = [0, 0]
        let (status, x_opt) =
            minimize(paraboloid, &mut x, &cons, (), 0.5, 0.0, 200, 1, (-10., 10.));
        println!("status = {}", status);
        println!("x = {:?}", x_opt);

        assert_abs_diff_eq!(x.as_slice(), [0., 0.].as_slice(), epsilon = 1e-3);
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
        let mut x = vec![1., 1.];

        let cons = vec![&cstr1 as &dyn Func<()>, &cstr2 as &dyn Func<()>];

        let (status, x_opt) = minimize(
            fletcher9115,
            &mut x,
            &cons,
            (),
            0.5,
            1e-4,
            200,
            1,
            (-10., 10.),
        );
        println!("status = {}", status);
        println!("x = {:?}", x_opt);

        let sqrt_0_5: f64 = 0.5_f64.sqrt();
        assert_abs_diff_eq!(x_opt, [sqrt_0_5, sqrt_0_5].as_slice(), epsilon = 1e-4);
    }
}
