//! cobyla
//!
//! COBYLA is an algorithm for minimizing a function of many variables. The method is derivatives free (only the function values are needed)
//! and take into account constraints on the variables. The algorithm is described in:
//!
//! > M.J.D. Powell, "A direct search optimization method that models the objective and constraint functions by linear interpolation," in
//! > Advances in Optimization and Numerical Analysis Mathematics and Its Applications, vol. 275 (eds. Susana Gomez and Jean-Pierre Hennart),
//! > Kluwer Academic Publishers, pp. 51-67 (1994).
//!
//! The objective function to be minimized has to implement the [`ObjFn`] trait, while constraints,
//! also defined as functions of the input variables have to implement the [`CstrFn`] trait.
//! Constraints values are intended to be positive at the end of the optimization.
//!
//! The algorithm can be run either using the [`fmin_cobyla`] function or using the [`CobylaSolver`]
//! which leverages the [argmin](https://www.argmin-rs.org/book/index.html) framework.
//!
//! Implementation Notes:
//!
//! 0.3.x : COBYLA is now also implemented as an argmin::Solver and benefits from [argmin framework](https://github.com/argmin-rs) tooling.
//!
//! 0.2.x : The C code is now translated in Rust using c2rust transpiler then manually edited to avoid FFI usage
//! to get Rust (unsafe) implementation.
//!  
//! 0.1.x : the C code is wrapped with with bindgen is visible as the `raw_cobyla` function using the callback type
//! `cobyla_calcfc` which is used to compute the objective function and the constraints.
//!
mod cobyla;
use nlopt_cobyla::nlopt_constraint;

use crate::cobyla::raw_cobyla;

mod nlopt_cobyla;
use crate::nlopt_cobyla::{
    cobyla_minimize, nlopt_constraint_raw_callback, nlopt_eval_constraint,
    nlopt_function_raw_callback, nlopt_stopping, NLoptConstraintCfg, NLoptFunctionCfg, NLoptObjFn,
};

mod cobyla_solver;
mod cobyla_state;
pub use crate::cobyla_solver::*;
pub use crate::cobyla_state::*;

use std::os::raw::c_void;
use std::slice;
/// A trait for an objective function to be minimized
///
/// An objective function takes the form of a closure `f(x: &[f64], user_data: &mut U) -> f64`
///
/// * `x` - n-dimensional array
/// * `user_data` - user defined data
pub trait ObjFn<U>: Fn(&[f64], &mut U) -> f64 {}
impl<F, U> ObjFn<U> for F where F: Fn(&[f64], &mut U) -> f64 {}

/// A trait for a constraint function which should be positive eventually
///
/// A constraint function takes the form of a closure `f(x: &[f64]) -> f64`
/// The algorithm makes the constraint positive eventually.
///
/// For instance if you want an upper bound MAX for x,
/// you have to define the constraint as `|x| MAX - x`.
/// Conversly for a lower bound you would define `|x| x - MIN`
///
/// * `x` - n-dimensional array
pub trait CstrFn: Fn(&[f64]) -> f64 {}
impl<F> CstrFn for F where F: Fn(&[f64]) -> f64 {}

/// Packs a function with a user defined parameter set of type `U`
/// and constraints to be made positive eventually by the optimizer
struct FunctionCfg<'a, F: ObjFn<U>, G: CstrFn, U> {
    pub func: F,
    pub cons: &'a [G],
    pub data: U,
}

/// Callback interface for Cobyla C code to evaluate objective and constraint functions
fn function_raw_callback<F: ObjFn<U>, G: CstrFn, U>(
    n: ::std::os::raw::c_long,
    m: ::std::os::raw::c_long,
    x: *const f64,
    con: *mut f64,
    data: *mut ::std::os::raw::c_void,
) -> f64 {
    // prepare args
    let argument = unsafe { slice::from_raw_parts(x, n as usize) };
    // recover FunctionCfg object from supplied params and call
    let f = unsafe { &mut *(data as *mut FunctionCfg<F, G, U>) };
    let res = (f.func)(argument, &mut f.data);

    for i in 0..m as isize {
        unsafe {
            *con.offset(i) = (f.cons[i as usize])(argument);
        }
    }

    // Important: we don't want f to get dropped at this point
    #[allow(forgetting_references)]
    std::mem::forget(f);
    res
}

/// Minimizes a function using the Constrained Optimization By Linear Approximation (COBYLA) method.
///
/// This interface is modeled after [scypi.optimize.fmin_cobyla](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_cobyla.html#scipy-optimize-fmin-cobyla)
///
/// # Example
/// ```
/// use cobyla::{fmin_cobyla, CstrFn};
///
/// fn paraboloid(x: &[f64], _data: &mut ()) -> f64 {
///     10. * (x[0] + 1.).powf(2.) + x[1].powf(2.)
/// }
///
/// let mut x = vec![1., 1.];
///
/// // Constraints definition to be positive eventually
/// let mut cons: Vec<&dyn CstrFn> = vec![];
/// cons.push(&|x: &[f64]| x[1] - x[0] * x[0]);
/// cons.push(&|x: &[f64]| 1. - x[0] * x[0] - x[1] * x[1]);
///
/// let (status, x_opt) = fmin_cobyla(paraboloid, &mut x, &cons, (), 0.5, 1e-4, 200, 0);
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
#[allow(clippy::useless_conversion)]
#[allow(clippy::too_many_arguments)]
pub fn fmin_cobyla<'a, F: ObjFn<U>, G: CstrFn, U>(
    func: F,
    x0: &'a mut [f64],
    cons: &[G],
    args: U,
    rhobeg: f64,
    rhoend: f64,
    maxfun: i32,
    iprint: i32,
) -> (i32, &'a [f64]) {
    let n: i32 = x0.len() as i32;
    let m: i32 = cons.len() as i32;

    // Trick from nlopt-rust: https://github.com/adwhit/rust-nlopt
    // The strategy is to pass the actual objective function as part of the
    // parameters to the callback. For this we pack it inside a FunctionCfg struct.
    // We allocate our FunctionCfg on the heap and pass a pointer to the C lib
    // (This is pretty unsafe but it works).
    let fn_cfg = Box::new(FunctionCfg {
        func,
        cons,
        data: args, // move user_data into FunctionCfg
    });
    let fn_cfg_ptr = Box::into_raw(fn_cfg) as *mut c_void;

    let x = x0;
    let mut maxfun = maxfun.into();
    let mut w = vec![0.; (n * (3 * n + 2 * m + 11) + 4 * m + 6) as usize];
    let mut iact = vec![0; (m + 1) as usize];
    let status = unsafe {
        raw_cobyla(
            n.into(),
            m.into(),
            Some(function_raw_callback::<F, G, U>),
            fn_cfg_ptr,
            x.as_mut_ptr(),
            rhobeg,
            rhoend,
            iprint.into(),
            &mut maxfun,
            w.as_mut_ptr(),
            iact.as_mut_ptr(),
        )
    };
    // Convert the raw pointer back into a Box with the Box::from_raw function,
    // allowing the Box destructor to perform the cleanup.
    unsafe {
        let _ = Box::from_raw(fn_cfg_ptr as *mut FunctionCfg<'_, F, G, U>);
    };
    (status, x)
}

/// Cobyla implementation generated from NLopt 2.7.1: https://github.com/stevengj/nlopt
/// and plugged as what is done in the NLopt rust binding: https://github.com/adwhit/rust-nlopt
/// but using the same API as fmin_cobyla (which is a bit backward as NLopt API is richer)
/// The idea here is to be able to easily to switch between the two implementations.
#[allow(clippy::useless_conversion)]
#[allow(clippy::too_many_arguments)]
pub fn nlopt_cobyla<'a, F: NLoptObjFn<U>, G: NLoptObjFn<U>, U: Clone>(
    func: F,
    x0: &'a mut [f64],
    cons: &[G],
    args: U,
    rhobeg: f64,
    rhoend: f64,
    maxfun: i32,
    _iprint: i32,
) -> (i32, &'a [f64]) {
    let fn_cfg = Box::new(NLoptFunctionCfg {
        objective_fn: func,
        user_data: args.clone(), // move user_data into FunctionCfg
    });
    let fn_cfg_ptr = Box::into_raw(fn_cfg) as *mut c_void;
    let mut cstr_tol = 1e-4;

    let mut cstr_cfg = cons
        .iter()
        .map(|c| {
            let c_cfg = Box::new(NLoptConstraintCfg {
                constraint_fn: c as &dyn NLoptObjFn<U>,
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

    let lb = vec![-f64::INFINITY; n as usize];
    let ub = vec![f64::INFINITY; n as usize];
    let xtol_abs = vec![0.; n as usize];
    let x_weights = vec![0.; n as usize];
    let mut dx = vec![rhobeg; n as usize];
    let mut minf = f64::INFINITY;
    let mut nevals_p = 0;
    let mut force_stop = 0;
    let mut stop = nlopt_stopping {
        n,
        minf_max: -f64::INFINITY,
        ftol_rel: 0.0,
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
    if !cons.is_empty() {
        let mut xtest = vec![1., 1.];
        unsafe {
            let mut result = -666.;
            let nc = cstr_cfg[0];
            let fc = &nc as *const nlopt_constraint;

            // It works: cstr1 is called
            let _res = nlopt_constraint_raw_callback::<&dyn NLoptObjFn<()>, ()>(
                2,
                xtest.as_mut_ptr(),
                std::ptr::null_mut::<libc::c_double>(),
                (nc).f_data,
            );
            // println!(
            //     "###################################### JUST direct nlopt_constraint_raw_callback is OK = {}",
            //     res,
            // );

            // XXX: Weird bug!
            // If fails : (*fc).f does call nlopt_constraint_raw_callback but
            // when unpacking (*fc).f_data with unsafe f = { &mut *(params as *mut NLoptConstraintCfg<F, T>) };
            // (*f).constraint_fn(...) calls the objective function instead of the constraint function!!!
            let _res = ((*fc).f.expect("func"))(
                2,
                xtest.as_mut_ptr(),
                std::ptr::null_mut::<libc::c_double>(),
                (*fc).f_data,
            );
            // println!(
            //     "###################################### JUST stored nlopt_constraint_raw_callback is NOT OK = {}",
            //     res,
            // );

            // It works: cstr1 is called
            // we use directly a copy of specialized nlopt_constraint_raw_callback
            nlopt_eval_constraint::<()>(
                &mut result,
                std::ptr::null_mut::<libc::c_double>(),
                fc,
                n,
                x.as_mut_ptr(),
            );
            // println!(
            //     "############################### TEST nlopt_eval_constraint (OK if opposite previous OK result) = {}",
            //     result
            // );
        }
    }

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
    // First problem (see cobyla.c case 1)

    fn paraboloid(x: &[f64], _data: &mut ()) -> f64 {
        10. * (x[0] + 1.).powf(2.) + x[1].powf(2.)
    }

    #[test]
    fn test_fmin_cobyla() {
        let mut x = vec![1., 1.];

        #[allow(bare_trait_objects)]
        let mut cons: Vec<&CstrFn> = vec![];
        let cstr1 = |x: &[f64]| x[0];
        cons.push(&cstr1);

        // x_opt = [0, 0]
        let (status, x_opt) = fmin_cobyla(paraboloid, &mut x, &cons, (), 0.5, 1e-4, 200, 0);
        println!("status = {}", status);
        println!("x = {:?}", x_opt);

        assert_abs_diff_eq!(x.as_slice(), [0., 0.].as_slice(), epsilon = 1e-4);
    }

    /// Direct usage Pb 1
    ///
    unsafe fn calcfc(
        _n: ::std::os::raw::c_long,
        _m: ::std::os::raw::c_long,
        x: *const f64,
        _con: *mut f64,
        _data: *mut ::std::os::raw::c_void,
    ) -> f64 {
        let r1 = *x.offset(0) + 1.0;
        let r2 = *x.offset(1);
        10.0 * (r1 * r1) + (r2 * r2)
    }

    #[test]
    fn test_cobyla() {
        let n = 2;
        let m = 0;
        let mut x = vec![1.0, 1.0];
        let rhobeg = 0.5;
        let rhoend = 1e-4;
        let iprint = 0;
        let mut maxfun = 2000;
        let mut w: Vec<_> = vec![0.; 3000];
        let mut iact: Vec<_> = vec![0; 51];
        let null = std::ptr::null_mut::<c_void>();

        // xopt = [-1., 0.]
        unsafe {
            raw_cobyla(
                n,
                m,
                Some(calcfc),
                null,
                x.as_mut_ptr(),
                rhobeg,
                rhoend,
                iprint,
                &mut maxfun,
                w.as_mut_ptr(),
                iact.as_mut_ptr(),
            );
        }
        assert_abs_diff_eq!(x.as_slice(), [-1., 0.].as_slice(), epsilon = 1e-3);
    }

    /////////////////////////////////////////////////////////////////////////
    // Second problem (see cobyla.c case 6)

    fn fletcher9115(x: &[f64], _user_data: &mut ()) -> f64 {
        -x[0] - x[1]
    }

    #[allow(clippy::vec_init_then_push)]
    #[test]
    fn test_fmin_cobyla2() {
        let mut x = vec![1., 1.];

        let mut cons: Vec<&dyn CstrFn> = vec![];
        cons.push(&|x: &[f64]| x[1] - x[0] * x[0]);
        cons.push(&|x: &[f64]| 1. - x[0] * x[0] - x[1] * x[1]);

        let (status, x_opt) = fmin_cobyla(fletcher9115, &mut x, &cons, (), 0.5, 1e-4, 200, 0);
        println!("status = {}", status);
        println!("x = {:?}", x_opt);

        let sqrt_0_5: f64 = 0.5_f64.sqrt();
        assert_abs_diff_eq!(x_opt, [sqrt_0_5, sqrt_0_5].as_slice(), epsilon = 1e-4);
    }

    /// Direct usage Pb 6
    ///
    unsafe fn calcfc_cstr(
        _n: ::std::os::raw::c_long,
        _m: ::std::os::raw::c_long,
        x: *const f64,
        con: *mut f64,
        _data: *mut ::std::os::raw::c_void,
    ) -> f64 {
        let r1 = *x.offset(0);
        let r2 = *x.offset(1);
        let fc = -r1 - r2;
        *con.offset(0) = r2 - r1 * r1;
        *con.offset(1) = 1.0 - r1 * r1 - r2 * r2;
        fc
    }

    #[test]
    fn test_cobyla_cstr() {
        let n = 2;
        let m = 2;
        let mut x = vec![1.0, 1.0];
        let rhobeg = 0.5;
        let rhoend = 1e-4;
        let iprint = 0;
        let mut maxfun = 2000;
        let mut w: Vec<_> = vec![0.; 3000];
        let mut iact: Vec<_> = vec![0; 51];
        let null = std::ptr::null_mut::<c_void>();

        // xopt = [sqrt(0.5), sqrt(0.5)]
        unsafe {
            raw_cobyla(
                n,
                m,
                Some(calcfc_cstr),
                null,
                x.as_mut_ptr(),
                rhobeg,
                rhoend,
                iprint,
                &mut maxfun,
                w.as_mut_ptr(),
                iact.as_mut_ptr(),
            );
        }

        let sqrt_0_5: f64 = 0.5_f64.sqrt();
        assert_abs_diff_eq!(
            x.as_slice(),
            [sqrt_0_5, sqrt_0_5].as_slice(),
            epsilon = 1e-4
        );
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// NlOpt cobyla

    fn nlopt_raw_paraboloid(
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
                Some(nlopt_raw_paraboloid),
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

    fn nlopt_paraboloid(x: &[f64], _g: Option<&mut [f64]>, _data: &mut ()) -> f64 {
        10. * (x[0] + 1.).powf(2.) + x[1].powf(2.)
    }

    #[test]
    fn test_nlopt_paraboloid() {
        let mut x = vec![1., 1.];

        // #[allow(bare_trait_objects)]
        let mut cons: Vec<&dyn NLoptObjFn<()>> = vec![];
        let cstr1 = |x: &[f64], _g: Option<&mut [f64]>, _user_data: &mut ()| x[0];
        cons.push(&cstr1 as &dyn NLoptObjFn<()>);

        // x_opt = [0, 0]
        let (status, x_opt) = nlopt_cobyla(nlopt_paraboloid, &mut x, &cons, (), 0.5, 1e-6, 200, 1);
        println!("status = {}", status);
        println!("x = {:?}", x_opt);

        assert_abs_diff_eq!(x.as_slice(), [0., 0.].as_slice(), epsilon = 1e-4);
    }

    fn nlopt_fletcher9115(x: &[f64], _g: Option<&mut [f64]>, _user_data: &mut ()) -> f64 {
        -x[0] - x[1]
    }

    fn cstr1(x: &[f64], _g: Option<&mut [f64]>, _user_data: &mut ()) -> f64 {
        x[1] - x[0] * x[0]
    }
    fn cstr2(x: &[f64], _g: Option<&mut [f64]>, _user_data: &mut ()) -> f64 {
        1. - x[0] * x[0] - x[1] * x[1]
    }

    #[test]
    fn test_nlopt_fletcher() {
        let mut x = vec![1., 1.];

        let cons = vec![&cstr1 as &dyn NLoptObjFn<()>, &cstr2 as &dyn NLoptObjFn<()>];

        let (status, x_opt) =
            nlopt_cobyla(nlopt_fletcher9115, &mut x, &cons, (), 0.5, 1e-4, 200, 1);
        println!("status = {}", status);
        println!("x = {:?}", x_opt);

        let sqrt_0_5: f64 = 0.5_f64.sqrt();
        assert_abs_diff_eq!(x_opt, [sqrt_0_5, sqrt_0_5].as_slice(), epsilon = 1e-4);
    }
}
