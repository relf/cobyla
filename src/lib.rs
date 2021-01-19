#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
use std::os::raw::c_void;
use std::slice;

/// A trait to represent objective function
/// where `T` is the type of user data
pub trait ObjFn<T>: Fn(&[f64], &mut T) -> f64 {}
impl<F, T> ObjFn<T> for F where F: Fn(&[f64], &mut T) -> f64 {}

/// A trait to represent constraint function
/// which should be positive eventually  
pub trait CstrFn: Fn(&[f64]) -> f64 {}
impl<F> CstrFn for F where F: Fn(&[f64]) -> f64 {}

/// Packs a function with a user defined parameter set of type `T`
/// and constraints to be made positive eventually by the optimizer
struct FunctionCfg<'a, F: ObjFn<T>, G: CstrFn, T> {
    pub func: F,
    pub cons: &'a [G],
    pub data: T,
}

/// Callback interface for Cobyla C code to evaluate objective and constraint functions
extern "C" fn function_raw_callback<F: ObjFn<T>, G: CstrFn, T>(
    n: ::std::os::raw::c_long,
    m: ::std::os::raw::c_long,
    x: *const f64,
    con: *mut f64,
    data: *mut ::std::os::raw::c_void,
) -> f64 {
    // prepare args
    let argument = unsafe { slice::from_raw_parts(x, n as usize) };
    // recover FunctionCfg object from supplied params and call
    let f = unsafe { &mut *(data as *mut FunctionCfg<F, G, T>) };
    let res = (f.func)(argument, &mut f.data);

    for i in 0..m as isize {
        unsafe {
            *con.offset(i) = (f.cons[i as usize])(argument);
        }
    }

    // Important: we don't want f to get dropped at this point
    #[allow(clippy::forget_ref)]
    std::mem::forget(f);
    res
}

/// Minimizes a function using the Constrained Optimization By Linear Approximation (COBYLA) method.
///
/// This method wraps a C implementation of the algorithm using bindgen.
/// This interface is modeled after [scypi.optimize.fmin_cobyla](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_cobyla.html#scipy-optimize-fmin-cobyla)
///
/// Example
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
/// let mut cons: Vec<&CstrFn> = vec![];
/// cons.push(&|x: &[f64]| x[1] - x[0] * x[0]);
/// cons.push(&|x: &[f64]| 1. - x[0] * x[0] - x[1] * x[1]);
///
/// let (status, x_opt) = fmin_cobyla(paraboloid, &mut x, &cons, (), 0.5, 1e-4, 200, 0);
/// println!("status = {}", status);
/// println!("x = {:?}", x_opt);
/// ```

#[allow(clippy::too_many_arguments)]
pub fn fmin_cobyla<'a, F: ObjFn<T>, G: CstrFn, T>(
    func: F,
    x0: &'a mut [f64],
    cons: &[G],
    args: T,
    rhobeg: f64,
    rhoend: f64,
    maxfun: i32,
    iprint: i32,
) -> (i32, &'a [f64]) {
    let n = x0.len();
    let m = cons.len();

    // Our strategy is to pass the actual objective function as part of the
    // parameters to the callback. For this we pack it inside a FunctionCfg struct.
    // We allocation our FunctionCfg on the heap and pass a pointer to the C lib
    // (This is pretty unsafe but works).
    // `into_raw` will leak the boxed object
    let fn_cfg = Box::new(FunctionCfg {
        func,
        cons,
        data: args, // move user_data into FunctionCfg
    });
    let fn_cfg_ptr = Box::into_raw(fn_cfg) as *mut c_void;

    let x = x0;
    let mut maxfun = maxfun;
    let mut w: Vec<f64> = vec![0.; n * (3 * n + 2 * m + 11) + 4 * m + 6];
    let mut iact: Vec<i32> = vec![0; m + 1];

    let status = unsafe {
        raw_cobyla(
            n as i32,
            m as i32,
            Some(function_raw_callback::<F, G, T>),
            fn_cfg_ptr,
            x.as_mut_ptr(),
            rhobeg,
            rhoend,
            iprint,
            &mut maxfun,
            w.as_mut_ptr(),
            iact.as_mut_ptr(),
        )
    };
    (status, x)
}

#[cfg(test)]
mod tests {
    use super::*;

    /////////////////////////////////////////////////////////////////////////
    // First problem

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

        let (status, x_opt) = fmin_cobyla(paraboloid, &mut x, &cons, (), 0.5, 1e-4, 200, 1);
        println!("status = {}", status);
        println!("x = {:?}", x_opt);
    }

    /// Direct usage Pb 1
    ///
    unsafe extern "C" fn calcfc(
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
        let iprint = 1;
        let mut maxfun = 2000;
        let mut w: Vec<f64> = vec![0.; 3000];
        let mut iact: Vec<i32> = vec![0; 51];
        let &null = &0;

        // xopt = [-1., 0.]
        unsafe {
            raw_cobyla(
                n,
                m,
                Some(calcfc),
                null as *mut _,
                x.as_mut_ptr(),
                rhobeg,
                rhoend,
                iprint,
                &mut maxfun,
                w.as_mut_ptr(),
                iact.as_mut_ptr(),
            );
        }
    }

    /////////////////////////////////////////////////////////////////////////
    // Second problem

    #[test]
    fn test_fmin_cobyla2() {
        let mut x = vec![1., 1.];

        #[allow(bare_trait_objects)]
        let mut cons: Vec<&CstrFn> = vec![];
        cons.push(&|x: &[f64]| x[1] - x[0] * x[0]);
        cons.push(&|x: &[f64]| 1. - x[0] * x[0] - x[1] * x[1]);

        let (status, x_opt) = fmin_cobyla(paraboloid, &mut x, &cons, (), 0.5, 1e-4, 200, 0);
        println!("status = {}", status);
        println!("x = {:?}", x_opt);
    }

    /// Direct usage Pb 2
    ///
    unsafe extern "C" fn calcfc_cstr(
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
        let iprint = 1;
        let mut maxfun = 2000;
        let mut w: Vec<f64> = vec![0.; 3000];
        let mut iact: Vec<i32> = vec![0; 51];
        let &null = &0;

        // xopt = [-1., 0.]
        unsafe {
            raw_cobyla(
                n,
                m,
                Some(calcfc_cstr),
                null as *mut _,
                x.as_mut_ptr(),
                rhobeg,
                rhoend,
                iprint,
                &mut maxfun,
                w.as_mut_ptr(),
                iact.as_mut_ptr(),
            );
        }
    }
}
