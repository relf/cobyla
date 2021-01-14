#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[cfg(test)]
mod tests {
    use super::*;
    use std::slice;
    use std::os::raw::{c_void};


    pub trait ObjFn<T>: Fn(&[f64], &mut T) -> f64 {}
    impl<F, T> ObjFn<T> for F where F: Fn(&[f64], &mut T) -> f64 {}
    
    /// Packs an objective function with a user defined parameter set of type `T`.
    struct FunctionCfg<F: ObjFn<T>, T> {
        pub objective_fn: F,
        pub user_data: T,
    }

    pub trait CstrFn<T>: Fn(&[f64], &mut T) -> f64 {}
    impl<F, T> CstrFn<T> for F where F: Fn(&[f64], &mut T) -> f64 {}

    extern "C" fn function_raw_callback<F: ObjFn<T>, T>(
        n: ::std::os::raw::c_long,
        _m: ::std::os::raw::c_long,
        x: *const f64,
        _con: *mut f64,
        data: *mut ::std::os::raw::c_void,
    ) -> f64 {
        // prepare args
        let argument = unsafe { slice::from_raw_parts(x, n as usize) };
        // recover FunctionCfg object from supplied params and call
        let f = unsafe { &mut *(data as *mut FunctionCfg<F, T>) };
        let res = (f.objective_fn)(argument, &mut f.user_data);
        // Important: we don't want f to get dropped at this point
        std::mem::forget(f);
        res
    }


    #[allow(clippy::too_many_arguments)]
    fn fmin_cobyla<'a, F: ObjFn<T1>, G: CstrFn<T2>, T1, T2>(
        func: F,
        x0: &'a mut [f64],
        cons: &[G],
        args: T1,
        _consargs: T2,
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
            objective_fn: func,
            user_data: args, // move user_data into FunctionCfg
        });
        let fn_cfg_ptr = Box::into_raw(fn_cfg) as *mut c_void;

        // TODO: manage constraints args
        // if let Some(cargs) = consargs {...}

        let x = x0;
        let mut maxfun = maxfun;
        let mut w: Vec<f64> = vec![0.; n * (3 * n + 2 * m + 11) + 4 * m + 6];
        let mut iact: Vec<i32> = vec![0; m + 1];

        let status = unsafe {
            cobyla(
                n as i32,
                m as i32,
                Some(function_raw_callback::<F, T1>),
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

    fn paraboloid(x: &[f64], _data: &mut ()) -> f64 {
        10.*(x[0]+1.).powf(2.) + x[1].powf(2.)
    }

    #[test]
    fn test_fmin_cobyla() {
        let mut x = vec![1., 1.];

        #[allow(bare_trait_objects)]
        let cons: Vec<&CstrFn<()>> = vec![];

        let (status, x_opt) = fmin_cobyla(paraboloid, &mut x, &cons, (), (), 0.5, 1e-4, 200, 0);
        println!("status = {}", status);
        println!("x = {:?}", x_opt);
    }

    /// Direct usage
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
        let fc = 10.0 * (r1 * r1) + (r2 * r2);
        fc
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
            cobyla(
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
}
