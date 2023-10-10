#![allow(
    dead_code,
    mutable_transmutes,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_assignments,
    unused_mut,
    clippy::needless_return,
    clippy::zero_ptr,
    clippy::toplevel_ref_arg,
    clippy::nonminimal_bool,
    clippy::assign_op_pattern,
    clippy::collapsible_if,
    clippy::neg_cmp_op_on_partial_ord,
    clippy::single_match,
    clippy::unnecessary_cast,
    clippy::excessive_precision,
    clippy::too_many_arguments
)]

use std::convert::TryFrom;
use std::time::{SystemTime, UNIX_EPOCH};

use std::slice;

pub fn nlopt_function_raw_callback<F: Func<T>, T>(
    n: libc::c_uint,
    x: *const f64,
    g: *mut f64,
    params: *mut libc::c_void,
) -> f64 {
    // prepare args
    let argument = unsafe { slice::from_raw_parts(x, n as usize) };
    let gradient = if g.is_null() {
        None
    } else {
        Some(unsafe { slice::from_raw_parts_mut(g, n as usize) })
    };

    // recover FunctionCfg object from supplied params and call
    let f = unsafe { &mut *(params as *mut NLoptFunctionCfg<F, T>) };
    let res = (f.objective_fn)(argument, gradient, &mut f.user_data);
    #[allow(forgetting_references)]
    std::mem::forget(f);
    res
}

pub fn nlopt_constraint_raw_callback<F: Func<T>, T>(
    n: libc::c_uint,
    x: *const f64,
    g: *mut f64,
    params: *mut libc::c_void,
) -> f64 {
    let f = unsafe { &mut *(params as *mut NLoptConstraintCfg<F, T>) };
    let argument = unsafe { slice::from_raw_parts(x, n as usize) };
    let gradient = if g.is_null() {
        None
    } else {
        Some(unsafe { slice::from_raw_parts_mut(g, n as usize) })
    };
    (f.constraint_fn)(argument, gradient, &mut f.user_data)
}

/// Packs an objective function with a user defined parameter set of type `T`.
pub struct NLoptFunctionCfg<F: Func<T>, T> {
    pub objective_fn: F,
    pub user_data: T,
}

pub struct NLoptConstraintCfg<F: Func<T>, T> {
    pub constraint_fn: F,
    pub user_data: T,
}

/// A trait representing an objective function.
///
/// An objective function takes the form of a closure `f(x: &[f64], gradient: Option<&mut [f64], user_data: &mut U) -> f64`
///
/// * `x` - `n`-dimensional array
/// * `gradient` - `n`-dimensional array to store the gradient `grad f(x)`. If `gradient` matches
/// `Some(x)`, the user is required to provide a gradient, otherwise the optimization will
/// probabely fail.
/// * `user_data` - user defined data
pub trait Func<U>: Fn(&[f64], Option<&mut [f64]>, &mut U) -> f64 {}
impl<T, U> Func<U> for T where T: Fn(&[f64], Option<&mut [f64]>, &mut U) -> f64 {}

enum Io {
    stderr,
    stdout,
}
fn fprintf(_io: Io, msg: &str) {
    println!("{}", msg);
}

// #![register_tool(c2rust)]
// #![feature(c_variadic, extern_types, register_tool)]
//  {
//     pub type _IO_wide_data;
//     pub type _IO_codecvt;
//     pub type _IO_marker;
//     fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
//     fn realloc(_: *mut libc::c_void, _: libc::c_ulong) -> *mut libc::c_void;
//     fn free(__ptr: *mut libc::c_void);
//     fn abort() -> !;
//     static mut stderr: *mut FILE;
//     fn fprintf(_: *mut FILE, _: *const libc::c_char, _: ...) -> libc::c_int;
//     fn vsnprintf(
//         _: *mut libc::c_char,
//         _: libc::c_ulong,
//         _: *const libc::c_char,
//         _: ::std::ffi::VaList,
//     ) -> libc::c_int;
//     fn sqrt(_: libc::c_double) -> libc::c_double;
//     fn fabs(_: libc::c_double) -> libc::c_double;
//     fn gettimeofday(__tv: *mut timeval, __tz: *mut libc::c_void) -> libc::c_int;
//     fn strlen(_: *const libc::c_char) -> libc::c_ulong;
// }
pub type __builtin_va_list = [__va_list_tag; 1];
#[derive(Copy, Clone)]
#[repr(C)]
pub struct __va_list_tag {
    pub gp_offset: libc::c_uint,
    pub fp_offset: libc::c_uint,
    pub overflow_arg_area: *mut libc::c_void,
    pub reg_save_area: *mut libc::c_void,
}
pub type size_t = libc::c_ulong;
pub type __uint32_t = libc::c_uint;
pub type __off_t = libc::c_long;
pub type __off64_t = libc::c_long;
pub type __time_t = libc::c_long;
pub type __suseconds_t = libc::c_long;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct timeval {
    pub tv_sec: __time_t,
    pub tv_usec: __suseconds_t,
}
pub type va_list = __builtin_va_list;
// #[derive(Copy, Clone)]
// #[repr(C)]
// pub struct _IO_FILE {
//     pub _flags: libc::c_int,
//     pub _IO_read_ptr: *mut libc::c_char,
//     pub _IO_read_end: *mut libc::c_char,
//     pub _IO_read_base: *mut libc::c_char,
//     pub _IO_write_base: *mut libc::c_char,
//     pub _IO_write_ptr: *mut libc::c_char,
//     pub _IO_write_end: *mut libc::c_char,
//     pub _IO_buf_base: *mut libc::c_char,
//     pub _IO_buf_end: *mut libc::c_char,
//     pub _IO_save_base: *mut libc::c_char,
//     pub _IO_backup_base: *mut libc::c_char,
//     pub _IO_save_end: *mut libc::c_char,
//     pub _markers: *mut _IO_marker,
//     pub _chain: *mut _IO_FILE,
//     pub _fileno: libc::c_int,
//     pub _flags2: libc::c_int,
//     pub _old_offset: __off_t,
//     pub _cur_column: libc::c_ushort,
//     pub _vtable_offset: libc::c_schar,
//     pub _shortbuf: [libc::c_char; 1],
//     pub _lock: *mut libc::c_void,
//     pub _offset: __off64_t,
//     pub _codecvt: *mut _IO_codecvt,
//     pub _wide_data: *mut _IO_wide_data,
//     pub _freeres_list: *mut _IO_FILE,
//     pub _freeres_buf: *mut libc::c_void,
//     pub __pad5: size_t,
//     pub _mode: libc::c_int,
//     pub _unused2: [libc::c_char; 20],
// }
// pub type _IO_lock_t = ();
// pub type FILE = _IO_FILE;

pub type nlopt_func = Option<
    fn(
        libc::c_uint,
        *const libc::c_double,
        *mut libc::c_double,
        *mut libc::c_void,
    ) -> libc::c_double,
>;
pub type nlopt_mfunc = Option<
    unsafe fn(
        libc::c_uint,
        *mut libc::c_double,
        libc::c_uint,
        *const libc::c_double,
        *mut libc::c_double,
        *mut libc::c_void,
    ) -> (),
>;
pub type nlopt_precond = Option<
    unsafe fn(
        libc::c_uint,
        *const libc::c_double,
        *const libc::c_double,
        *mut libc::c_double,
        *mut libc::c_void,
    ) -> (),
>;
pub type nlopt_result = libc::c_int;
pub const NLOPT_NUM_RESULTS: nlopt_result = 7;
pub const NLOPT_MAXTIME_REACHED: nlopt_result = 6;
pub const NLOPT_MAXEVAL_REACHED: nlopt_result = 5;
pub const NLOPT_XTOL_REACHED: nlopt_result = 4;
pub const NLOPT_FTOL_REACHED: nlopt_result = 3;
pub const NLOPT_STOPVAL_REACHED: nlopt_result = 2;
pub const NLOPT_SUCCESS: nlopt_result = 1;
pub const NLOPT_NUM_FAILURES: nlopt_result = -6;
pub const NLOPT_FORCED_STOP: nlopt_result = -5;
pub const NLOPT_ROUNDOFF_LIMITED: nlopt_result = -4;
pub const NLOPT_OUT_OF_MEMORY: nlopt_result = -3;
pub const NLOPT_INVALID_ARGS: nlopt_result = -2;
pub const NLOPT_FAILURE: nlopt_result = -1;
#[derive(Clone)]
#[repr(C)]
pub struct nlopt_stopping {
    pub n: libc::c_uint,
    pub minf_max: libc::c_double,
    pub ftol_rel: libc::c_double,
    pub ftol_abs: libc::c_double,
    pub xtol_rel: libc::c_double,
    pub xtol_abs: *const libc::c_double,
    pub x_weights: *const libc::c_double,
    pub nevals_p: *mut libc::c_int,
    pub maxeval: libc::c_int,
    pub maxtime: libc::c_double,
    pub start: libc::c_double,
    pub force_stop: *mut libc::c_int,
    pub stop_msg: String,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct nlopt_constraint {
    pub m: libc::c_uint,
    pub f: nlopt_func,
    pub mf: nlopt_mfunc,
    pub pre: nlopt_precond,
    pub f_data: *mut libc::c_void,
    pub tol: *mut libc::c_double,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct func_wrap_state {
    pub f: nlopt_func,
    pub f_data: *mut libc::c_void,
    pub m_orig: libc::c_uint,
    pub fc: *mut nlopt_constraint,
    pub p: libc::c_uint,
    pub h: *mut nlopt_constraint,
    pub xtmp: *mut libc::c_double,
    pub lb: *mut libc::c_double,
    pub ub: *mut libc::c_double,
    pub con_tol: *mut libc::c_double,
    pub scale: *mut libc::c_double,
    pub stop: *mut nlopt_stopping,
}
pub const COBYLA_MSG_NONE: C2RustUnnamed = 0;
pub type cobyla_function = unsafe fn(
    libc::c_int,
    libc::c_int,
    *mut libc::c_double,
    *mut libc::c_double,
    *mut libc::c_double,
    *mut func_wrap_state,
) -> libc::c_int;
pub type uint32_t = __uint32_t;
pub type C2RustUnnamed = libc::c_uint;
pub const COBYLA_MSG_INFO: C2RustUnnamed = 3;
pub const COBYLA_MSG_ITER: C2RustUnnamed = 2;
pub const COBYLA_MSG_EXIT: C2RustUnnamed = 1;
#[no_mangle]
pub unsafe fn nlopt_time_seed() -> libc::c_ulong {
    // let mut tv = libc::timeval {
    //     tv_sec: 0,
    //     tv_usec: 0,
    // };
    // libc::gettimeofday(&mut tv, 0 as *mut libc::timezone);
    //return (tv.tv_sec ^ tv.tv_usec) as libc::c_ulong;
    let start = SystemTime::now();
    let since_the_epoch = start.duration_since(UNIX_EPOCH).expect("Time flies");
    since_the_epoch.as_millis() as libc::c_ulong
}
#[no_mangle]
pub unsafe fn nlopt_seconds() -> libc::c_double {
    // static mut start_inited: libc::c_int = 0 as libc::c_int;
    // static mut start: libc::timeval = libc::timeval {
    //     tv_sec: 0,
    //     tv_usec: 0,
    // };
    // let mut tv: libc::timeval = libc::timeval {
    //     tv_sec: 0,
    //     tv_usec: 0,
    // };
    // if start_inited == 0 {
    //     start_inited = 1 as libc::c_int;
    //     libc::gettimeofday(&mut start, 0 as *mut libc::timezone);
    // }
    // libc::gettimeofday(&mut tv, 0 as *mut libc::timezone);
    // return (tv.tv_sec - start.tv_sec) as libc::c_double
    //     + 1.0e-6f64 * (tv.tv_usec - start.tv_usec) as libc::c_double;
    static mut start_inited: bool = false;
    static mut start: SystemTime = UNIX_EPOCH;
    if !start_inited {
        start_inited = true;
        start = SystemTime::now();
    }
    start
        .duration_since(UNIX_EPOCH)
        .expect("Time flies")
        .as_secs_f64()
}
unsafe fn sc(
    mut x: libc::c_double,
    mut smin: libc::c_double,
    mut smax: libc::c_double,
) -> libc::c_double {
    return smin + x * (smax - smin);
}
unsafe fn vector_norm(
    mut n: libc::c_uint,
    mut vec: *const libc::c_double,
    mut w: *const libc::c_double,
    mut scale_min: *const libc::c_double,
    mut scale_max: *const libc::c_double,
) -> libc::c_double {
    let mut i: libc::c_uint = 0;
    let mut ret: libc::c_double = 0 as libc::c_int as libc::c_double;
    if !scale_min.is_null() && !scale_max.is_null() {
        if !w.is_null() {
            i = 0 as libc::c_int as libc::c_uint;
            while i < n {
                ret += *w.offset(i as isize)
                    * (sc(
                        *vec.offset(i as isize),
                        *scale_min.offset(i as isize),
                        *scale_max.offset(i as isize),
                    )
                    .abs());
                i = i.wrapping_add(1);
            }
        } else {
            i = 0 as libc::c_int as libc::c_uint;
            while i < n {
                ret += (sc(
                    *vec.offset(i as isize),
                    *scale_min.offset(i as isize),
                    *scale_max.offset(i as isize),
                ))
                .abs();
                i = i.wrapping_add(1);
            }
        }
    } else if !w.is_null() {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            ret += *w.offset(i as isize) * (*vec.offset(i as isize)).abs();
            i = i.wrapping_add(1);
        }
    } else {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            ret += (*vec.offset(i as isize)).abs();
            i = i.wrapping_add(1);
        }
    }
    return ret;
}
unsafe fn diff_norm(
    mut n: libc::c_uint,
    mut x: *const libc::c_double,
    mut oldx: *const libc::c_double,
    mut w: *const libc::c_double,
    mut scale_min: *const libc::c_double,
    mut scale_max: *const libc::c_double,
) -> libc::c_double {
    let mut i: libc::c_uint = 0;
    let mut ret: libc::c_double = 0 as libc::c_int as libc::c_double;
    if !scale_min.is_null() && !scale_max.is_null() {
        if !w.is_null() {
            i = 0 as libc::c_int as libc::c_uint;
            while i < n {
                ret += *w.offset(i as isize)
                    * (sc(
                        *x.offset(i as isize),
                        *scale_min.offset(i as isize),
                        *scale_max.offset(i as isize),
                    ) - sc(
                        *oldx.offset(i as isize),
                        *scale_min.offset(i as isize),
                        *scale_max.offset(i as isize),
                    ))
                    .abs();
                i = i.wrapping_add(1);
            }
        } else {
            i = 0 as libc::c_int as libc::c_uint;
            while i < n {
                ret += (sc(
                    *x.offset(i as isize),
                    *scale_min.offset(i as isize),
                    *scale_max.offset(i as isize),
                ) - sc(
                    *oldx.offset(i as isize),
                    *scale_min.offset(i as isize),
                    *scale_max.offset(i as isize),
                ))
                .abs();
                i = i.wrapping_add(1);
            }
        }
    } else if !w.is_null() {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            ret += *w.offset(i as isize) * (*x.offset(i as isize) - *oldx.offset(i as isize)).abs();
            i = i.wrapping_add(1);
        }
    } else {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            ret += (*x.offset(i as isize) - *oldx.offset(i as isize)).abs();
            i = i.wrapping_add(1);
        }
    }
    return ret;
}
unsafe fn relstop(
    mut vold: libc::c_double,
    mut vnew: libc::c_double,
    mut reltol: libc::c_double,
    mut abstol: libc::c_double,
) -> libc::c_int {
    if nlopt_isinf(vold) != 0 {
        return 0 as libc::c_int;
    }
    return ((vnew - vold).abs() < abstol
        || (vnew - vold).abs() < reltol * ((vnew).abs() + (vold)).abs() * 0.5f64
        || reltol > 0 as libc::c_int as libc::c_double && vnew == vold) as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_ftol(
    mut s: *const nlopt_stopping,
    mut f: libc::c_double,
    mut oldf: libc::c_double,
) -> libc::c_int {
    return relstop(oldf, f, (*s).ftol_rel, (*s).ftol_abs);
}
#[no_mangle]
pub unsafe fn nlopt_stop_f(
    mut s: *const nlopt_stopping,
    mut f: libc::c_double,
    mut oldf: libc::c_double,
) -> libc::c_int {
    return (f <= (*s).minf_max || nlopt_stop_ftol(s, f, oldf) != 0) as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_x(
    mut s: *const nlopt_stopping,
    mut x: *const libc::c_double,
    mut oldx: *const libc::c_double,
) -> libc::c_int {
    let mut i: libc::c_uint = 0;
    if diff_norm(
        (*s).n,
        x,
        oldx,
        (*s).x_weights,
        0 as *const libc::c_double,
        0 as *const libc::c_double,
    ) < (*s).xtol_rel
        * vector_norm(
            (*s).n,
            x,
            (*s).x_weights,
            0 as *const libc::c_double,
            0 as *const libc::c_double,
        )
    {
        return 1 as libc::c_int;
    }
    if ((*s).xtol_abs).is_null() {
        return 0 as libc::c_int;
    }
    i = 0 as libc::c_int as libc::c_uint;
    while i < (*s).n {
        if (*x.offset(i as isize) - *oldx.offset(i as isize)).abs()
            >= *((*s).xtol_abs).offset(i as isize)
        {
            return 0 as libc::c_int;
        }
        i = i.wrapping_add(1);
    }
    return 1 as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_dx(
    mut s: *const nlopt_stopping,
    mut x: *const libc::c_double,
    mut dx: *const libc::c_double,
) -> libc::c_int {
    let mut i: libc::c_uint = 0;
    if vector_norm(
        (*s).n,
        dx,
        (*s).x_weights,
        0 as *const libc::c_double,
        0 as *const libc::c_double,
    ) < (*s).xtol_rel
        * vector_norm(
            (*s).n,
            x,
            (*s).x_weights,
            0 as *const libc::c_double,
            0 as *const libc::c_double,
        )
    {
        return 1 as libc::c_int;
    }
    if ((*s).xtol_abs).is_null() {
        return 0 as libc::c_int;
    }
    i = 0 as libc::c_int as libc::c_uint;
    while i < (*s).n {
        if (*dx.offset(i as isize)).abs() >= *((*s).xtol_abs).offset(i as isize) {
            return 0 as libc::c_int;
        }
        i = i.wrapping_add(1);
    }
    return 1 as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_xs(
    mut s: *const nlopt_stopping,
    mut xs: *const libc::c_double,
    mut oldxs: *const libc::c_double,
    mut scale_min: *const libc::c_double,
    mut scale_max: *const libc::c_double,
) -> libc::c_int {
    let mut i: libc::c_uint = 0;
    if diff_norm((*s).n, xs, oldxs, (*s).x_weights, scale_min, scale_max)
        < (*s).xtol_rel * vector_norm((*s).n, xs, (*s).x_weights, scale_min, scale_max)
    {
        return 1 as libc::c_int;
    }
    if ((*s).xtol_abs).is_null() {
        return 0 as libc::c_int;
    }
    i = 0 as libc::c_int as libc::c_uint;
    while i < (*s).n {
        if (sc(
            *xs.offset(i as isize),
            *scale_min.offset(i as isize),
            *scale_max.offset(i as isize),
        ) - sc(
            *oldxs.offset(i as isize),
            *scale_min.offset(i as isize),
            *scale_max.offset(i as isize),
        ))
        .abs()
            >= *((*s).xtol_abs).offset(i as isize)
        {
            return 0 as libc::c_int;
        }
        i = i.wrapping_add(1);
    }
    return 1 as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_isinf(mut x: libc::c_double) -> libc::c_int {
    return ((x).abs() >= ::std::f64::INFINITY * 0.99f64
        || if x.is_infinite() {
            if x.is_sign_positive() {
                1
            } else {
                -1
            }
        } else {
            0
        } != 0) as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_isfinite(mut x: libc::c_double) -> libc::c_int {
    return (x.abs() <= 1.7976931348623157e+308f64) as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_istiny(mut x: libc::c_double) -> libc::c_int {
    if x == 0.0f64 {
        return 1 as libc::c_int;
    } else {
        return (x.abs() < 2.2250738585072014e-308f64) as libc::c_int;
    };
}
#[no_mangle]
pub unsafe fn nlopt_isnan(mut x: libc::c_double) -> libc::c_int {
    return x.is_nan() as i32;
}
#[no_mangle]
pub unsafe fn nlopt_stop_evals(mut s: *const nlopt_stopping) -> libc::c_int {
    return ((*s).maxeval > 0 as libc::c_int && *(*s).nevals_p >= (*s).maxeval) as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_time_(
    mut start: libc::c_double,
    mut maxtime: libc::c_double,
) -> libc::c_int {
    return (maxtime > 0 as libc::c_int as libc::c_double && nlopt_seconds() - start >= maxtime)
        as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_time(mut s: *const nlopt_stopping) -> libc::c_int {
    return nlopt_stop_time_((*s).start, (*s).maxtime);
}
#[no_mangle]
pub unsafe fn nlopt_stop_evalstime(mut stop: *const nlopt_stopping) -> libc::c_int {
    return (nlopt_stop_evals(stop) != 0 || nlopt_stop_time(stop) != 0) as libc::c_int;
}
#[no_mangle]
pub unsafe fn nlopt_stop_forced(mut stop: *const nlopt_stopping) -> libc::c_int {
    return (!((*stop).force_stop).is_null() && *(*stop).force_stop != 0) as libc::c_int;
}
// #[no_mangle]
// pub unsafe  fn nlopt_vsprintf(
//     mut p: *mut libc::c_char,
//     mut format: *const libc::c_char,
//     mut ap: ::std::ffi::VaList,
// ) -> *mut libc::c_char {
//     let mut len: size_t = (strlen(format)).wrapping_add(128 as libc::c_int as libc::c_ulong);
//     let mut ret: libc::c_int = 0;
//     p = realloc(p as *mut libc::c_void, len) as *mut libc::c_char;
//     if p.is_null() {
//         abort();
//     }
//     loop {
//         ret = vsnprintf(p, len, format, ap.as_va_list());
//         if !(ret < 0 as libc::c_int || ret as size_t >= len) {
//             break;
//         }
//         len = if ret >= 0 as libc::c_int {
//             (ret + 1 as libc::c_int) as size_t
//         } else {
//             len.wrapping_mul(3 as libc::c_int as libc::c_ulong) >> 1 as libc::c_int
//         };
//         p = realloc(p as *mut libc::c_void, len) as *mut libc::c_char;
//         if p.is_null() {
//             abort();
//         }
//     }
//     return p;
// }
#[no_mangle]
pub unsafe fn nlopt_stop_msg(mut s: *mut nlopt_stopping, msg: &str) {
    (*s).stop_msg = msg.to_string();
}
#[no_mangle]
pub unsafe fn nlopt_count_constraints(
    mut p: libc::c_uint,
    mut c: *const nlopt_constraint,
) -> libc::c_uint {
    let mut i: libc::c_uint = 0;
    let mut count: libc::c_uint = 0 as libc::c_int as libc::c_uint;
    i = 0 as libc::c_int as libc::c_uint;
    while i < p {
        count = count.wrapping_add((*c.offset(i as isize)).m);
        i = i.wrapping_add(1);
    }
    return count;
}
#[no_mangle]
pub unsafe fn nlopt_max_constraint_dim(
    mut p: libc::c_uint,
    mut c: *const nlopt_constraint,
) -> libc::c_uint {
    let mut i: libc::c_uint = 0;
    let mut max_dim: libc::c_uint = 0 as libc::c_int as libc::c_uint;
    i = 0 as libc::c_int as libc::c_uint;
    while i < p {
        if (*c.offset(i as isize)).m > max_dim {
            max_dim = (*c.offset(i as isize)).m;
        }
        i = i.wrapping_add(1);
    }
    return max_dim;
}

pub unsafe fn nlopt_eval_constraint<U>(
    mut result: *mut libc::c_double,
    mut grad: *mut libc::c_double,
    mut c: *const nlopt_constraint,
    mut n: libc::c_uint,
    mut x: *const libc::c_double,
) {
    if ((*c).f).is_some() {
        *result.offset(0 as libc::c_int as isize) =
        // PATCH Weird bug ((*c).f).expect("non-null function pointer") calls the objective function!!!
        // even if (*c), nlopt_constraint object was correctly built with a nlopt_constraint_raw_callback!!! 
        //    ((*c).f).expect("non-null function pointer")(n, x, grad, (*c).f_data);
        // Maybe the U generic parameter required explains it cannot work like with C ???
        nlopt_constraint_raw_callback::<&dyn Func<U>, U>(n, x, grad, (*c).f_data);
        // RLA: Take the opposite to manage cstr as being nonnegative in the end like the original cobyla
        *result.offset(0 as libc::c_int as isize) = -*result.offset(0 as libc::c_int as isize)
    } else {
        ((*c).mf).expect("non-null function pointer")((*c).m, result, n, x, grad, (*c).f_data);
    };
}
#[no_mangle]
pub unsafe fn nlopt_compute_rescaling(
    mut n: libc::c_uint,
    mut dx: *const libc::c_double,
) -> *mut libc::c_double {
    // let mut s: *mut libc::c_double = malloc(
    //     (::std::mem::size_of::<libc::c_double>() as libc::c_ulong).wrapping_mul(n as libc::c_ulong),
    // ) as *mut libc::c_double;

    let mut space: Box<Vec<libc::c_double>> = Box::new(vec![0.; usize::try_from(n).unwrap()]);
    let s = space.as_mut_ptr() as *mut libc::c_double;
    std::mem::forget(space);

    let mut i: libc::c_uint = 0;
    if s.is_null() {
        return 0 as *mut libc::c_double;
    }
    i = 0 as libc::c_int as libc::c_uint;
    while i < n {
        *s.offset(i as isize) = 1.0f64;
        i = i.wrapping_add(1);
    }
    if n == 1 as libc::c_int as libc::c_uint {
        return s;
    }
    i = 1 as libc::c_int as libc::c_uint;
    while i < n
        && *dx.offset(i as isize)
            == *dx.offset(i.wrapping_sub(1 as libc::c_int as libc::c_uint) as isize)
    {
        i = i.wrapping_add(1);
    }
    if i < n {
        i = 1 as libc::c_int as libc::c_uint;
        while i < n {
            *s.offset(i as isize) = *dx.offset(i as isize) / *dx.offset(0 as libc::c_int as isize);
            i = i.wrapping_add(1);
        }
    }
    return s;
}
#[no_mangle]
pub unsafe fn nlopt_rescale(
    mut n: libc::c_uint,
    mut s: *const libc::c_double,
    mut x: *const libc::c_double,
    mut xs: *mut libc::c_double,
) {
    let mut i: libc::c_uint = 0;
    if s.is_null() {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            *xs.offset(i as isize) = *x.offset(i as isize);
            i = i.wrapping_add(1);
        }
    } else {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            *xs.offset(i as isize) = *x.offset(i as isize) / *s.offset(i as isize);
            i = i.wrapping_add(1);
        }
    };
}
#[no_mangle]
pub unsafe fn nlopt_unscale(
    mut n: libc::c_uint,
    mut s: *const libc::c_double,
    mut x: *const libc::c_double,
    mut xs: *mut libc::c_double,
) {
    let mut i: libc::c_uint = 0;
    if s.is_null() {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            *xs.offset(i as isize) = *x.offset(i as isize);
            i = i.wrapping_add(1);
        }
    } else {
        i = 0 as libc::c_int as libc::c_uint;
        while i < n {
            *xs.offset(i as isize) = *x.offset(i as isize) * *s.offset(i as isize);
            i = i.wrapping_add(1);
        }
    };
}
#[no_mangle]
pub unsafe fn nlopt_new_rescaled(
    mut n: libc::c_uint,
    mut s: *const libc::c_double,
    mut x: *const libc::c_double,
) -> *mut libc::c_double {
    // let mut xs: *mut libc::c_double = malloc(
    //     (::std::mem::size_of::<libc::c_double>() as libc::c_ulong).wrapping_mul(n as libc::c_ulong),
    // ) as *mut libc::c_double;

    let mut space: Box<Vec<libc::c_double>> = Box::new(vec![0.; usize::try_from(n).unwrap()]);
    let xs = space.as_mut_ptr() as *mut libc::c_double;
    std::mem::forget(space);

    if xs.is_null() {
        return 0 as *mut libc::c_double;
    }
    nlopt_rescale(n, s, x, xs);
    return xs;
}
#[no_mangle]
pub unsafe fn nlopt_reorder_bounds(
    mut n: libc::c_uint,
    mut lb: *mut libc::c_double,
    mut ub: *mut libc::c_double,
) {
    let mut i: libc::c_uint = 0;
    i = 0 as libc::c_int as libc::c_uint;
    while i < n {
        if *lb.offset(i as isize) > *ub.offset(i as isize) {
            let mut t: libc::c_double = *lb.offset(i as isize);
            *lb.offset(i as isize) = *ub.offset(i as isize);
            *ub.offset(i as isize) = t;
        }
        i = i.wrapping_add(1);
    }
}
unsafe fn func_wrap<U>(
    mut ni: libc::c_int,
    mut _mi: libc::c_int,
    mut x: *mut libc::c_double,
    mut f: *mut libc::c_double,
    mut con: *mut libc::c_double,
    mut s: *mut func_wrap_state,
) -> libc::c_int {
    let mut n: libc::c_uint = ni as libc::c_uint;
    let mut i: libc::c_uint = 0;
    let mut j: libc::c_uint = 0;
    let mut k: libc::c_uint = 0;
    let mut xtmp: *mut libc::c_double = (*s).xtmp;
    let mut lb: *const libc::c_double = (*s).lb;
    let mut ub: *const libc::c_double = (*s).ub;
    j = 0 as libc::c_int as libc::c_uint;
    while j < n {
        if *x.offset(j as isize) < *lb.offset(j as isize) {
            *xtmp.offset(j as isize) = *lb.offset(j as isize);
        } else if *x.offset(j as isize) > *ub.offset(j as isize) {
            *xtmp.offset(j as isize) = *ub.offset(j as isize);
        } else {
            *xtmp.offset(j as isize) = *x.offset(j as isize);
        }
        j = j.wrapping_add(1);
    }
    nlopt_unscale(n, (*s).scale, xtmp, xtmp);
    *f = ((*s).f).expect("non-null function pointer")(
        n,
        xtmp,
        0 as *mut libc::c_double,
        (*s).f_data,
    );
    if nlopt_stop_forced((*s).stop) != 0 {
        return 1 as libc::c_int;
    }
    i = 0 as libc::c_int as libc::c_uint;
    j = 0 as libc::c_int as libc::c_uint;
    while j < (*s).m_orig {
        nlopt_eval_constraint::<U>(
            con.offset(i as isize),
            0 as *mut libc::c_double,
            ((*s).fc).offset(j as isize),
            n,
            xtmp,
        );
        if nlopt_stop_forced((*s).stop) != 0 {
            return 1 as libc::c_int;
        }
        k = 0 as libc::c_int as libc::c_uint;
        while k < (*((*s).fc).offset(j as isize)).m {
            *con.offset(i.wrapping_add(k) as isize) = -*con.offset(i.wrapping_add(k) as isize);
            k = k.wrapping_add(1);
        }
        i = i.wrapping_add((*((*s).fc).offset(j as isize)).m);
        j = j.wrapping_add(1);
    }
    j = 0 as libc::c_int as libc::c_uint;
    while j < (*s).p {
        nlopt_eval_constraint::<U>(
            con.offset(i as isize),
            0 as *mut libc::c_double,
            ((*s).h).offset(j as isize),
            n,
            xtmp,
        );
        if nlopt_stop_forced((*s).stop) != 0 {
            return 1 as libc::c_int;
        }
        k = 0 as libc::c_int as libc::c_uint;
        while k < (*((*s).h).offset(j as isize)).m {
            *con.offset(
                i.wrapping_add((*((*s).h).offset(j as isize)).m)
                    .wrapping_add(k) as isize,
            ) = -*con.offset(i.wrapping_add(k) as isize);
            k = k.wrapping_add(1);
        }
        i = i.wrapping_add(
            (2 as libc::c_int as libc::c_uint).wrapping_mul((*((*s).h).offset(j as isize)).m),
        );
        j = j.wrapping_add(1);
    }
    j = 0 as libc::c_int as libc::c_uint;
    while j < n {
        if nlopt_isinf(*lb.offset(j as isize)) == 0 {
            let fresh1 = i;
            i = i.wrapping_add(1);
            *con.offset(fresh1 as isize) = *x.offset(j as isize) - *lb.offset(j as isize);
        }
        if nlopt_isinf(*ub.offset(j as isize)) == 0 {
            let fresh2 = i;
            i = i.wrapping_add(1);
            *con.offset(fresh2 as isize) = *ub.offset(j as isize) - *x.offset(j as isize);
        }
        j = j.wrapping_add(1);
    }
    return 0 as libc::c_int;
}
pub unsafe fn cobyla_minimize<U>(
    mut n: libc::c_uint,
    mut f: nlopt_func,
    mut f_data: *mut libc::c_void,
    mut m: libc::c_uint,
    mut fc: *mut nlopt_constraint,
    mut p: libc::c_uint,
    mut h: *mut nlopt_constraint,
    mut lb: *const libc::c_double,
    mut ub: *const libc::c_double,
    mut x: *mut libc::c_double,
    mut minf: *mut libc::c_double,
    mut stop: *mut nlopt_stopping,
    mut dx: *const libc::c_double,
) -> nlopt_result {
    let mut current_block: u64;
    let mut i: libc::c_uint = 0;
    let mut j: libc::c_uint = 0;
    let mut s: func_wrap_state = func_wrap_state {
        f: None,
        f_data: 0 as *mut libc::c_void,
        m_orig: 0,
        fc: 0 as *mut nlopt_constraint,
        p: 0,
        h: 0 as *mut nlopt_constraint,
        xtmp: 0 as *mut libc::c_double,
        lb: 0 as *mut libc::c_double,
        ub: 0 as *mut libc::c_double,
        con_tol: 0 as *mut libc::c_double,
        scale: 0 as *mut libc::c_double,
        stop: 0 as *mut nlopt_stopping,
    };
    let mut ret: nlopt_result = 0 as nlopt_result;
    let mut rhobeg: libc::c_double = 0.;
    let mut rhoend: libc::c_double = 0.;
    s.f = f;
    s.f_data = f_data;
    s.m_orig = m;
    s.fc = fc;
    s.p = p;
    s.h = h;
    s.stop = stop;
    s.scale = 0 as *mut libc::c_double;
    s.con_tol = s.scale;
    s.xtmp = s.con_tol;
    s.ub = s.xtmp;
    s.lb = s.ub;
    s.scale = nlopt_compute_rescaling(n, dx);
    if (s.scale).is_null() {
        ret = NLOPT_OUT_OF_MEMORY;
    } else {
        j = 0 as libc::c_int as libc::c_uint;
        loop {
            if !(j < n) {
                current_block = 15652330335145281839;
                break;
            }
            if *(s.scale).offset(j as isize) == 0 as libc::c_int as libc::c_double
                || nlopt_isfinite(*(s.scale).offset(j as isize)) == 0
            {
                nlopt_stop_msg(
                    stop,
                    &format!(
                        "invalid scaling {} of dimension {}: possible over/underflow?",
                        *(s.scale).offset(j as isize),
                        j
                    ),
                );
                ret = NLOPT_INVALID_ARGS;
                current_block = 762786280471819104;
                break;
            } else {
                j = j.wrapping_add(1);
            }
        }
        match current_block {
            762786280471819104 => {}
            _ => {
                s.lb = nlopt_new_rescaled(n, s.scale, lb);
                if (s.lb).is_null() {
                    ret = NLOPT_OUT_OF_MEMORY;
                } else {
                    s.ub = nlopt_new_rescaled(n, s.scale, ub);
                    if (s.ub).is_null() {
                        ret = NLOPT_OUT_OF_MEMORY;
                    } else {
                        nlopt_reorder_bounds(n, s.lb, s.ub);
                        // s.xtmp = malloc(
                        //     (::std::mem::size_of::<libc::c_double>() as libc::c_ulong)
                        //         .wrapping_mul(n as libc::c_ulong),
                        // ) as *mut libc::c_double;

                        let mut space: Box<Vec<libc::c_double>> =
                            Box::new(vec![0.; usize::try_from(n).unwrap()]);
                        s.xtmp = space.as_mut_ptr() as *mut libc::c_double;
                        std::mem::forget(space);

                        if (s.xtmp).is_null() {
                            ret = NLOPT_OUT_OF_MEMORY;
                        } else {
                            rhobeg = (*dx.offset(0 as libc::c_int as isize)
                                / *(s.scale).offset(0 as libc::c_int as isize))
                            .abs();
                            rhoend = (*stop).xtol_rel * rhobeg;
                            if !((*stop).xtol_abs).is_null() {
                                j = 0 as libc::c_int as libc::c_uint;
                                while j < n {
                                    if rhoend
                                        < *((*stop).xtol_abs).offset(j as isize)
                                            / (*(s.scale).offset(j as isize)).abs()
                                    {
                                        rhoend = *((*stop).xtol_abs).offset(j as isize)
                                            / (*(s.scale).offset(j as isize)).abs();
                                    }
                                    j = j.wrapping_add(1);
                                }
                            }
                            m = (nlopt_count_constraints(m, fc)).wrapping_add(
                                (2 as libc::c_int as libc::c_uint)
                                    .wrapping_mul(nlopt_count_constraints(p, h)),
                            );
                            j = 0 as libc::c_int as libc::c_uint;
                            while j < n {
                                if nlopt_isinf(*lb.offset(j as isize)) == 0 {
                                    m = m.wrapping_add(1);
                                }
                                if nlopt_isinf(*ub.offset(j as isize)) == 0 {
                                    m = m.wrapping_add(1);
                                }
                                j = j.wrapping_add(1);
                            }
                            // s.con_tol = malloc(
                            //     (::std::mem::size_of::<libc::c_double>() as libc::c_ulong)
                            //         .wrapping_mul(m as libc::c_ulong),
                            // ) as *mut libc::c_double;

                            if m > 0 {
                                let mut space: Box<Vec<libc::c_double>> =
                                    Box::new(vec![0.; usize::try_from(m).unwrap()]);
                                s.con_tol = space.as_mut_ptr() as *mut libc::c_double;
                                std::mem::forget(space);
                            }

                            if m != 0 && (s.con_tol).is_null() {
                                ret = NLOPT_OUT_OF_MEMORY;
                            } else {
                                j = 0 as libc::c_int as libc::c_uint;
                                while j < m {
                                    *(s.con_tol).offset(j as isize) =
                                        0 as libc::c_int as libc::c_double;
                                    j = j.wrapping_add(1);
                                }
                                i = 0 as libc::c_int as libc::c_uint;
                                j = i;
                                while i < s.m_orig {
                                    let mut ji: libc::c_uint = j;
                                    let mut jnext: libc::c_uint =
                                        j.wrapping_add((*fc.offset(i as isize)).m);
                                    while j < jnext {
                                        *(s.con_tol).offset(j as isize) =
                                            *((*fc.offset(i as isize)).tol)
                                                .offset(j.wrapping_sub(ji) as isize);
                                        j = j.wrapping_add(1);
                                    }
                                    i = i.wrapping_add(1);
                                }
                                i = 0 as libc::c_int as libc::c_uint;
                                while i < s.p {
                                    let mut ji_0: libc::c_uint = j;
                                    let mut jnext_0: libc::c_uint =
                                        j.wrapping_add((*h.offset(i as isize)).m);
                                    while j < jnext_0 {
                                        *(s.con_tol).offset(j as isize) =
                                            *((*h.offset(i as isize)).tol)
                                                .offset(j.wrapping_sub(ji_0) as isize);
                                        j = j.wrapping_add(1);
                                    }
                                    ji_0 = j;
                                    jnext_0 = j.wrapping_add((*h.offset(i as isize)).m);
                                    while j < jnext_0 {
                                        *(s.con_tol).offset(j as isize) =
                                            *((*h.offset(i as isize)).tol)
                                                .offset(j.wrapping_sub(ji_0) as isize);
                                        j = j.wrapping_add(1);
                                    }
                                    i = i.wrapping_add(1);
                                }
                                nlopt_rescale(n, s.scale, x, x);
                                ret = cobyla(
                                    n as libc::c_int,
                                    m as libc::c_int,
                                    x,
                                    minf,
                                    rhobeg,
                                    rhoend,
                                    stop,
                                    s.lb,
                                    s.ub,
                                    COBYLA_MSG_NONE as libc::c_int,
                                    Some(
                                        func_wrap::<U>
                                            as unsafe fn(
                                                libc::c_int,
                                                libc::c_int,
                                                *mut libc::c_double,
                                                *mut libc::c_double,
                                                *mut libc::c_double,
                                                *mut func_wrap_state,
                                            )
                                                -> libc::c_int,
                                    ),
                                    &mut s,
                                );
                                nlopt_unscale(n, s.scale, x, x);
                                j = 0 as libc::c_int as libc::c_uint;
                                while j < n {
                                    if *x.offset(j as isize) < *lb.offset(j as isize) {
                                        *x.offset(j as isize) = *lb.offset(j as isize);
                                    }
                                    if *x.offset(j as isize) > *ub.offset(j as isize) {
                                        *x.offset(j as isize) = *ub.offset(j as isize);
                                    }
                                    j = j.wrapping_add(1);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // free(s.con_tol as *mut libc::c_void);
    // free(s.xtmp as *mut libc::c_void);
    // free(s.ub as *mut libc::c_void);
    // free(s.lb as *mut libc::c_void);
    // free(s.scale as *mut libc::c_void);

    if m > 0 {
        let _ = Box::from_raw(s.con_tol);
    }
    let _ = Box::from_raw(s.xtmp);
    let _ = Box::from_raw(s.ub);
    let _ = Box::from_raw(s.lb);
    let _ = Box::from_raw(s.scale);

    return ret;
}
unsafe fn lcg_rand(mut seed: *mut uint32_t) -> uint32_t {
    *seed = (*seed)
        .wrapping_mul(1103515245 as libc::c_int as libc::c_uint)
        .wrapping_add(12345 as libc::c_int as libc::c_uint);
    return *seed;
}
unsafe fn lcg_urand(
    mut seed: *mut uint32_t,
    mut a: libc::c_double,
    mut b: libc::c_double,
) -> libc::c_double {
    return a + lcg_rand(seed) as libc::c_double * (b - a)
        / -(1 as libc::c_int) as uint32_t as libc::c_double;
}
#[no_mangle]
pub unsafe fn cobyla(
    mut n: libc::c_int,
    mut m: libc::c_int,
    mut x: *mut libc::c_double,
    mut minf: *mut libc::c_double,
    mut rhobeg: libc::c_double,
    mut rhoend: libc::c_double,
    mut stop: *mut nlopt_stopping,
    mut lb: *const libc::c_double,
    mut ub: *const libc::c_double,
    mut iprint: libc::c_int,
    mut calcfc: Option<cobyla_function>,
    mut state: *mut func_wrap_state,
) -> nlopt_result {
    let mut icon: libc::c_int = 0;
    let mut isim: libc::c_int = 0;
    let mut isigb: libc::c_int = 0;
    let mut idatm: libc::c_int = 0;
    let mut iveta: libc::c_int = 0;
    let mut isimi: libc::c_int = 0;
    let mut ivsig: libc::c_int = 0;
    let mut iwork: libc::c_int = 0;
    let mut ia: libc::c_int = 0;
    let mut idx: libc::c_int = 0;
    let mut mpp: libc::c_int = 0;
    let mut _iact: *mut libc::c_int = 0 as *mut libc::c_int;
    let mut _w: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut rc: nlopt_result = 0 as nlopt_result;
    *(*stop).nevals_p = 0 as libc::c_int;
    if n == 0 as libc::c_int {
        if iprint >= 1 as libc::c_int {
            fprintf(Io::stderr, "cobyla: N==0.");
        }
        return NLOPT_SUCCESS;
    }
    if n < 0 as libc::c_int || m < 0 as libc::c_int {
        if iprint >= 1 as libc::c_int {
            fprintf(Io::stderr, "cobyla: N<0 or M<0.");
        }
        return NLOPT_INVALID_ARGS;
    }
    // w = malloc(
    //     ((n * (3 as libc::c_int * n + 2 as libc::c_int * m + 11 as libc::c_int)
    //         + 4 as libc::c_int * m
    //         + 6 as libc::c_int) as libc::c_uint as libc::c_ulong)
    //         .wrapping_mul(::std::mem::size_of::<libc::c_double>() as libc::c_ulong),
    // ) as *mut libc::c_double;

    let space_size = n * (3 as libc::c_int * n + 2 as libc::c_int * m + 11 as libc::c_int)
        + 4 as libc::c_int * m
        + 6 as libc::c_int;
    let mut space: Box<Vec<libc::c_double>> =
        Box::new(vec![0.; usize::try_from(space_size).unwrap()]);
    let mut w = space.as_mut_ptr() as *mut libc::c_double;
    std::mem::forget(space);

    if w.is_null() {
        if iprint >= 1 as libc::c_int {
            fprintf(Io::stderr, "cobyla: memory allocation error");
        }
        return NLOPT_OUT_OF_MEMORY;
    }
    // iact = malloc(
    //     ((m + 1 as libc::c_int) as libc::c_uint as libc::c_ulong)
    //         .wrapping_mul(::std::mem::size_of::<libc::c_int>() as libc::c_ulong),
    // ) as *mut libc::c_int;

    let space_size = m + 1;
    let mut space: Box<Vec<libc::c_double>> =
        Box::new(vec![0.; usize::try_from(space_size).unwrap()]);
    let mut iact = space.as_mut_ptr() as *mut libc::c_int;
    std::mem::forget(space);

    if iact.is_null() {
        if iprint >= 1 as libc::c_int {
            fprintf(Io::stderr, "cobyla: memory allocation error.");
        }
        //free(w as *mut libc::c_void);
        let _ = Box::from_raw(w);
        return NLOPT_OUT_OF_MEMORY;
    }
    iact = iact.offset(-1);
    w = w.offset(-1);
    x = x.offset(-1);
    lb = lb.offset(-1);
    ub = ub.offset(-1);
    mpp = m + 2 as libc::c_int;
    icon = 1 as libc::c_int;
    isim = icon + mpp;
    isimi = isim + n * n + n;
    idatm = isimi + n * n;
    ia = idatm + n * mpp + mpp;
    ivsig = ia + m * n + n;
    iveta = ivsig + n;
    isigb = iveta + n;
    idx = isigb + n;
    iwork = idx + n;
    rc = cobylb(
        &mut n,
        &mut m,
        &mut mpp,
        &mut *x.offset(1 as libc::c_int as isize),
        minf,
        &mut rhobeg,
        rhoend,
        stop,
        &*lb.offset(1 as libc::c_int as isize),
        &*ub.offset(1 as libc::c_int as isize),
        &mut iprint,
        &mut *w.offset(icon as isize),
        &mut *w.offset(isim as isize),
        &mut *w.offset(isimi as isize),
        &mut *w.offset(idatm as isize),
        &mut *w.offset(ia as isize),
        &mut *w.offset(ivsig as isize),
        &mut *w.offset(iveta as isize),
        &mut *w.offset(isigb as isize),
        &mut *w.offset(idx as isize),
        &mut *w.offset(iwork as isize),
        &mut *iact.offset(1 as libc::c_int as isize),
        calcfc,
        state,
    );
    iact = iact.offset(1);
    w = w.offset(1);
    // free(w as *mut libc::c_void);
    // free(iact as *mut libc::c_void);
    let _ = Box::from_raw(w);
    let _ = Box::from_raw(iact);
    return rc;
}
unsafe fn cobylb(
    mut n: *mut libc::c_int,
    mut m: *mut libc::c_int,
    mut mpp: *mut libc::c_int,
    mut x: *mut libc::c_double,
    mut minf: *mut libc::c_double,
    mut rhobeg: *mut libc::c_double,
    mut rhoend: libc::c_double,
    mut stop: *mut nlopt_stopping,
    mut lb: *const libc::c_double,
    mut ub: *const libc::c_double,
    mut iprint: *mut libc::c_int,
    mut con: *mut libc::c_double,
    mut sim: *mut libc::c_double,
    mut simi: *mut libc::c_double,
    mut datmat: *mut libc::c_double,
    mut a: *mut libc::c_double,
    mut vsig: *mut libc::c_double,
    mut veta: *mut libc::c_double,
    mut sigbar: *mut libc::c_double,
    mut dx: *mut libc::c_double,
    mut w: *mut libc::c_double,
    mut iact: *mut libc::c_int,
    mut calcfc: Option<cobyla_function>,
    mut state: *mut func_wrap_state,
) -> nlopt_result {
    let mut current_block: u64;
    let mut sim_dim1: libc::c_int = 0;
    let mut sim_offset: libc::c_int = 0;
    let mut simi_dim1: libc::c_int = 0;
    let mut simi_offset: libc::c_int = 0;
    let mut datmat_dim1: libc::c_int = 0;
    let mut datmat_offset: libc::c_int = 0;
    let mut a_dim1: libc::c_int = 0;
    let mut a_offset: libc::c_int = 0;
    let mut i__1: libc::c_int = 0;
    let mut i__2: libc::c_int = 0;
    let mut i__3: libc::c_int = 0;
    let mut d__1: libc::c_double = 0.;
    let mut d__2: libc::c_double = 0.;
    let mut alpha: libc::c_double = 0.;
    let mut delta: libc::c_double = 0.;
    let mut denom: libc::c_double = 0.;
    let mut tempa: libc::c_double = 0.;
    let mut barmu: libc::c_double = 0.;
    let mut beta: libc::c_double = 0.;
    let mut cmin: libc::c_double = 0.0f64;
    let mut cmax: libc::c_double = 0.0f64;
    let mut cvmaxm: libc::c_double = 0.;
    let mut dxsign: libc::c_double = 0.;
    let mut prerem: libc::c_double = 0.0f64;
    let mut edgmax: libc::c_double = 0.;
    let mut pareta: libc::c_double = 0.;
    let mut prerec: libc::c_double = 0.0f64;
    let mut phimin: libc::c_double = 0.;
    let mut parsig: libc::c_double = 0.0f64;
    let mut gamma_: libc::c_double = 0.;
    let mut phi: libc::c_double = 0.;
    let mut rho: libc::c_double = 0.;
    let mut sum: libc::c_double = 0.0f64;
    let mut ratio: libc::c_double = 0.;
    let mut vmold: libc::c_double = 0.;
    let mut parmu: libc::c_double = 0.;
    let mut error: libc::c_double = 0.;
    let mut vmnew: libc::c_double = 0.;
    let mut resmax: libc::c_double = 0.;
    let mut cvmaxp: libc::c_double = 0.;
    let mut resnew: libc::c_double = 0.;
    let mut trured: libc::c_double = 0.;
    let mut temp: libc::c_double = 0.;
    let mut wsig: libc::c_double = 0.;
    let mut f: libc::c_double = 0.;
    let mut weta: libc::c_double = 0.;
    let mut i__: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut k: libc::c_int = 0;
    let mut l: libc::c_int = 0;
    let mut idxnew: libc::c_int = 0;
    let mut iflag: libc::c_int = 0 as libc::c_int;
    let mut iptemp: libc::c_int = 0;
    let mut isdirn: libc::c_int = 0;
    let mut izdota: libc::c_int = 0;
    let mut ivmc: libc::c_int = 0;
    let mut ivmd: libc::c_int = 0;
    let mut mp: libc::c_int = 0;
    let mut np: libc::c_int = 0;
    let mut iz: libc::c_int = 0;
    let mut ibrnch: libc::c_int = 0;
    let mut nbest: libc::c_int = 0;
    let mut ifull: libc::c_int = 0 as libc::c_int;
    let mut iptem: libc::c_int = 0;
    let mut jdrop: libc::c_int = 0;
    let mut rc: nlopt_result = NLOPT_SUCCESS;
    let mut seed: uint32_t = (*n + *m) as uint32_t;
    let mut feasible: libc::c_int = 0;
    *minf = ::std::f64::INFINITY;
    a_dim1 = *n;
    a_offset = 1 as libc::c_int + a_dim1 * 1 as libc::c_int;
    a = a.offset(-(a_offset as isize));
    simi_dim1 = *n;
    simi_offset = 1 as libc::c_int + simi_dim1 * 1 as libc::c_int;
    simi = simi.offset(-(simi_offset as isize));
    sim_dim1 = *n;
    sim_offset = 1 as libc::c_int + sim_dim1 * 1 as libc::c_int;
    sim = sim.offset(-(sim_offset as isize));
    datmat_dim1 = *mpp;
    datmat_offset = 1 as libc::c_int + datmat_dim1 * 1 as libc::c_int;
    datmat = datmat.offset(-(datmat_offset as isize));
    x = x.offset(-1);
    con = con.offset(-1);
    vsig = vsig.offset(-1);
    veta = veta.offset(-1);
    sigbar = sigbar.offset(-1);
    dx = dx.offset(-1);
    w = w.offset(-1);
    iact = iact.offset(-1);
    lb = lb.offset(-1);
    ub = ub.offset(-1);
    iptem = if *n <= 4 as libc::c_int {
        *n
    } else {
        4 as libc::c_int
    };
    iptemp = iptem + 1 as libc::c_int;
    np = *n + 1 as libc::c_int;
    mp = *m + 1 as libc::c_int;
    alpha = 0.25f64;
    beta = 2.1f64;
    gamma_ = 0.5f64;
    delta = 1.1f64;
    rho = *rhobeg;
    parmu = 0.0f64;
    if *iprint >= 2 as libc::c_int {
        fprintf(
            Io::stderr,
            &format!(
                "cobyla: the initial value of RHO is {} and PARMU is set to zero.",
                rho
            ),
        );
    }
    temp = 1.0f64 / rho;
    i__1 = *n;
    i__ = 1 as libc::c_int;
    while i__ <= i__1 {
        let mut rhocur: libc::c_double = 0.;
        *sim.offset((i__ + np * sim_dim1) as isize) = *x.offset(i__ as isize);
        i__2 = *n;
        j = 1 as libc::c_int;
        while j <= i__2 {
            *sim.offset((i__ + j * sim_dim1) as isize) = 0.0f64;
            *simi.offset((i__ + j * simi_dim1) as isize) = 0.0f64;
            j += 1;
        }
        rhocur = rho;
        if *x.offset(i__ as isize) + rhocur > *ub.offset(i__ as isize) {
            if *x.offset(i__ as isize) - rhocur >= *lb.offset(i__ as isize) {
                rhocur = -rhocur;
            } else if *ub.offset(i__ as isize) - *x.offset(i__ as isize)
                > *x.offset(i__ as isize) - *lb.offset(i__ as isize)
            {
                rhocur = 0.5f64 * (*ub.offset(i__ as isize) - *x.offset(i__ as isize));
            } else {
                rhocur = 0.5f64 * (*x.offset(i__ as isize) - *lb.offset(i__ as isize));
            }
        }
        *sim.offset((i__ + i__ * sim_dim1) as isize) = rhocur;
        *simi.offset((i__ + i__ * simi_dim1) as isize) = 1.0f64 / rhocur;
        i__ += 1;
    }
    jdrop = np;
    ibrnch = 0 as libc::c_int;
    'c_6104: loop {
        if nlopt_stop_forced(stop) != 0 {
            rc = NLOPT_FORCED_STOP;
        } else if *(*stop).nevals_p > 0 as libc::c_int {
            if nlopt_stop_evals(stop) != 0 {
                rc = NLOPT_MAXEVAL_REACHED;
            } else if nlopt_stop_time(stop) != 0 {
                rc = NLOPT_MAXTIME_REACHED;
            }
        }
        if rc as libc::c_int != NLOPT_SUCCESS as libc::c_int {
            current_block = 16949430136398296108;
            break;
        }
        let ref mut fresh3 = *(*stop).nevals_p;
        *fresh3 += 1;
        if calcfc.expect("non-null function pointer")(
            *n,
            *m,
            &mut *x.offset(1 as libc::c_int as isize),
            &mut f,
            &mut *con.offset(1 as libc::c_int as isize),
            state,
        ) != 0
        {
            if *iprint >= 1 as libc::c_int {
                fprintf(Io::stderr, "cobyla: user requested end of minimization");
            }
            rc = NLOPT_FORCED_STOP;
            current_block = 16949430136398296108;
            break;
        } else {
            resmax = 0.0f64;
            feasible = 1 as libc::c_int;
            if *m > 0 as libc::c_int {
                i__1 = *m;
                k = 1 as libc::c_int;
                while k <= i__1 {
                    d__1 = resmax;
                    d__2 = -*con.offset(k as isize);
                    resmax = if d__1 >= d__2 { d__1 } else { d__2 };
                    if d__2 > *((*state).con_tol).offset((k - 1 as libc::c_int) as isize) {
                        feasible = 0 as libc::c_int;
                    }
                    k += 1;
                }
            }
            if f < (*stop).minf_max && feasible != 0 {
                rc = NLOPT_STOPVAL_REACHED;
                current_block = 10710279849762345920;
                break;
            } else {
                if *(*stop).nevals_p == *iprint - 1 as libc::c_int || *iprint == 3 as libc::c_int {
                    fprintf(
                        Io::stderr,
                        &format!(
                            "cobyla: NFVALS = {}, F ={}, MAXCV ={}",
                            *(*stop).nevals_p,
                            f,
                            resmax
                        ),
                    );
                    i__1 = iptem;
                    fprintf(Io::stderr, "cobyla: X =");
                    i__ = 1 as libc::c_int;
                    while i__ <= i__1 {
                        if i__ > 1 as libc::c_int {
                            fprintf(Io::stderr, "  ");
                        }
                        fprintf(Io::stderr, &format!("{}", *x.offset(i__ as isize)));
                        i__ += 1;
                    }
                    if iptem < *n {
                        i__1 = *n;
                        i__ = iptemp;
                        while i__ <= i__1 {
                            if (i__ - 1 as libc::c_int) % 4 as libc::c_int == 0 {
                                fprintf(Io::stderr, "\ncobyla:  ");
                            }
                            fprintf(Io::stderr, &format!("{}", *x.offset(i__ as isize)));
                            i__ += 1;
                        }
                    }
                    fprintf(Io::stderr, "");
                }
                *con.offset(mp as isize) = f;
                *con.offset(*mpp as isize) = resmax;
                if ibrnch == 1 as libc::c_int {
                    vmold = *datmat.offset((mp + np * datmat_dim1) as isize)
                        + parmu * *datmat.offset((*mpp + np * datmat_dim1) as isize);
                    vmnew = f + parmu * resmax;
                    trured = vmold - vmnew;
                    if parmu == 0.0f64 && f == *datmat.offset((mp + np * datmat_dim1) as isize) {
                        prerem = prerec;
                        trured = *datmat.offset((*mpp + np * datmat_dim1) as isize) - resmax;
                    }
                    ratio = 0.0f64;
                    if trured <= 0.0f32 as libc::c_double {
                        ratio = 1.0f32 as libc::c_double;
                    }
                    jdrop = 0 as libc::c_int;
                    i__1 = *n;
                    j = 1 as libc::c_int;
                    while j <= i__1 {
                        temp = 0.0f64;
                        i__2 = *n;
                        i__ = 1 as libc::c_int;
                        while i__ <= i__2 {
                            temp += *simi.offset((j + i__ * simi_dim1) as isize)
                                * *dx.offset(i__ as isize);
                            i__ += 1;
                        }
                        temp = temp.abs();
                        if temp > ratio {
                            jdrop = j;
                            ratio = temp;
                        }
                        *sigbar.offset(j as isize) = temp * *vsig.offset(j as isize);
                        j += 1;
                    }
                    edgmax = delta * rho;
                    l = 0 as libc::c_int;
                    i__1 = *n;
                    j = 1 as libc::c_int;
                    while j <= i__1 {
                        if *sigbar.offset(j as isize) >= parsig
                            || *sigbar.offset(j as isize) >= *vsig.offset(j as isize)
                        {
                            temp = *veta.offset(j as isize);
                            if trured > 0.0f64 {
                                temp = 0.0f64;
                                i__2 = *n;
                                i__ = 1 as libc::c_int;
                                while i__ <= i__2 {
                                    d__1 = *dx.offset(i__ as isize)
                                        - *sim.offset((i__ + j * sim_dim1) as isize);
                                    temp += d__1 * d__1;
                                    i__ += 1;
                                }
                                temp = temp.sqrt();
                            }
                            if temp > edgmax {
                                l = j;
                                edgmax = temp;
                            }
                        }
                        j += 1;
                    }
                    if l > 0 as libc::c_int {
                        jdrop = l;
                    }
                    if jdrop == 0 as libc::c_int {
                        current_block = 17974563553836679504;
                    } else {
                        temp = 0.0f64;
                        i__1 = *n;
                        i__ = 1 as libc::c_int;
                        while i__ <= i__1 {
                            *sim.offset((i__ + jdrop * sim_dim1) as isize) =
                                *dx.offset(i__ as isize);
                            temp += *simi.offset((jdrop + i__ * simi_dim1) as isize)
                                * *dx.offset(i__ as isize);
                            i__ += 1;
                        }
                        i__1 = *n;
                        i__ = 1 as libc::c_int;
                        while i__ <= i__1 {
                            *simi.offset((jdrop + i__ * simi_dim1) as isize) /= temp;
                            i__ += 1;
                        }
                        i__1 = *n;
                        j = 1 as libc::c_int;
                        while j <= i__1 {
                            if j != jdrop {
                                temp = 0.0f64;
                                i__2 = *n;
                                i__ = 1 as libc::c_int;
                                while i__ <= i__2 {
                                    temp += *simi.offset((j + i__ * simi_dim1) as isize)
                                        * *dx.offset(i__ as isize);
                                    i__ += 1;
                                }
                                i__2 = *n;
                                i__ = 1 as libc::c_int;
                                while i__ <= i__2 {
                                    *simi.offset((j + i__ * simi_dim1) as isize) -=
                                        temp * *simi.offset((jdrop + i__ * simi_dim1) as isize);
                                    i__ += 1;
                                }
                            }
                            j += 1;
                        }
                        i__1 = *mpp;
                        k = 1 as libc::c_int;
                        while k <= i__1 {
                            *datmat.offset((k + jdrop * datmat_dim1) as isize) =
                                *con.offset(k as isize);
                            k += 1;
                        }
                        if trured > 0.0f64 && trured >= prerem * 0.1f64 {
                            if trured >= prerem * 0.9f64 && trured <= prerem * 1.1f64 && iflag != 0
                            {
                                rho *= 2.0f64;
                            }
                            current_block = 16207618807156029286;
                        } else {
                            current_block = 17974563553836679504;
                        }
                    }
                } else {
                    i__1 = *mpp;
                    k = 1 as libc::c_int;
                    while k <= i__1 {
                        *datmat.offset((k + jdrop * datmat_dim1) as isize) =
                            *con.offset(k as isize);
                        k += 1;
                    }
                    if !(*(*stop).nevals_p > np) {
                        if jdrop <= *n {
                            if *datmat.offset((mp + np * datmat_dim1) as isize) <= f {
                                *x.offset(jdrop as isize) =
                                    *sim.offset((jdrop + np * sim_dim1) as isize);
                            } else {
                                let mut rhocur_0: libc::c_double = *x.offset(jdrop as isize)
                                    - *sim.offset((jdrop + np * sim_dim1) as isize);
                                *sim.offset((jdrop + np * sim_dim1) as isize) =
                                    *x.offset(jdrop as isize);
                                i__1 = *mpp;
                                k = 1 as libc::c_int;
                                while k <= i__1 {
                                    *datmat.offset((k + jdrop * datmat_dim1) as isize) =
                                        *datmat.offset((k + np * datmat_dim1) as isize);
                                    *datmat.offset((k + np * datmat_dim1) as isize) =
                                        *con.offset(k as isize);
                                    k += 1;
                                }
                                i__1 = jdrop;
                                k = 1 as libc::c_int;
                                while k <= i__1 {
                                    *sim.offset((jdrop + k * sim_dim1) as isize) = -rhocur_0;
                                    temp = 0.0f32 as libc::c_double;
                                    i__2 = jdrop;
                                    i__ = k;
                                    while i__ <= i__2 {
                                        temp -= *simi.offset((i__ + k * simi_dim1) as isize);
                                        i__ += 1;
                                    }
                                    *simi.offset((jdrop + k * simi_dim1) as isize) = temp;
                                    k += 1;
                                }
                            }
                        }
                        if *(*stop).nevals_p <= *n {
                            jdrop = *(*stop).nevals_p;
                            *x.offset(jdrop as isize) +=
                                *sim.offset((jdrop + jdrop * sim_dim1) as isize);
                            continue;
                        }
                    }
                    ibrnch = 1 as libc::c_int;
                    current_block = 16207618807156029286;
                }
                'c_6122: loop {
                    match current_block {
                        17974563553836679504 => {
                            if iflag == 0 as libc::c_int {
                                ibrnch = 0 as libc::c_int;
                                current_block = 16207618807156029286;
                            } else {
                                let mut fbest: libc::c_double = if ifull == 1 as libc::c_int {
                                    f
                                } else {
                                    *datmat.offset((mp + np * datmat_dim1) as isize)
                                };
                                if fbest < *minf && nlopt_stop_ftol(stop, fbest, *minf) != 0 {
                                    rc = NLOPT_FTOL_REACHED;
                                    current_block = 16949430136398296108;
                                    break 'c_6104;
                                } else {
                                    *minf = fbest;
                                    if rho > rhoend {
                                        rho *= 0.5f64;
                                        if rho <= rhoend * 1.5f64 {
                                            rho = rhoend;
                                        }
                                        if parmu > 0.0f64 {
                                            denom = 0.0f64;
                                            i__1 = mp;
                                            k = 1 as libc::c_int;
                                            while k <= i__1 {
                                                cmin =
                                                    *datmat.offset((k + np * datmat_dim1) as isize);
                                                cmax = cmin;
                                                i__2 = *n;
                                                i__ = 1 as libc::c_int;
                                                while i__ <= i__2 {
                                                    d__1 = cmin;
                                                    d__2 = *datmat
                                                        .offset((k + i__ * datmat_dim1) as isize);
                                                    cmin = if d__1 <= d__2 { d__1 } else { d__2 };
                                                    d__1 = cmax;
                                                    d__2 = *datmat
                                                        .offset((k + i__ * datmat_dim1) as isize);
                                                    cmax = if d__1 >= d__2 { d__1 } else { d__2 };
                                                    i__ += 1;
                                                }
                                                if k <= *m && cmin < cmax * 0.5f64 {
                                                    temp = (if cmax >= 0.0f64 {
                                                        cmax
                                                    } else {
                                                        0.0f64
                                                    }) - cmin;
                                                    if denom <= 0.0f64 {
                                                        denom = temp;
                                                    } else {
                                                        denom = if denom <= temp {
                                                            denom
                                                        } else {
                                                            temp
                                                        };
                                                    }
                                                }
                                                k += 1;
                                            }
                                            if denom == 0.0f64 {
                                                parmu = 0.0f64;
                                            } else if cmax - cmin < parmu * denom {
                                                parmu = (cmax - cmin) / denom;
                                            }
                                        }
                                        if *iprint >= 2 as libc::c_int {
                                            fprintf(
                                                Io::stderr,
                                                &format!(
                                                    "cobyla: reduction in RHO to {} and PARMU ={}",
                                                    rho, parmu
                                                ),
                                            );
                                        }
                                        if *iprint == 2 as libc::c_int {
                                            fprintf(
                                                Io::stderr,
                                                &format!(
                                                    "cobyla: NFVALS = {}, F ={}, MAXCV ={}",
                                                    *(*stop).nevals_p,
                                                    *datmat
                                                        .offset((mp + np * datmat_dim1) as isize),
                                                    *datmat
                                                        .offset((*mpp + np * datmat_dim1) as isize)
                                                ),
                                            );
                                            fprintf(Io::stderr, "cobyla: X =");
                                            i__1 = iptem;
                                            i__ = 1 as libc::c_int;
                                            while i__ <= i__1 {
                                                if i__ > 1 as libc::c_int {
                                                    fprintf(Io::stderr, "  ");
                                                }
                                                fprintf(
                                                    Io::stderr,
                                                    &format!(
                                                        "{}",
                                                        *sim.offset((i__ + np * sim_dim1) as isize)
                                                    ),
                                                );
                                                i__ += 1;
                                            }
                                            if iptem < *n {
                                                i__1 = *n;
                                                i__ = iptemp;
                                                while i__ <= i__1 {
                                                    if (i__ - 1 as libc::c_int) % 4 as libc::c_int
                                                        == 0
                                                    {
                                                        fprintf(Io::stderr, "\ncobyla:  ");
                                                    }
                                                    fprintf(
                                                        Io::stderr,
                                                        &format!("{}", *x.offset(i__ as isize)),
                                                    );
                                                    i__ += 1;
                                                }
                                            }
                                            fprintf(Io::stderr, "");
                                        }
                                        current_block = 16207618807156029286;
                                    } else {
                                        rc = (if rhoend > 0 as libc::c_int as libc::c_double {
                                            NLOPT_XTOL_REACHED as libc::c_int
                                        } else {
                                            NLOPT_ROUNDOFF_LIMITED as libc::c_int
                                        })
                                            as nlopt_result;
                                        if *iprint >= 1 as libc::c_int {
                                            fprintf(Io::stderr, "cobyla: normal return.");
                                        }
                                        if ifull == 1 as libc::c_int {
                                            current_block = 10710279849762345920;
                                            break 'c_6104;
                                        } else {
                                            current_block = 16949430136398296108;
                                            break 'c_6104;
                                        }
                                    }
                                }
                            }
                        }
                        _ => {
                            phimin = *datmat.offset((mp + np * datmat_dim1) as isize)
                                + parmu * *datmat.offset((*mpp + np * datmat_dim1) as isize);
                            nbest = np;
                            i__1 = *n;
                            j = 1 as libc::c_int;
                            while j <= i__1 {
                                temp = *datmat.offset((mp + j * datmat_dim1) as isize)
                                    + parmu * *datmat.offset((*mpp + j * datmat_dim1) as isize);
                                if temp < phimin {
                                    nbest = j;
                                    phimin = temp;
                                } else if temp == phimin && parmu == 0.0f64 {
                                    if *datmat.offset((*mpp + j * datmat_dim1) as isize)
                                        < *datmat.offset((*mpp + nbest * datmat_dim1) as isize)
                                    {
                                        nbest = j;
                                    }
                                }
                                j += 1;
                            }
                            if nbest <= *n {
                                i__1 = *mpp;
                                i__ = 1 as libc::c_int;
                                while i__ <= i__1 {
                                    temp = *datmat.offset((i__ + np * datmat_dim1) as isize);
                                    *datmat.offset((i__ + np * datmat_dim1) as isize) =
                                        *datmat.offset((i__ + nbest * datmat_dim1) as isize);
                                    *datmat.offset((i__ + nbest * datmat_dim1) as isize) = temp;
                                    i__ += 1;
                                }
                                i__1 = *n;
                                i__ = 1 as libc::c_int;
                                while i__ <= i__1 {
                                    temp = *sim.offset((i__ + nbest * sim_dim1) as isize);
                                    *sim.offset((i__ + nbest * sim_dim1) as isize) = 0.0f64;
                                    *sim.offset((i__ + np * sim_dim1) as isize) += temp;
                                    tempa = 0.0f64;
                                    i__2 = *n;
                                    k = 1 as libc::c_int;
                                    while k <= i__2 {
                                        *sim.offset((i__ + k * sim_dim1) as isize) -= temp;
                                        tempa -= *simi.offset((k + i__ * simi_dim1) as isize);
                                        k += 1;
                                    }
                                    *simi.offset((nbest + i__ * simi_dim1) as isize) = tempa;
                                    i__ += 1;
                                }
                            }
                            error = 0.0f64;
                            i__1 = *n;
                            i__ = 1 as libc::c_int;
                            while i__ <= i__1 {
                                i__2 = *n;
                                j = 1 as libc::c_int;
                                while j <= i__2 {
                                    temp = 0.0f64;
                                    if i__ == j {
                                        temp += -1.0f64;
                                    }
                                    i__3 = *n;
                                    k = 1 as libc::c_int;
                                    while k <= i__3 {
                                        if *sim.offset((k + j * sim_dim1) as isize)
                                            != 0 as libc::c_int as libc::c_double
                                        {
                                            temp += *simi.offset((i__ + k * simi_dim1) as isize)
                                                * *sim.offset((k + j * sim_dim1) as isize);
                                        }
                                        k += 1;
                                    }
                                    d__1 = error;
                                    d__2 = temp.abs();
                                    error = if d__1 >= d__2 { d__1 } else { d__2 };
                                    j += 1;
                                }
                                i__ += 1;
                            }
                            if error > 0.1f64 {
                                if *iprint >= 1 as libc::c_int {
                                    fprintf(
                                        Io::stderr,
                                        "cobyla: rounding errors are becoming damaging.",
                                    );
                                }
                                rc = NLOPT_ROUNDOFF_LIMITED;
                                current_block = 16949430136398296108;
                                break 'c_6104;
                            } else {
                                i__2 = mp;
                                k = 1 as libc::c_int;
                                while k <= i__2 {
                                    *con.offset(k as isize) =
                                        -*datmat.offset((k + np * datmat_dim1) as isize);
                                    i__1 = *n;
                                    j = 1 as libc::c_int;
                                    while j <= i__1 {
                                        *w.offset(j as isize) = *datmat
                                            .offset((k + j * datmat_dim1) as isize)
                                            + *con.offset(k as isize);
                                        j += 1;
                                    }
                                    i__1 = *n;
                                    i__ = 1 as libc::c_int;
                                    while i__ <= i__1 {
                                        temp = 0.0f64;
                                        i__3 = *n;
                                        j = 1 as libc::c_int;
                                        while j <= i__3 {
                                            temp += *w.offset(j as isize)
                                                * *simi.offset((j + i__ * simi_dim1) as isize);
                                            j += 1;
                                        }
                                        if k == mp {
                                            temp = -temp;
                                        }
                                        *a.offset((i__ + k * a_dim1) as isize) = temp;
                                        i__ += 1;
                                    }
                                    k += 1;
                                }
                                iflag = 1 as libc::c_int;
                                parsig = alpha * rho;
                                pareta = beta * rho;
                                i__1 = *n;
                                j = 1 as libc::c_int;
                                while j <= i__1 {
                                    wsig = 0.0f64;
                                    weta = 0.0f64;
                                    i__2 = *n;
                                    i__ = 1 as libc::c_int;
                                    while i__ <= i__2 {
                                        d__1 = *simi.offset((j + i__ * simi_dim1) as isize);
                                        wsig += d__1 * d__1;
                                        d__1 = *sim.offset((i__ + j * sim_dim1) as isize);
                                        weta += d__1 * d__1;
                                        i__ += 1;
                                    }
                                    *vsig.offset(j as isize) = 1.0f64 / wsig.sqrt();
                                    *veta.offset(j as isize) = weta.sqrt();
                                    if *vsig.offset(j as isize) < parsig
                                        || *veta.offset(j as isize) > pareta
                                    {
                                        iflag = 0 as libc::c_int;
                                    }
                                    j += 1;
                                }
                                if ibrnch == 1 as libc::c_int || iflag == 1 as libc::c_int {
                                    iz = 1 as libc::c_int;
                                    izdota = iz + *n * *n;
                                    ivmc = izdota + *n;
                                    isdirn = ivmc + mp;
                                    idxnew = isdirn + *n;
                                    ivmd = idxnew + *n;
                                    rc = trstlp(
                                        n,
                                        m,
                                        &mut *a.offset(a_offset as isize),
                                        &mut *con.offset(1 as libc::c_int as isize),
                                        &mut rho,
                                        &mut *dx.offset(1 as libc::c_int as isize),
                                        &mut ifull,
                                        &mut *iact.offset(1 as libc::c_int as isize),
                                        &mut *w.offset(iz as isize),
                                        &mut *w.offset(izdota as isize),
                                        &mut *w.offset(ivmc as isize),
                                        &mut *w.offset(isdirn as isize),
                                        &mut *w.offset(idxnew as isize),
                                        &mut *w.offset(ivmd as isize),
                                    );
                                    if rc as libc::c_int != NLOPT_SUCCESS as libc::c_int {
                                        current_block = 16949430136398296108;
                                        break 'c_6104;
                                    }
                                    i__1 = *n;
                                    i__ = 1 as libc::c_int;
                                    while i__ <= i__1 {
                                        let mut xi_0: libc::c_double =
                                            *sim.offset((i__ + np * sim_dim1) as isize);
                                        if xi_0 + *dx.offset(i__ as isize)
                                            > *ub.offset(i__ as isize)
                                        {
                                            *dx.offset(i__ as isize) =
                                                *ub.offset(i__ as isize) - xi_0;
                                        }
                                        if xi_0 + *dx.offset(i__ as isize)
                                            < *lb.offset(i__ as isize)
                                        {
                                            *dx.offset(i__ as isize) =
                                                xi_0 - *lb.offset(i__ as isize);
                                        }
                                        i__ += 1;
                                    }
                                    if ifull == 0 as libc::c_int {
                                        temp = 0.0f64;
                                        i__1 = *n;
                                        i__ = 1 as libc::c_int;
                                        while i__ <= i__1 {
                                            d__1 = *dx.offset(i__ as isize);
                                            temp += d__1 * d__1;
                                            i__ += 1;
                                        }
                                        if temp < rho * 0.25f64 * rho {
                                            ibrnch = 1 as libc::c_int;
                                            current_block = 17974563553836679504;
                                            continue;
                                        }
                                    }
                                    resnew = 0.0f64;
                                    *con.offset(mp as isize) = 0.0f64;
                                    i__1 = mp;
                                    k = 1 as libc::c_int;
                                    while k <= i__1 {
                                        sum = *con.offset(k as isize);
                                        i__2 = *n;
                                        i__ = 1 as libc::c_int;
                                        while i__ <= i__2 {
                                            sum -= *a.offset((i__ + k * a_dim1) as isize)
                                                * *dx.offset(i__ as isize);
                                            i__ += 1;
                                        }
                                        if k < mp {
                                            resnew = if resnew >= sum { resnew } else { sum };
                                        }
                                        k += 1;
                                    }
                                    barmu = 0.0f64;
                                    prerec =
                                        *datmat.offset((*mpp + np * datmat_dim1) as isize) - resnew;
                                    if prerec > 0.0f64 {
                                        barmu = sum / prerec;
                                    }
                                    if !(parmu < barmu * 1.5f64) {
                                        break;
                                    }
                                    parmu = barmu * 2.0f64;
                                    if *iprint >= 2 as libc::c_int {
                                        fprintf(
                                            Io::stderr,
                                            &format!("cobyla: increase in PARMU to {}", parmu),
                                        );
                                    }
                                    phi = *datmat.offset((mp + np * datmat_dim1) as isize)
                                        + parmu
                                            * *datmat.offset((*mpp + np * datmat_dim1) as isize);
                                    i__1 = *n;
                                    j = 1 as libc::c_int;
                                    loop {
                                        if !(j <= i__1) {
                                            break 'c_6122;
                                        }
                                        temp = *datmat.offset((mp + j * datmat_dim1) as isize)
                                            + parmu
                                                * *datmat.offset((*mpp + j * datmat_dim1) as isize);
                                        if temp < phi {
                                            current_block = 16207618807156029286;
                                            break;
                                        }
                                        if temp == phi && parmu == 0.0f32 as libc::c_double {
                                            if *datmat.offset((*mpp + j * datmat_dim1) as isize)
                                                < *datmat.offset((*mpp + np * datmat_dim1) as isize)
                                            {
                                                current_block = 16207618807156029286;
                                                break;
                                            }
                                        }
                                        j += 1;
                                    }
                                } else {
                                    jdrop = 0 as libc::c_int;
                                    temp = pareta;
                                    i__1 = *n;
                                    j = 1 as libc::c_int;
                                    while j <= i__1 {
                                        if *veta.offset(j as isize) > temp {
                                            jdrop = j;
                                            temp = *veta.offset(j as isize);
                                        }
                                        j += 1;
                                    }
                                    if jdrop == 0 as libc::c_int {
                                        i__1 = *n;
                                        j = 1 as libc::c_int;
                                        while j <= i__1 {
                                            if *vsig.offset(j as isize) < temp {
                                                jdrop = j;
                                                temp = *vsig.offset(j as isize);
                                            }
                                            j += 1;
                                        }
                                    }
                                    temp = gamma_ * rho * *vsig.offset(jdrop as isize);
                                    i__1 = *n;
                                    i__ = 1 as libc::c_int;
                                    while i__ <= i__1 {
                                        *dx.offset(i__ as isize) =
                                            temp * *simi.offset((jdrop + i__ * simi_dim1) as isize);
                                        i__ += 1;
                                    }
                                    cvmaxp = 0.0f64;
                                    cvmaxm = 0.0f64;
                                    i__1 = mp;
                                    k = 1 as libc::c_int;
                                    while k <= i__1 {
                                        sum = 0.0f64;
                                        i__2 = *n;
                                        i__ = 1 as libc::c_int;
                                        while i__ <= i__2 {
                                            sum += *a.offset((i__ + k * a_dim1) as isize)
                                                * *dx.offset(i__ as isize);
                                            i__ += 1;
                                        }
                                        if k < mp {
                                            temp = *datmat.offset((k + np * datmat_dim1) as isize);
                                            d__1 = cvmaxp;
                                            d__2 = -sum - temp;
                                            cvmaxp = if d__1 >= d__2 { d__1 } else { d__2 };
                                            d__1 = cvmaxm;
                                            d__2 = sum - temp;
                                            cvmaxm = if d__1 >= d__2 { d__1 } else { d__2 };
                                        }
                                        k += 1;
                                    }
                                    dxsign = 1.0f64;
                                    if parmu * (cvmaxp - cvmaxm) > sum + sum {
                                        dxsign = -1.0f64;
                                    }
                                    temp = 0.0f64;
                                    i__1 = *n;
                                    i__ = 1 as libc::c_int;
                                    while i__ <= i__1 {
                                        *dx.offset(i__ as isize) = dxsign
                                            * *dx.offset(i__ as isize)
                                            * lcg_urand(
                                                &mut seed,
                                                0.01f64,
                                                1 as libc::c_int as libc::c_double,
                                            );
                                        let mut xi: libc::c_double =
                                            *sim.offset((i__ + np * sim_dim1) as isize);
                                        loop {
                                            if xi + *dx.offset(i__ as isize)
                                                > *ub.offset(i__ as isize)
                                            {
                                                *dx.offset(i__ as isize) =
                                                    -*dx.offset(i__ as isize);
                                            }
                                            if !(xi + *dx.offset(i__ as isize)
                                                < *lb.offset(i__ as isize))
                                            {
                                                break;
                                            }
                                            if xi - *dx.offset(i__ as isize)
                                                <= *ub.offset(i__ as isize)
                                            {
                                                *dx.offset(i__ as isize) =
                                                    -*dx.offset(i__ as isize);
                                                break;
                                            } else {
                                                *dx.offset(i__ as isize) *= 0.5f64;
                                            }
                                        }
                                        *sim.offset((i__ + jdrop * sim_dim1) as isize) =
                                            *dx.offset(i__ as isize);
                                        temp += *simi.offset((jdrop + i__ * simi_dim1) as isize)
                                            * *dx.offset(i__ as isize);
                                        i__ += 1;
                                    }
                                    i__1 = *n;
                                    i__ = 1 as libc::c_int;
                                    while i__ <= i__1 {
                                        *simi.offset((jdrop + i__ * simi_dim1) as isize) /= temp;
                                        i__ += 1;
                                    }
                                    i__1 = *n;
                                    j = 1 as libc::c_int;
                                    while j <= i__1 {
                                        if j != jdrop {
                                            temp = 0.0f64;
                                            i__2 = *n;
                                            i__ = 1 as libc::c_int;
                                            while i__ <= i__2 {
                                                temp += *simi
                                                    .offset((j + i__ * simi_dim1) as isize)
                                                    * *dx.offset(i__ as isize);
                                                i__ += 1;
                                            }
                                            i__2 = *n;
                                            i__ = 1 as libc::c_int;
                                            while i__ <= i__2 {
                                                *simi.offset((j + i__ * simi_dim1) as isize) -= temp
                                                    * *simi
                                                        .offset((jdrop + i__ * simi_dim1) as isize);
                                                i__ += 1;
                                            }
                                        }
                                        *x.offset(j as isize) = *sim
                                            .offset((j + np * sim_dim1) as isize)
                                            + *dx.offset(j as isize);
                                        j += 1;
                                    }
                                    continue 'c_6104;
                                }
                            }
                        }
                    }
                }
                prerem = parmu * prerec - sum;
                i__1 = *n;
                i__ = 1 as libc::c_int;
                while i__ <= i__1 {
                    *x.offset(i__ as isize) =
                        *sim.offset((i__ + np * sim_dim1) as isize) + *dx.offset(i__ as isize);
                    i__ += 1;
                }
                ibrnch = 1 as libc::c_int;
            }
        }
    }
    match current_block {
        16949430136398296108 => {
            i__1 = *n;
            i__ = 1 as libc::c_int;
            while i__ <= i__1 {
                *x.offset(i__ as isize) = *sim.offset((i__ + np * sim_dim1) as isize);
                i__ += 1;
            }
            f = *datmat.offset((mp + np * datmat_dim1) as isize);
            resmax = *datmat.offset((*mpp + np * datmat_dim1) as isize);
        }
        _ => {}
    }
    *minf = f;
    if *iprint >= 1 as libc::c_int {
        fprintf(
            Io::stderr,
            &format!(
                "cobyla: NFVALS = {}, F ={}, MAXCV ={}\n",
                *(*stop).nevals_p,
                f,
                resmax
            ),
        );
        i__1 = iptem;
        fprintf(Io::stderr, "cobyla: X =");
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            if i__ > 1 as libc::c_int {
                fprintf(Io::stderr, "  ");
            }
            fprintf(Io::stderr, &format!("{}", *x.offset(i__ as isize)));
            i__ += 1;
        }
        if iptem < *n {
            i__1 = *n;
            i__ = iptemp;
            while i__ <= i__1 {
                if (i__ - 1 as libc::c_int) % 4 as libc::c_int == 0 {
                    fprintf(Io::stderr, "\ncobyla:  ");
                }
                fprintf(Io::stderr, &format!("{}", *x.offset(i__ as isize)));
                i__ += 1;
            }
        }
        fprintf(Io::stderr, "");
    }
    rc
}
unsafe fn trstlp(
    mut n: *mut libc::c_int,
    mut m: *mut libc::c_int,
    mut a: *mut libc::c_double,
    mut b: *mut libc::c_double,
    mut rho: *mut libc::c_double,
    mut dx: *mut libc::c_double,
    mut ifull: *mut libc::c_int,
    mut iact: *mut libc::c_int,
    mut z__: *mut libc::c_double,
    mut zdota: *mut libc::c_double,
    mut vmultc: *mut libc::c_double,
    mut sdirn: *mut libc::c_double,
    mut dxnew: *mut libc::c_double,
    mut vmultd: *mut libc::c_double,
) -> nlopt_result {
    let mut current_block: u64;
    let mut a_dim1: libc::c_int = 0;
    let mut a_offset: libc::c_int = 0;
    let mut z_dim1: libc::c_int = 0;
    let mut z_offset: libc::c_int = 0;
    let mut i__1: libc::c_int = 0;
    let mut i__2: libc::c_int = 0;
    let mut d__1: libc::c_double = 0.;
    let mut d__2: libc::c_double = 0.;
    let mut alpha: libc::c_double = 0.;
    let mut tempa: libc::c_double = 0.;
    let mut beta: libc::c_double = 0.;
    let mut optnew: libc::c_double = 0.;
    let mut stpful: libc::c_double = 0.;
    let mut sum: libc::c_double = 0.;
    let mut tot: libc::c_double = 0.;
    let mut acca: libc::c_double = 0.;
    let mut accb: libc::c_double = 0.;
    let mut ratio: libc::c_double = 0.;
    let mut vsave: libc::c_double = 0.;
    let mut zdotv: libc::c_double = 0.;
    let mut zdotw: libc::c_double = 0.;
    let mut dd: libc::c_double = 0.;
    let mut sd: libc::c_double = 0.;
    let mut sp: libc::c_double = 0.;
    let mut ss: libc::c_double = 0.;
    let mut resold: libc::c_double = 0.0f64;
    let mut zdvabs: libc::c_double = 0.;
    let mut zdwabs: libc::c_double = 0.;
    let mut sumabs: libc::c_double = 0.;
    let mut resmax: libc::c_double = 0.;
    let mut optold: libc::c_double = 0.;
    let mut spabs: libc::c_double = 0.;
    let mut temp: libc::c_double = 0.;
    let mut step: libc::c_double = 0.;
    let mut icount: libc::c_int = 0;
    let mut i__: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut k: libc::c_int = 0;
    let mut isave: libc::c_int = 0;
    let mut kk: libc::c_int = 0;
    let mut kl: libc::c_int = 0;
    let mut kp: libc::c_int = 0;
    let mut kw: libc::c_int = 0;
    let mut nact: libc::c_int = 0;
    let mut icon: libc::c_int = 0 as libc::c_int;
    let mut mcon: libc::c_int = 0;
    let mut nactx: libc::c_int = 0 as libc::c_int;
    z_dim1 = *n;
    z_offset = 1 as libc::c_int + z_dim1 * 1 as libc::c_int;
    z__ = z__.offset(-(z_offset as isize));
    a_dim1 = *n;
    a_offset = 1 as libc::c_int + a_dim1 * 1 as libc::c_int;
    a = a.offset(-(a_offset as isize));
    b = b.offset(-1);
    dx = dx.offset(-1);
    iact = iact.offset(-1);
    zdota = zdota.offset(-1);
    vmultc = vmultc.offset(-1);
    sdirn = sdirn.offset(-1);
    dxnew = dxnew.offset(-1);
    vmultd = vmultd.offset(-1);
    *ifull = 1 as libc::c_int;
    mcon = *m;
    nact = 0 as libc::c_int;
    resmax = 0.0f64;
    i__1 = *n;
    i__ = 1 as libc::c_int;
    while i__ <= i__1 {
        i__2 = *n;
        j = 1 as libc::c_int;
        while j <= i__2 {
            *z__.offset((i__ + j * z_dim1) as isize) = 0.0f64;
            j += 1;
        }
        *z__.offset((i__ + i__ * z_dim1) as isize) = 1.0f64;
        *dx.offset(i__ as isize) = 0.0f64;
        i__ += 1;
    }
    if *m >= 1 as libc::c_int {
        i__1 = *m;
        k = 1 as libc::c_int;
        while k <= i__1 {
            if *b.offset(k as isize) > resmax {
                resmax = *b.offset(k as isize);
                icon = k;
            }
            k += 1;
        }
        i__1 = *m;
        k = 1 as libc::c_int;
        while k <= i__1 {
            *iact.offset(k as isize) = k;
            *vmultc.offset(k as isize) = resmax - *b.offset(k as isize);
            k += 1;
        }
    }
    if resmax == 0.0f64 {
        current_block = 11188143500741601598;
    } else {
        i__1 = *n;
        i__ = 1 as libc::c_int;
        while i__ <= i__1 {
            *sdirn.offset(i__ as isize) = 0.0f64;
            i__ += 1;
        }
        current_block = 13859042411183768487;
    }
    'c_8601: loop {
        match current_block {
            11188143500741601598 => {
                mcon = *m + 1 as libc::c_int;
                icon = mcon;
                *iact.offset(mcon as isize) = mcon;
                *vmultc.offset(mcon as isize) = 0.0f64;
                current_block = 13859042411183768487;
            }
            _ => {
                optold = 0.0f64;
                icount = 0 as libc::c_int;
                loop {
                    if mcon == *m {
                        optnew = resmax;
                    } else {
                        optnew = 0.0f64;
                        i__1 = *n;
                        i__ = 1 as libc::c_int;
                        while i__ <= i__1 {
                            optnew -= *dx.offset(i__ as isize)
                                * *a.offset((i__ + mcon * a_dim1) as isize);
                            i__ += 1;
                        }
                    }
                    if icount == 0 as libc::c_int || optnew < optold {
                        optold = optnew;
                        nactx = nact;
                        icount = 3 as libc::c_int;
                    } else if nact > nactx {
                        nactx = nact;
                        icount = 3 as libc::c_int;
                    } else {
                        icount -= 1;
                        if icount == 0 as libc::c_int {
                            break;
                        }
                    }
                    if icon <= nact {
                        if icon < nact {
                            isave = *iact.offset(icon as isize);
                            vsave = *vmultc.offset(icon as isize);
                            k = icon;
                            loop {
                                kp = k + 1 as libc::c_int;
                                kk = *iact.offset(kp as isize);
                                sp = 0.0f64;
                                i__1 = *n;
                                i__ = 1 as libc::c_int;
                                while i__ <= i__1 {
                                    sp += *z__.offset((i__ + k * z_dim1) as isize)
                                        * *a.offset((i__ + kk * a_dim1) as isize);
                                    i__ += 1;
                                }
                                d__1 = *zdota.offset(kp as isize);
                                temp = (sp * sp + d__1 * d__1).sqrt();
                                alpha = *zdota.offset(kp as isize) / temp;
                                beta = sp / temp;
                                *zdota.offset(kp as isize) = alpha * *zdota.offset(k as isize);
                                *zdota.offset(k as isize) = temp;
                                i__1 = *n;
                                i__ = 1 as libc::c_int;
                                while i__ <= i__1 {
                                    temp = alpha * *z__.offset((i__ + kp * z_dim1) as isize)
                                        + beta * *z__.offset((i__ + k * z_dim1) as isize);
                                    *z__.offset((i__ + kp * z_dim1) as isize) = alpha
                                        * *z__.offset((i__ + k * z_dim1) as isize)
                                        - beta * *z__.offset((i__ + kp * z_dim1) as isize);
                                    *z__.offset((i__ + k * z_dim1) as isize) = temp;
                                    i__ += 1;
                                }
                                *iact.offset(k as isize) = kk;
                                *vmultc.offset(k as isize) = *vmultc.offset(kp as isize);
                                k = kp;
                                if !(k < nact) {
                                    break;
                                }
                            }
                            *iact.offset(k as isize) = isave;
                            *vmultc.offset(k as isize) = vsave;
                        }
                        nact -= 1;
                        if mcon > *m {
                            current_block = 15623375721314334080;
                        } else {
                            temp = 0.0f64;
                            i__1 = *n;
                            i__ = 1 as libc::c_int;
                            while i__ <= i__1 {
                                temp += *sdirn.offset(i__ as isize)
                                    * *z__.offset(
                                        (i__ + (nact + 1 as libc::c_int) * z_dim1) as isize,
                                    );
                                i__ += 1;
                            }
                            i__1 = *n;
                            i__ = 1 as libc::c_int;
                            while i__ <= i__1 {
                                *sdirn.offset(i__ as isize) -= temp
                                    * *z__.offset(
                                        (i__ + (nact + 1 as libc::c_int) * z_dim1) as isize,
                                    );
                                i__ += 1;
                            }
                            current_block = 9824026811267781195;
                        }
                    } else {
                        kk = *iact.offset(icon as isize);
                        i__1 = *n;
                        i__ = 1 as libc::c_int;
                        while i__ <= i__1 {
                            *dxnew.offset(i__ as isize) = *a.offset((i__ + kk * a_dim1) as isize);
                            i__ += 1;
                        }
                        tot = 0.0f64;
                        k = *n;
                        while k > nact {
                            sp = 0.0f64;
                            spabs = 0.0f64;
                            i__1 = *n;
                            i__ = 1 as libc::c_int;
                            while i__ <= i__1 {
                                temp = *z__.offset((i__ + k * z_dim1) as isize)
                                    * *dxnew.offset(i__ as isize);
                                sp += temp;
                                spabs += (temp).abs();
                                i__ += 1;
                            }
                            acca = spabs + (sp).abs() * 0.1f64;
                            accb = spabs + (sp).abs() * 0.2f64;
                            if spabs >= acca || acca >= accb {
                                sp = 0.0f64;
                            }
                            if tot == 0.0f64 {
                                tot = sp;
                            } else {
                                kp = k + 1 as libc::c_int;
                                temp = (sp * sp + tot * tot).sqrt();
                                alpha = sp / temp;
                                beta = tot / temp;
                                tot = temp;
                                i__1 = *n;
                                i__ = 1 as libc::c_int;
                                while i__ <= i__1 {
                                    temp = alpha * *z__.offset((i__ + k * z_dim1) as isize)
                                        + beta * *z__.offset((i__ + kp * z_dim1) as isize);
                                    *z__.offset((i__ + kp * z_dim1) as isize) = alpha
                                        * *z__.offset((i__ + kp * z_dim1) as isize)
                                        - beta * *z__.offset((i__ + k * z_dim1) as isize);
                                    *z__.offset((i__ + k * z_dim1) as isize) = temp;
                                    i__ += 1;
                                }
                            }
                            k -= 1;
                        }
                        if tot != 0.0f64 {
                            nact += 1;
                            *zdota.offset(nact as isize) = tot;
                            *vmultc.offset(icon as isize) = *vmultc.offset(nact as isize);
                            *vmultc.offset(nact as isize) = 0.0f64;
                        } else {
                            ratio = -1.0f64;
                            k = nact;
                            loop {
                                zdotv = 0.0f64;
                                zdvabs = 0.0f64;
                                i__1 = *n;
                                i__ = 1 as libc::c_int;
                                while i__ <= i__1 {
                                    temp = *z__.offset((i__ + k * z_dim1) as isize)
                                        * *dxnew.offset(i__ as isize);
                                    zdotv += temp;
                                    zdvabs += (temp).abs();
                                    i__ += 1;
                                }
                                acca = zdvabs + (zdotv).abs() * 0.1f64;
                                accb = zdvabs + (zdotv).abs() * 0.2f64;
                                if zdvabs < acca && acca < accb {
                                    temp = zdotv / *zdota.offset(k as isize);
                                    if temp > 0.0f64 && *iact.offset(k as isize) <= *m {
                                        tempa = *vmultc.offset(k as isize) / temp;
                                        if ratio < 0.0f64 || tempa < ratio {
                                            ratio = tempa;
                                        }
                                    }
                                    if k >= 2 as libc::c_int {
                                        kw = *iact.offset(k as isize);
                                        i__1 = *n;
                                        i__ = 1 as libc::c_int;
                                        while i__ <= i__1 {
                                            *dxnew.offset(i__ as isize) -=
                                                temp * *a.offset((i__ + kw * a_dim1) as isize);
                                            i__ += 1;
                                        }
                                    }
                                    *vmultd.offset(k as isize) = temp;
                                } else {
                                    *vmultd.offset(k as isize) = 0.0f64;
                                }
                                k -= 1;
                                if !(k > 0 as libc::c_int) {
                                    break;
                                }
                            }
                            if ratio < 0.0f64 {
                                break;
                            }
                            i__1 = nact;
                            k = 1 as libc::c_int;
                            while k <= i__1 {
                                d__1 = 0.0f64;
                                d__2 =
                                    *vmultc.offset(k as isize) - ratio * *vmultd.offset(k as isize);
                                *vmultc.offset(k as isize) = if d__1 >= d__2 { d__1 } else { d__2 };
                                k += 1;
                            }
                            if icon < nact {
                                isave = *iact.offset(icon as isize);
                                vsave = *vmultc.offset(icon as isize);
                                k = icon;
                                loop {
                                    kp = k + 1 as libc::c_int;
                                    kw = *iact.offset(kp as isize);
                                    sp = 0.0f64;
                                    i__1 = *n;
                                    i__ = 1 as libc::c_int;
                                    while i__ <= i__1 {
                                        sp += *z__.offset((i__ + k * z_dim1) as isize)
                                            * *a.offset((i__ + kw * a_dim1) as isize);
                                        i__ += 1;
                                    }
                                    d__1 = *zdota.offset(kp as isize);
                                    temp = (sp * sp + d__1 * d__1).sqrt();
                                    alpha = *zdota.offset(kp as isize) / temp;
                                    beta = sp / temp;
                                    *zdota.offset(kp as isize) = alpha * *zdota.offset(k as isize);
                                    *zdota.offset(k as isize) = temp;
                                    i__1 = *n;
                                    i__ = 1 as libc::c_int;
                                    while i__ <= i__1 {
                                        temp = alpha * *z__.offset((i__ + kp * z_dim1) as isize)
                                            + beta * *z__.offset((i__ + k * z_dim1) as isize);
                                        *z__.offset((i__ + kp * z_dim1) as isize) = alpha
                                            * *z__.offset((i__ + k * z_dim1) as isize)
                                            - beta * *z__.offset((i__ + kp * z_dim1) as isize);
                                        *z__.offset((i__ + k * z_dim1) as isize) = temp;
                                        i__ += 1;
                                    }
                                    *iact.offset(k as isize) = kw;
                                    *vmultc.offset(k as isize) = *vmultc.offset(kp as isize);
                                    k = kp;
                                    if !(k < nact) {
                                        break;
                                    }
                                }
                                *iact.offset(k as isize) = isave;
                                *vmultc.offset(k as isize) = vsave;
                            }
                            temp = 0.0f64;
                            i__1 = *n;
                            i__ = 1 as libc::c_int;
                            while i__ <= i__1 {
                                temp += *z__.offset((i__ + nact * z_dim1) as isize)
                                    * *a.offset((i__ + kk * a_dim1) as isize);
                                i__ += 1;
                            }
                            if temp == 0.0f64 {
                                break;
                            }
                            *zdota.offset(nact as isize) = temp;
                            *vmultc.offset(icon as isize) = 0.0f64;
                            *vmultc.offset(nact as isize) = ratio;
                        }
                        *iact.offset(icon as isize) = *iact.offset(nact as isize);
                        *iact.offset(nact as isize) = kk;
                        if mcon > *m && kk != mcon {
                            k = nact - 1 as libc::c_int;
                            sp = 0.0f64;
                            i__1 = *n;
                            i__ = 1 as libc::c_int;
                            while i__ <= i__1 {
                                sp += *z__.offset((i__ + k * z_dim1) as isize)
                                    * *a.offset((i__ + kk * a_dim1) as isize);
                                i__ += 1;
                            }
                            d__1 = *zdota.offset(nact as isize);
                            temp = (sp * sp + d__1 * d__1).sqrt();
                            alpha = *zdota.offset(nact as isize) / temp;
                            beta = sp / temp;
                            *zdota.offset(nact as isize) = alpha * *zdota.offset(k as isize);
                            *zdota.offset(k as isize) = temp;
                            i__1 = *n;
                            i__ = 1 as libc::c_int;
                            while i__ <= i__1 {
                                temp = alpha * *z__.offset((i__ + nact * z_dim1) as isize)
                                    + beta * *z__.offset((i__ + k * z_dim1) as isize);
                                *z__.offset((i__ + nact * z_dim1) as isize) = alpha
                                    * *z__.offset((i__ + k * z_dim1) as isize)
                                    - beta * *z__.offset((i__ + nact * z_dim1) as isize);
                                *z__.offset((i__ + k * z_dim1) as isize) = temp;
                                i__ += 1;
                            }
                            *iact.offset(nact as isize) = *iact.offset(k as isize);
                            *iact.offset(k as isize) = kk;
                            temp = *vmultc.offset(k as isize);
                            *vmultc.offset(k as isize) = *vmultc.offset(nact as isize);
                            *vmultc.offset(nact as isize) = temp;
                        }
                        if mcon > *m {
                            current_block = 15623375721314334080;
                        } else {
                            kk = *iact.offset(nact as isize);
                            temp = 0.0f64;
                            i__1 = *n;
                            i__ = 1 as libc::c_int;
                            while i__ <= i__1 {
                                temp += *sdirn.offset(i__ as isize)
                                    * *a.offset((i__ + kk * a_dim1) as isize);
                                i__ += 1;
                            }
                            temp += -1.0f64;
                            temp /= *zdota.offset(nact as isize);
                            i__1 = *n;
                            i__ = 1 as libc::c_int;
                            while i__ <= i__1 {
                                *sdirn.offset(i__ as isize) -=
                                    temp * *z__.offset((i__ + nact * z_dim1) as isize);
                                i__ += 1;
                            }
                            current_block = 9824026811267781195;
                        }
                    }
                    match current_block {
                        15623375721314334080 => {
                            temp = 1.0f64 / *zdota.offset(nact as isize);
                            i__1 = *n;
                            i__ = 1 as libc::c_int;
                            while i__ <= i__1 {
                                *sdirn.offset(i__ as isize) =
                                    temp * *z__.offset((i__ + nact * z_dim1) as isize);
                                i__ += 1;
                            }
                        }
                        _ => {}
                    }
                    dd = *rho * *rho;
                    sd = 0.0f64;
                    ss = 0.0f64;
                    i__1 = *n;
                    i__ = 1 as libc::c_int;
                    while i__ <= i__1 {
                        d__1 = *dx.offset(i__ as isize);
                        if (d__1).abs() >= *rho * 1e-6f32 as libc::c_double {
                            d__2 = *dx.offset(i__ as isize);
                            dd -= d__2 * d__2;
                        }
                        sd += *dx.offset(i__ as isize) * *sdirn.offset(i__ as isize);
                        d__1 = *sdirn.offset(i__ as isize);
                        ss += d__1 * d__1;
                        i__ += 1;
                    }
                    if dd <= 0.0f64 {
                        break;
                    }
                    temp = (ss * dd).sqrt();
                    if (sd).abs() >= temp * 1e-6f32 as libc::c_double {
                        temp = (ss * dd + sd * sd).sqrt();
                    }
                    stpful = dd / (temp + sd);
                    step = stpful;
                    if mcon == *m {
                        acca = step + resmax * 0.1f64;
                        accb = step + resmax * 0.2f64;
                        if step >= acca || acca >= accb {
                            current_block = 11188143500741601598;
                            continue 'c_8601;
                        }
                        step = if step <= resmax { step } else { resmax };
                    }
                    if nlopt_isinf(step) != 0 {
                        return NLOPT_ROUNDOFF_LIMITED;
                    }
                    i__1 = *n;
                    i__ = 1 as libc::c_int;
                    while i__ <= i__1 {
                        *dxnew.offset(i__ as isize) =
                            *dx.offset(i__ as isize) + step * *sdirn.offset(i__ as isize);
                        i__ += 1;
                    }
                    if mcon == *m {
                        resold = resmax;
                        resmax = 0.0f64;
                        i__1 = nact;
                        k = 1 as libc::c_int;
                        while k <= i__1 {
                            kk = *iact.offset(k as isize);
                            temp = *b.offset(kk as isize);
                            i__2 = *n;
                            i__ = 1 as libc::c_int;
                            while i__ <= i__2 {
                                temp -= *a.offset((i__ + kk * a_dim1) as isize)
                                    * *dxnew.offset(i__ as isize);
                                i__ += 1;
                            }
                            resmax = if resmax >= temp { resmax } else { temp };
                            k += 1;
                        }
                    }
                    k = nact;
                    loop {
                        zdotw = 0.0f64;
                        zdwabs = 0.0f64;
                        i__1 = *n;
                        i__ = 1 as libc::c_int;
                        while i__ <= i__1 {
                            temp = *z__.offset((i__ + k * z_dim1) as isize)
                                * *dxnew.offset(i__ as isize);
                            zdotw += temp;
                            zdwabs += (temp).abs();
                            i__ += 1;
                        }
                        acca = zdwabs + (zdotw).abs() * 0.1f64;
                        accb = zdwabs + (zdotw).abs() * 0.2f64;
                        if zdwabs >= acca || acca >= accb {
                            zdotw = 0.0f64;
                        }
                        *vmultd.offset(k as isize) = zdotw / *zdota.offset(k as isize);
                        if !(k >= 2 as libc::c_int) {
                            break;
                        }
                        kk = *iact.offset(k as isize);
                        i__1 = *n;
                        i__ = 1 as libc::c_int;
                        while i__ <= i__1 {
                            *dxnew.offset(i__ as isize) -= *vmultd.offset(k as isize)
                                * *a.offset((i__ + kk * a_dim1) as isize);
                            i__ += 1;
                        }
                        k -= 1;
                    }
                    if mcon > *m {
                        d__1 = 0.0f64;
                        d__2 = *vmultd.offset(nact as isize);
                        *vmultd.offset(nact as isize) = if d__1 >= d__2 { d__1 } else { d__2 };
                    }
                    i__1 = *n;
                    i__ = 1 as libc::c_int;
                    while i__ <= i__1 {
                        *dxnew.offset(i__ as isize) =
                            *dx.offset(i__ as isize) + step * *sdirn.offset(i__ as isize);
                        i__ += 1;
                    }
                    if mcon > nact {
                        kl = nact + 1 as libc::c_int;
                        i__1 = mcon;
                        k = kl;
                        while k <= i__1 {
                            kk = *iact.offset(k as isize);
                            sum = resmax - *b.offset(kk as isize);
                            d__1 = *b.offset(kk as isize);
                            sumabs = resmax + (d__1).abs();
                            i__2 = *n;
                            i__ = 1 as libc::c_int;
                            while i__ <= i__2 {
                                temp = *a.offset((i__ + kk * a_dim1) as isize)
                                    * *dxnew.offset(i__ as isize);
                                sum += temp;
                                sumabs += (temp).abs();
                                i__ += 1;
                            }
                            acca = sumabs + (sum).abs() * 0.1f32 as libc::c_double;
                            accb = sumabs + (sum).abs() * 0.2f32 as libc::c_double;
                            if sumabs >= acca || acca >= accb {
                                sum = 0.0f32 as libc::c_double;
                            }
                            *vmultd.offset(k as isize) = sum;
                            k += 1;
                        }
                    }
                    ratio = 1.0f64;
                    icon = 0 as libc::c_int;
                    i__1 = mcon;
                    k = 1 as libc::c_int;
                    while k <= i__1 {
                        if *vmultd.offset(k as isize) < 0.0f64 {
                            temp = *vmultc.offset(k as isize)
                                / (*vmultc.offset(k as isize) - *vmultd.offset(k as isize));
                            if temp < ratio {
                                ratio = temp;
                                icon = k;
                            }
                        }
                        k += 1;
                    }
                    temp = 1.0f64 - ratio;
                    i__1 = *n;
                    i__ = 1 as libc::c_int;
                    while i__ <= i__1 {
                        *dx.offset(i__ as isize) =
                            temp * *dx.offset(i__ as isize) + ratio * *dxnew.offset(i__ as isize);
                        i__ += 1;
                    }
                    i__1 = mcon;
                    k = 1 as libc::c_int;
                    while k <= i__1 {
                        d__1 = 0.0f64;
                        d__2 =
                            temp * *vmultc.offset(k as isize) + ratio * *vmultd.offset(k as isize);
                        *vmultc.offset(k as isize) = if d__1 >= d__2 { d__1 } else { d__2 };
                        k += 1;
                    }
                    if mcon == *m {
                        resmax = resold + ratio * (resmax - resold);
                    }
                    if icon > 0 as libc::c_int {
                        continue;
                    }
                    if step == stpful {
                        break 'c_8601;
                    } else {
                        current_block = 11188143500741601598;
                        continue 'c_8601;
                    }
                }
                if mcon == *m {
                    current_block = 11188143500741601598;
                    continue;
                }
                *ifull = 0 as libc::c_int;
                break;
            }
        }
    }
    return NLOPT_SUCCESS;
}
