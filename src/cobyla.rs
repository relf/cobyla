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
    clippy::unnecessary_cast
)]

//  {
//     // pub type _IO_wide_data;
//     // pub type _IO_codecvt;
//     // pub type _IO_marker;
//     // fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
//     // fn free(__ptr: *mut libc::c_void);
//     // fn memset(_: *mut libc::c_void, _: libc::c_int, _: libc::c_ulong) -> *mut libc::c_void;
//     // fn __errno_location() -> *mut libc::c_int;
//     // fn sqrt(_: libc::c_double) -> libc::c_double;
//     // fn fabs(_: libc::c_double) -> libc::c_double;
//     // static mut stdout: *mut FILE;
//     // static mut stderr: *mut FILE;
//     // fn fprintf(_: *mut FILE, _: *const libc::c_char, _: ...) -> libc::c_int;
// }
pub type size_t = libc::c_ulong;
pub type __off_t = libc::c_long;
pub type __off64_t = libc::c_long;
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
pub type cobyla_calcfc = unsafe fn(
    libc::c_long,
    libc::c_long,
    *const libc::c_double,
    *mut libc::c_double,
    *mut libc::c_void,
) -> libc::c_double;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _cobyla_context {
    pub n: libc::c_long,
    pub m: libc::c_long,
    pub iprint: libc::c_long,
    pub maxfun: libc::c_long,
    pub nfvals: libc::c_long,
    pub rhobeg: libc::c_double,
    pub rhoend: libc::c_double,
    pub iact: *mut libc::c_long,
    pub con: *mut libc::c_double,
    pub sim: *mut libc::c_double,
    pub simi: *mut libc::c_double,
    pub datmat: *mut libc::c_double,
    pub a: *mut libc::c_double,
    pub vsig: *mut libc::c_double,
    pub veta: *mut libc::c_double,
    pub sigbar: *mut libc::c_double,
    pub dx: *mut libc::c_double,
    pub w: *mut libc::c_double,
    pub parmu: libc::c_double,
    pub parsig: libc::c_double,
    pub prerec: libc::c_double,
    pub prerem: libc::c_double,
    pub rho: libc::c_double,
    pub f: libc::c_double,
    pub ibrnch: libc::c_long,
    pub iflag: libc::c_long,
    pub ifull: libc::c_long,
    pub jdrop: libc::c_long,
    pub status: libc::c_int,
}

pub type cobyla_context_t = _cobyla_context;

impl Default for cobyla_context_t {
    fn default() -> Self {
        cobyla_context_t {
            n: 0,
            m: 0,
            iprint: 0,
            maxfun: 0,
            nfvals: 0,
            rhobeg: 0.,
            rhoend: 0.,
            iact: 0 as *mut libc::c_long,
            con: 0 as *mut libc::c_double,
            sim: 0 as *mut libc::c_double,
            simi: 0 as *mut libc::c_double,
            datmat: 0 as *mut libc::c_double,
            a: 0 as *mut libc::c_double,
            vsig: 0 as *mut libc::c_double,
            veta: 0 as *mut libc::c_double,
            sigbar: 0 as *mut libc::c_double,
            dx: 0 as *mut libc::c_double,
            w: 0 as *mut libc::c_double,
            parmu: 0.,
            parsig: 0.,
            prerec: 0.,
            prerem: 0.,
            rho: 0.,
            f: 0.,
            ibrnch: 0,
            iflag: 0,
            ifull: 0,
            jdrop: 0,
            status: 0,
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub unsafe fn raw_cobyla(
    mut n: libc::c_long,
    mut m: libc::c_long,
    mut calcfc: Option<cobyla_calcfc>,
    mut calcfc_data: *mut libc::c_void,
    mut x: *mut libc::c_double,
    mut rhobeg: libc::c_double,
    mut rhoend: libc::c_double,
    mut iprint: libc::c_long,
    mut maxfun: *mut libc::c_long,
    mut w: *mut libc::c_double,
    mut iact: *mut libc::c_long,
) -> libc::c_int {
    let mut mpp: libc::c_long = m + 2 as libc::c_int as libc::c_long;
    let mut con: *mut libc::c_double = w;
    let mut sim: *mut libc::c_double = con.offset(mpp as isize);
    let mut simi: *mut libc::c_double = sim.offset((n * n) as isize).offset(n as isize);
    let mut datmat: *mut libc::c_double = simi.offset((n * n) as isize);
    let mut a: *mut libc::c_double = datmat.offset((n * mpp) as isize).offset(mpp as isize);
    let mut vsig: *mut libc::c_double = a.offset((m * n) as isize).offset(n as isize);
    let mut veta: *mut libc::c_double = vsig.offset(n as isize);
    let mut sigbar: *mut libc::c_double = veta.offset(n as isize);
    let mut dx: *mut libc::c_double = sigbar.offset(n as isize);
    let mut work: *mut libc::c_double = dx.offset(n as isize);
    return cobylb(
        n,
        m,
        calcfc,
        calcfc_data,
        x,
        rhobeg,
        rhoend,
        iprint,
        maxfun,
        con,
        sim,
        simi,
        datmat,
        a,
        vsig,
        veta,
        sigbar,
        dx,
        work,
        iact,
    );
}
#[no_mangle]
pub unsafe fn cobyla_create(
    mut n: libc::c_long,
    mut m: libc::c_long,
    mut rhobeg: libc::c_double,
    mut rhoend: libc::c_double,
    mut iprint: libc::c_long,
    mut maxfun: libc::c_long,
) -> *mut cobyla_context_t {
    let mut ctx: *mut cobyla_context_t = 0 as *mut cobyla_context_t;
    let mut size: libc::c_long = 0;
    let mut offset1: libc::c_long = 0;
    let mut offset2: libc::c_long = 0;
    let mut mpp: libc::c_long = 0;
    if n < 1 as libc::c_int as libc::c_long
        || m < 0 as libc::c_int as libc::c_long
        || rhobeg < rhoend
        || rhoend <= 0 as libc::c_int as libc::c_double
        || maxfun < 1 as libc::c_int as libc::c_long
    {
        // *__errno_location() = 22 as libc::c_int;
        return 0 as *mut cobyla_context_t;
    }
    size = ::std::mem::size_of::<cobyla_context_t>() as libc::c_ulong as libc::c_long;
    offset1 = (::std::mem::size_of::<libc::c_long>() as libc::c_ulong)
        .wrapping_sub(1 as libc::c_int as libc::c_ulong)
        .wrapping_add(size as libc::c_ulong)
        .wrapping_div(::std::mem::size_of::<libc::c_long>() as libc::c_ulong)
        .wrapping_mul(::std::mem::size_of::<libc::c_long>() as libc::c_ulong)
        as libc::c_long;
    size = (offset1 as libc::c_ulong).wrapping_add(
        ((m + 1 as libc::c_int as libc::c_long) as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<libc::c_long>() as libc::c_ulong),
    ) as libc::c_long;
    offset2 = (::std::mem::size_of::<libc::c_double>() as libc::c_ulong)
        .wrapping_sub(1 as libc::c_int as libc::c_ulong)
        .wrapping_add(size as libc::c_ulong)
        .wrapping_div(::std::mem::size_of::<libc::c_double>() as libc::c_ulong)
        .wrapping_mul(::std::mem::size_of::<libc::c_double>() as libc::c_ulong)
        as libc::c_long;
    size = (offset2 as libc::c_ulong).wrapping_add(
        ((n * (3 as libc::c_int as libc::c_long * n
            + 2 as libc::c_int as libc::c_long * m
            + 11 as libc::c_int as libc::c_long)
            + 4 as libc::c_int as libc::c_long * m
            + 6 as libc::c_int as libc::c_long) as libc::c_ulong)
            .wrapping_mul(::std::mem::size_of::<libc::c_double>() as libc::c_ulong),
    ) as libc::c_long;
    // ctx = malloc(size as libc::c_ulong) as *mut cobyla_context_t;
    ctx = Box::into_raw(Box::default());
    if ctx.is_null() {
        return 0 as *mut cobyla_context_t;
    }
    // memset(
    //     ctx as *mut libc::c_void,
    //     0 as libc::c_int,
    //     size as libc::c_ulong,
    // );
    (*ctx).n = n;
    (*ctx).m = m;
    (*ctx).nfvals = 0 as libc::c_int as libc::c_long;
    (*ctx).status = 1 as libc::c_int;
    (*ctx).iprint = iprint;
    (*ctx).maxfun = maxfun;
    (*ctx).rhobeg = rhobeg;
    (*ctx).rhoend = rhoend;
    mpp = m + 2 as libc::c_int as libc::c_long;
    let ref mut fresh0 = (*ctx).iact;
    *fresh0 = (ctx as *mut libc::c_char).offset(offset1 as isize) as *mut libc::c_long;
    let ref mut fresh1 = (*ctx).con;
    *fresh1 = (ctx as *mut libc::c_char).offset(offset2 as isize) as *mut libc::c_double;
    let ref mut fresh2 = (*ctx).sim;
    *fresh2 = ((*ctx).con).offset(mpp as isize);
    let ref mut fresh3 = (*ctx).simi;
    *fresh3 = ((*ctx).sim).offset((n * n) as isize).offset(n as isize);
    let ref mut fresh4 = (*ctx).datmat;
    *fresh4 = ((*ctx).simi).offset((n * n) as isize);
    let ref mut fresh5 = (*ctx).a;
    *fresh5 = ((*ctx).datmat)
        .offset((n * mpp) as isize)
        .offset(mpp as isize);
    let ref mut fresh6 = (*ctx).vsig;
    *fresh6 = ((*ctx).a).offset((m * n) as isize).offset(n as isize);
    let ref mut fresh7 = (*ctx).veta;
    *fresh7 = ((*ctx).vsig).offset(n as isize);
    let ref mut fresh8 = (*ctx).sigbar;
    *fresh8 = ((*ctx).veta).offset(n as isize);
    let ref mut fresh9 = (*ctx).dx;
    *fresh9 = ((*ctx).sigbar).offset(n as isize);
    let ref mut fresh10 = (*ctx).w;
    *fresh10 = ((*ctx).dx).offset(n as isize);
    return ctx;
}
#[no_mangle]
pub unsafe fn cobyla_delete(mut ctx: *mut cobyla_context_t) {
    if !ctx.is_null() {
        //free(ctx as *mut libc::c_void);
        drop(Box::from_raw(ctx));
    }
}
#[no_mangle]
pub unsafe fn cobyla_restart(mut ctx: *mut cobyla_context_t) -> libc::c_int {
    if ctx.is_null() {
        // *__errno_location() = 14 as libc::c_int;
        return -(3 as libc::c_int);
    }
    (*ctx).nfvals = 0 as libc::c_int as libc::c_long;
    (*ctx).status = 1 as libc::c_int;
    return (*ctx).status;
}
#[no_mangle]
pub unsafe fn cobyla_get_status(mut ctx: *const cobyla_context_t) -> libc::c_int {
    if ctx.is_null() {
        // *__errno_location() = 14 as libc::c_int;
        return -(3 as libc::c_int);
    }
    return (*ctx).status;
}
#[no_mangle]
pub unsafe fn cobyla_get_nevals(mut ctx: *const cobyla_context_t) -> libc::c_long {
    if ctx.is_null() {
        // *__errno_location() = 14 as libc::c_int;
        return -(1 as libc::c_int) as libc::c_long;
    }
    return (*ctx).nfvals;
}
#[no_mangle]
pub unsafe fn cobyla_get_rho(mut ctx: *const cobyla_context_t) -> libc::c_double {
    if ctx.is_null() {
        // *__errno_location() = 14 as libc::c_int;
        return -(1 as libc::c_int) as libc::c_double;
    }
    return (*ctx).rho;
}
#[no_mangle]
pub unsafe fn cobyla_get_last_f(mut ctx: *const cobyla_context_t) -> libc::c_double {
    if ctx.is_null() {
        // *__errno_location() = 14 as libc::c_int;
        return -(1 as libc::c_int) as libc::c_double;
    }
    return (*ctx).f;
}
#[no_mangle]
pub unsafe fn cobyla_iterate(
    mut ctx: *mut cobyla_context_t,
    mut f: libc::c_double,
    mut x: *mut libc::c_double,
    mut c: *mut libc::c_double,
) -> libc::c_int {
    let mut current_block: u64;
    let zero: libc::c_double = 0.0f64;
    let one: libc::c_double = 1.0f64;
    let alpha: libc::c_double = 0.25f64;
    let beta: libc::c_double = 2.1f64;
    let gamma: libc::c_double = 0.5f64;
    let delta: libc::c_double = 1.1f64;
    let mut parmu: libc::c_double = 0.;
    let mut parsig: libc::c_double = 0.;
    let mut prerec: libc::c_double = 0.;
    let mut prerem: libc::c_double = 0.;
    let mut rho: libc::c_double = 0.;
    let mut barmu: libc::c_double = 0.;
    let mut cvmaxm: libc::c_double = 0.;
    let mut cvmaxp: libc::c_double = 0.;
    let mut dxsign: libc::c_double = 0.;
    let mut edgmax: libc::c_double = 0.;
    let mut error: libc::c_double = 0.;
    let mut pareta: libc::c_double = 0.;
    let mut phi: libc::c_double = 0.;
    let mut phimin: libc::c_double = 0.;
    let mut ratio: libc::c_double = 0.;
    let mut resmax: libc::c_double = 0.;
    let mut resnew: libc::c_double = 0.;
    let mut sum: libc::c_double = 0.;
    let mut temp: libc::c_double = 0.;
    let mut tempa: libc::c_double = 0.;
    let mut trured: libc::c_double = 0.;
    let mut vmnew: libc::c_double = 0.;
    let mut vmold: libc::c_double = 0.;
    let mut ibrnch: libc::c_long = 0;
    let mut iflag: libc::c_long = 0;
    let mut ifull: libc::c_long = 0;
    let mut jdrop: libc::c_long = 0;
    let mut nfvals: libc::c_long = 0;
    let mut maxfun: libc::c_long = 0;
    let mut i: libc::c_long = 0;
    let mut j: libc::c_long = 0;
    let mut k: libc::c_long = 0;
    let mut l: libc::c_long = 0;
    let mut mp: libc::c_long = 0;
    let mut mpp: libc::c_long = 0;
    let mut np: libc::c_long = 0;
    let mut nbest: libc::c_long = 0;
    let mut n: libc::c_long = 0;
    let mut m: libc::c_long = 0;
    let mut iprint: libc::c_long = 0;
    let mut rhobeg: libc::c_double = 0.;
    let mut rhoend: libc::c_double = 0.;
    let mut con: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut sim: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut simi: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut datmat: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut a: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut vsig: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut veta: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut sigbar: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut dx: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut w: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut iact: *mut libc::c_long = 0 as *mut libc::c_long;
    let mut status: libc::c_int = 0;
    if ctx.is_null() {
        // *__errno_location() = 14 as libc::c_int;
        return -(3 as libc::c_int);
    }
    n = (*ctx).n;
    m = (*ctx).m;
    iprint = (*ctx).iprint;
    maxfun = (*ctx).maxfun;
    nfvals = (*ctx).nfvals;
    rhobeg = (*ctx).rhobeg;
    rhoend = (*ctx).rhoend;
    iact = (*ctx).iact;
    con = (*ctx).con;
    sim = (*ctx).sim;
    simi = (*ctx).simi;
    datmat = (*ctx).datmat;
    a = (*ctx).a;
    vsig = (*ctx).vsig;
    veta = (*ctx).veta;
    sigbar = (*ctx).sigbar;
    dx = (*ctx).dx;
    w = (*ctx).w;
    status = (*ctx).status;
    if x.is_null() || c.is_null() && m > 0 as libc::c_int as libc::c_long {
        // *__errno_location() = 14 as libc::c_int;
        (*ctx).status = -(3 as libc::c_int);
        return -(3 as libc::c_int);
    }
    if status != 1 as libc::c_int || nfvals < 0 as libc::c_int as libc::c_long {
        // *__errno_location() = 22 as libc::c_int;
        (*ctx).status = -(4 as libc::c_int);
        return (*ctx).status;
    }
    np = n + 1 as libc::c_int as libc::c_long;
    mp = m + 1 as libc::c_int as libc::c_long;
    mpp = m + 2 as libc::c_int as libc::c_long;
    if nfvals == 0 as libc::c_int as libc::c_long {
        status = 2 as libc::c_int;
        prerec = zero;
        prerem = zero;
        parsig = zero;
        iflag = 0 as libc::c_int as libc::c_long;
        rho = rhobeg;
        parmu = zero;
        jdrop = np;
        ibrnch = 0 as libc::c_int as libc::c_long;
        temp = one / rho;
        i = 1 as libc::c_int as libc::c_long;
        while i <= n {
            *sim.offset(
                (i - 1 as libc::c_int as libc::c_long + n * (np - 1 as libc::c_int as libc::c_long))
                    as isize,
            ) = *x.offset((i - 1 as libc::c_int as libc::c_long) as isize);
            j = 1 as libc::c_int as libc::c_long;
            while j <= n {
                *sim.offset(
                    (i - 1 as libc::c_int as libc::c_long
                        + n * (j - 1 as libc::c_int as libc::c_long)) as isize,
                ) = zero;
                *simi.offset(
                    (i - 1 as libc::c_int as libc::c_long
                        + n * (j - 1 as libc::c_int as libc::c_long)) as isize,
                ) = zero;
                j += 1;
            }
            *sim.offset(
                (i - 1 as libc::c_int as libc::c_long + n * (i - 1 as libc::c_int as libc::c_long))
                    as isize,
            ) = rho;
            *simi.offset(
                (i - 1 as libc::c_int as libc::c_long + n * (i - 1 as libc::c_int as libc::c_long))
                    as isize,
            ) = temp;
            i += 1;
        }
        if iprint >= 2 as libc::c_int as libc::c_long {
            println!(
                "\n   The initial value of RHO is {}  and PARMU is set to zero.\n",
                rho,
            );
        }
        current_block = 14453151562619017203;
    } else {
        parmu = (*ctx).parmu;
        parsig = (*ctx).parsig;
        prerec = (*ctx).prerec;
        prerem = (*ctx).prerem;
        rho = (*ctx).rho;
        ibrnch = (*ctx).ibrnch;
        iflag = (*ctx).iflag;
        ifull = (*ctx).ifull;
        jdrop = (*ctx).jdrop;
        current_block = 2561692241959298350;
    }
    'c_12387: loop {
        match current_block {
            14453151562619017203 => {
                if nfvals >= maxfun
                    && nfvals > 0 as libc::c_int as libc::c_long
                    && iprint > 0 as libc::c_int as libc::c_long
                {
                    status = -(2 as libc::c_int);
                    println!(
                        "Return from subroutine COBYLA because {:?}.\n",
                        cobyla_reason(status),
                    );
                    current_block = 2880604979707412946;
                    break;
                } else if status == 2 as libc::c_int {
                    status = 1 as libc::c_int;
                    current_block = 2561692241959298350;
                } else {
                    status = 1 as libc::c_int;
                    current_block = 14015749458106396694;
                    break;
                }
            }
            _ => {
                nfvals += 1;
                resmax = zero;
                k = 0 as libc::c_int as libc::c_long;
                while k < m {
                    temp = *c.offset(k as isize);
                    *con.offset(k as isize) = temp;
                    if resmax < -temp {
                        resmax = -temp;
                    }
                    k += 1;
                }
                if nfvals == iprint - 1 as libc::c_int as libc::c_long
                    || iprint == 3 as libc::c_int as libc::c_long
                {
                    print_calcfc(
                        n,
                        nfvals,
                        f,
                        resmax,
                        &mut *x.offset(0) as *mut libc::c_double as *const libc::c_double,
                    );
                }
                *con.offset((mp - 1 as libc::c_int as libc::c_long) as isize) = f;
                *con.offset((mpp - 1 as libc::c_int as libc::c_long) as isize) = resmax;
                if ibrnch == 1 as libc::c_int as libc::c_long {
                    vmold = *datmat.offset(
                        (mp - 1 as libc::c_int as libc::c_long
                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                            as isize,
                    ) + parmu
                        * *datmat.offset(
                            (mpp - 1 as libc::c_int as libc::c_long
                                + mpp * (np - 1 as libc::c_int as libc::c_long))
                                as isize,
                        );
                    vmnew = f + parmu * resmax;
                    trured = vmold - vmnew;
                    if parmu == zero
                        && f == *datmat.offset(
                            (mp - 1 as libc::c_int as libc::c_long
                                + mpp * (np - 1 as libc::c_int as libc::c_long))
                                as isize,
                        )
                    {
                        prerem = prerec;
                        trured = *datmat.offset(
                            (mpp - 1 as libc::c_int as libc::c_long
                                + mpp * (np - 1 as libc::c_int as libc::c_long))
                                as isize,
                        ) - resmax;
                    }
                    ratio = if trured <= zero { one } else { zero };
                    jdrop = 0 as libc::c_int as libc::c_long;
                    j = 1 as libc::c_int as libc::c_long;
                    while j <= n {
                        temp = zero;
                        i = 1 as libc::c_int as libc::c_long;
                        while i <= n {
                            temp += *simi.offset(
                                (j - 1 as libc::c_int as libc::c_long
                                    + n * (i - 1 as libc::c_int as libc::c_long))
                                    as isize,
                            ) * *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                            i += 1;
                        }
                        temp = temp.abs();
                        if temp > ratio {
                            jdrop = j;
                            ratio = temp;
                        }
                        *sigbar.offset((j - 1 as libc::c_int as libc::c_long) as isize) =
                            temp * *vsig.offset((j - 1 as libc::c_int as libc::c_long) as isize);
                        j += 1;
                    }
                    edgmax = delta * rho;
                    l = 0 as libc::c_int as libc::c_long;
                    j = 1 as libc::c_int as libc::c_long;
                    while j <= n {
                        if *sigbar.offset((j - 1 as libc::c_int as libc::c_long) as isize) >= parsig
                            || *sigbar.offset((j - 1 as libc::c_int as libc::c_long) as isize)
                                >= *vsig.offset((j - 1 as libc::c_int as libc::c_long) as isize)
                        {
                            temp = *veta.offset((j - 1 as libc::c_int as libc::c_long) as isize);
                            if trured > zero {
                                temp = zero;
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    let mut tempb: libc::c_double = *dx
                                        .offset((i - 1 as libc::c_int as libc::c_long) as isize)
                                        - *sim.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (j - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                    temp += tempb * tempb;
                                    i += 1;
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
                    if l > 0 as libc::c_int as libc::c_long {
                        jdrop = l;
                    }
                    if jdrop == 0 as libc::c_int as libc::c_long {
                        current_block = 12414752556692412193;
                    } else {
                        temp = zero;
                        i = 1 as libc::c_int as libc::c_long;
                        while i <= n {
                            *sim.offset(
                                (i - 1 as libc::c_int as libc::c_long
                                    + n * (jdrop - 1 as libc::c_int as libc::c_long))
                                    as isize,
                            ) = *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                            temp += *simi.offset(
                                (jdrop - 1 as libc::c_int as libc::c_long
                                    + n * (i - 1 as libc::c_int as libc::c_long))
                                    as isize,
                            ) * *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                            i += 1;
                        }
                        i = 1 as libc::c_int as libc::c_long;
                        while i <= n {
                            *simi.offset(
                                (jdrop - 1 as libc::c_int as libc::c_long
                                    + n * (i - 1 as libc::c_int as libc::c_long))
                                    as isize,
                            ) = *simi.offset(
                                (jdrop - 1 as libc::c_int as libc::c_long
                                    + n * (i - 1 as libc::c_int as libc::c_long))
                                    as isize,
                            ) / temp;
                            i += 1;
                        }
                        j = 1 as libc::c_int as libc::c_long;
                        while j <= n {
                            if j != jdrop {
                                temp = zero;
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    temp += *simi.offset(
                                        (j - 1 as libc::c_int as libc::c_long
                                            + n * (i - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) * *dx
                                        .offset((i - 1 as libc::c_int as libc::c_long) as isize);
                                    i += 1;
                                }
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    *simi.offset(
                                        (j - 1 as libc::c_int as libc::c_long
                                            + n * (i - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) -= temp
                                        * *simi.offset(
                                            (jdrop - 1 as libc::c_int as libc::c_long
                                                + n * (i - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                    i += 1;
                                }
                            }
                            j += 1;
                        }
                        k = 1 as libc::c_int as libc::c_long;
                        while k <= mpp {
                            *datmat.offset(
                                (k - 1 as libc::c_int as libc::c_long
                                    + mpp * (jdrop - 1 as libc::c_int as libc::c_long))
                                    as isize,
                            ) = *con.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                            k += 1;
                        }
                        if trured > zero && trured >= 0.1f64 * prerem {
                            current_block = 10213206415157612706;
                        } else {
                            current_block = 12414752556692412193;
                        }
                    }
                } else {
                    k = 1 as libc::c_int as libc::c_long;
                    while k <= mpp {
                        *datmat.offset(
                            (k - 1 as libc::c_int as libc::c_long
                                + mpp * (jdrop - 1 as libc::c_int as libc::c_long))
                                as isize,
                        ) = *con.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                        k += 1;
                    }
                    if !(nfvals > np) {
                        if jdrop <= n {
                            if *datmat.offset(
                                (mp - 1 as libc::c_int as libc::c_long
                                    + mpp * (np - 1 as libc::c_int as libc::c_long))
                                    as isize,
                            ) <= f
                            {
                                *x.offset((jdrop - 1 as libc::c_int as libc::c_long) as isize) =
                                    *sim.offset(
                                        (jdrop - 1 as libc::c_int as libc::c_long
                                            + n * (np - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                            } else {
                                *sim.offset(
                                    (jdrop - 1 as libc::c_int as libc::c_long
                                        + n * (np - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) = *x.offset((jdrop - 1 as libc::c_int as libc::c_long) as isize);
                                k = 1 as libc::c_int as libc::c_long;
                                while k <= mpp {
                                    *datmat.offset(
                                        (k - 1 as libc::c_int as libc::c_long
                                            + mpp * (jdrop - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = *datmat.offset(
                                        (k - 1 as libc::c_int as libc::c_long
                                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                    *datmat.offset(
                                        (k - 1 as libc::c_int as libc::c_long
                                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = *con
                                        .offset((k - 1 as libc::c_int as libc::c_long) as isize);
                                    k += 1;
                                }
                                k = 1 as libc::c_int as libc::c_long;
                                while k <= jdrop {
                                    temp = zero;
                                    *sim.offset(
                                        (jdrop - 1 as libc::c_int as libc::c_long
                                            + n * (k - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = -rho;
                                    i = k;
                                    while i <= jdrop {
                                        temp -= *simi.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (k - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                        i += 1;
                                    }
                                    *simi.offset(
                                        (jdrop - 1 as libc::c_int as libc::c_long
                                            + n * (k - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = temp;
                                    k += 1;
                                }
                            }
                        }
                        if nfvals <= n {
                            jdrop = nfvals;
                            *x.offset((jdrop - 1 as libc::c_int as libc::c_long) as isize) += rho;
                            current_block = 14453151562619017203;
                            continue;
                        }
                    }
                    ibrnch = 1 as libc::c_int as libc::c_long;
                    current_block = 10213206415157612706;
                }
                'c_12399: loop {
                    match current_block {
                        12414752556692412193 => {
                            if iflag == 0 as libc::c_int as libc::c_long {
                                ibrnch = 0 as libc::c_int as libc::c_long;
                                current_block = 10213206415157612706;
                            } else if rho > rhoend {
                                rho = 0.5f64 * rho;
                                if rho <= 1.5f64 * rhoend {
                                    rho = rhoend;
                                }
                                if parmu > zero {
                                    let mut cmin: libc::c_double = 0.;
                                    let mut cmax: libc::c_double = 0.;
                                    let mut denom: libc::c_double = 0.;
                                    denom = zero;
                                    k = 1 as libc::c_int as libc::c_long;
                                    while k <= mp {
                                        cmin = *datmat.offset(
                                            (k - 1 as libc::c_int as libc::c_long
                                                + mpp * (np - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                        cmax = cmin;
                                        i = 1 as libc::c_int as libc::c_long;
                                        while i <= n {
                                            temp = *datmat.offset(
                                                (k - 1 as libc::c_int as libc::c_long
                                                    + mpp * (i - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            );
                                            cmin = if cmin <= temp { cmin } else { temp };
                                            cmax = if cmax >= temp { cmax } else { temp };
                                            i += 1;
                                        }
                                        if k <= m && cmin < 0.5f64 * cmax {
                                            temp = (if cmax >= zero { cmax } else { zero }) - cmin;
                                            if denom <= zero {
                                                denom = temp;
                                            } else {
                                                denom = if denom <= temp { denom } else { temp };
                                            }
                                        }
                                        k += 1;
                                    }
                                    if denom == zero {
                                        parmu = zero;
                                    } else if cmax - cmin < parmu * denom {
                                        parmu = (cmax - cmin) / denom;
                                    }
                                }
                                if iprint >= 2 as libc::c_int as libc::c_long {
                                    println!(
                                        "\n   Reduction in RHO to {}  and PARMU ={}\n",
                                        rho, parmu,
                                    );
                                    if iprint == 2 as libc::c_int as libc::c_long {
                                        print_calcfc(
                                            n,
                                            nfvals,
                                            *datmat.offset(
                                                (mp - 1 as libc::c_int as libc::c_long
                                                    + mpp * (np - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            ),
                                            *datmat.offset(
                                                (mpp - 1 as libc::c_int as libc::c_long
                                                    + mpp * (np - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            ),
                                            &mut *sim.offset(
                                                (n * (np - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            )
                                                as *mut libc::c_double
                                                as *const libc::c_double,
                                        );
                                    }
                                }
                                current_block = 10213206415157612706;
                            } else {
                                if iprint >= 1 as libc::c_int as libc::c_long {
                                    println!("\n   Normal return from subroutine COBYLA\n");
                                }
                                status = 0 as libc::c_int;
                                if ifull == 1 as libc::c_int as libc::c_long {
                                    current_block = 18000475564953935197;
                                    break 'c_12387;
                                } else {
                                    current_block = 2880604979707412946;
                                    break 'c_12387;
                                }
                            }
                        }
                        _ => {
                            phimin = *datmat.offset(
                                (mp - 1 as libc::c_int as libc::c_long
                                    + mpp * (np - 1 as libc::c_int as libc::c_long))
                                    as isize,
                            ) + parmu
                                * *datmat.offset(
                                    (mpp - 1 as libc::c_int as libc::c_long
                                        + mpp * (np - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                );
                            nbest = np;
                            j = 1 as libc::c_int as libc::c_long;
                            while j <= n {
                                temp = *datmat.offset(
                                    (mp - 1 as libc::c_int as libc::c_long
                                        + mpp * (j - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) + parmu
                                    * *datmat.offset(
                                        (mpp - 1 as libc::c_int as libc::c_long
                                            + mpp * (j - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                if temp < phimin {
                                    nbest = j;
                                    phimin = temp;
                                } else if temp == phimin && parmu == zero {
                                    if *datmat.offset(
                                        (mpp - 1 as libc::c_int as libc::c_long
                                            + mpp * (j - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) < *datmat.offset(
                                        (mpp - 1 as libc::c_int as libc::c_long
                                            + mpp * (nbest - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) {
                                        nbest = j;
                                    }
                                }
                                j += 1;
                            }
                            if nbest <= n {
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= mpp {
                                    temp = *datmat.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                    *datmat.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = *datmat.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + mpp * (nbest - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                    *datmat.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + mpp * (nbest - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = temp;
                                    i += 1;
                                }
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    temp = *sim.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (nbest - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                    *sim.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (nbest - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = zero;
                                    *sim.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (np - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) += temp;
                                    tempa = zero;
                                    k = 1 as libc::c_int as libc::c_long;
                                    while k <= n {
                                        *sim.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (k - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) -= temp;
                                        tempa -= *simi.offset(
                                            (k - 1 as libc::c_int as libc::c_long
                                                + n * (i - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                        k += 1;
                                    }
                                    *simi.offset(
                                        (nbest - 1 as libc::c_int as libc::c_long
                                            + n * (i - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = tempa;
                                    i += 1;
                                }
                            }
                            error = zero;
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                j = 1 as libc::c_int as libc::c_long;
                                while j <= n {
                                    temp = if i == j { -one } else { zero };
                                    k = 1 as libc::c_int as libc::c_long;
                                    while k <= n {
                                        temp += *simi.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (k - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) * *sim.offset(
                                            (k - 1 as libc::c_int as libc::c_long
                                                + n * (j - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                        k += 1;
                                    }
                                    temp = temp.abs();
                                    error = if error >= temp { error } else { temp };
                                    j += 1;
                                }
                                i += 1;
                            }
                            if error > 0.1f64 {
                                status = -(1 as libc::c_int);
                                if iprint >= 1 as libc::c_int as libc::c_long {
                                    println!(
                                        "Return from subroutine COBYLA because {:?}.\n",
                                        cobyla_reason(status),
                                    );
                                }
                                current_block = 2880604979707412946;
                                break 'c_12387;
                            } else {
                                k = 1 as libc::c_int as libc::c_long;
                                while k <= mp {
                                    *con.offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                        -*datmat.offset(
                                            (k - 1 as libc::c_int as libc::c_long
                                                + mpp * (np - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                    j = 1 as libc::c_int as libc::c_long;
                                    while j <= n {
                                        *w.offset(
                                            (j - 1 as libc::c_int as libc::c_long) as isize,
                                        ) = *datmat.offset(
                                            (k - 1 as libc::c_int as libc::c_long
                                                + mpp * (j - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) + *con.offset(
                                            (k - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                        j += 1;
                                    }
                                    i = 1 as libc::c_int as libc::c_long;
                                    while i <= n {
                                        temp = zero;
                                        j = 1 as libc::c_int as libc::c_long;
                                        while j <= n {
                                            temp += *w.offset(
                                                (j - 1 as libc::c_int as libc::c_long) as isize,
                                            ) * *simi.offset(
                                                (j - 1 as libc::c_int as libc::c_long
                                                    + n * (i - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            );
                                            j += 1;
                                        }
                                        *a.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (k - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) = if k == mp { -temp } else { temp };
                                        i += 1;
                                    }
                                    k += 1;
                                }
                                iflag = 1 as libc::c_int as libc::c_long;
                                parsig = alpha * rho;
                                pareta = beta * rho;
                                j = 1 as libc::c_int as libc::c_long;
                                while j <= n {
                                    let mut wsig: libc::c_double = zero;
                                    let mut weta: libc::c_double = zero;
                                    i = 1 as libc::c_int as libc::c_long;
                                    while i <= n {
                                        wsig += *simi.offset(
                                            (j - 1 as libc::c_int as libc::c_long
                                                + n * (i - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) * *simi.offset(
                                            (j - 1 as libc::c_int as libc::c_long
                                                + n * (i - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                        weta += *sim.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (j - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) * *sim.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (j - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                        i += 1;
                                    }
                                    *vsig.offset((j - 1 as libc::c_int as libc::c_long) as isize) =
                                        one / wsig.sqrt();
                                    *veta.offset((j - 1 as libc::c_int as libc::c_long) as isize) =
                                        weta.sqrt();
                                    if *vsig.offset((j - 1 as libc::c_int as libc::c_long) as isize)
                                        < parsig
                                        || *veta
                                            .offset((j - 1 as libc::c_int as libc::c_long) as isize)
                                            > pareta
                                    {
                                        iflag = 0 as libc::c_int as libc::c_long;
                                    }
                                    j += 1;
                                }
                                if ibrnch == 1 as libc::c_int as libc::c_long
                                    || iflag == 1 as libc::c_int as libc::c_long
                                {
                                    let mut z: *mut libc::c_double = w;
                                    let mut zdota: *mut libc::c_double = z.offset((n * n) as isize);
                                    let mut vmc: *mut libc::c_double = zdota.offset(n as isize);
                                    let mut sdirn: *mut libc::c_double = vmc.offset(mp as isize);
                                    let mut dxnew: *mut libc::c_double = sdirn.offset(n as isize);
                                    let mut vmd: *mut libc::c_double = dxnew.offset(n as isize);
                                    trstlp(
                                        n,
                                        m,
                                        a as *const libc::c_double,
                                        con as *const libc::c_double,
                                        rho,
                                        dx,
                                        &mut ifull,
                                        iact,
                                        z,
                                        zdota,
                                        vmc,
                                        sdirn,
                                        dxnew,
                                        vmd,
                                    );
                                    if ifull == 0 as libc::c_int as libc::c_long {
                                        temp = zero;
                                        i = 1 as libc::c_int as libc::c_long;
                                        while i <= n {
                                            temp += *dx.offset(
                                                (i - 1 as libc::c_int as libc::c_long) as isize,
                                            ) * *dx.offset(
                                                (i - 1 as libc::c_int as libc::c_long) as isize,
                                            );
                                            i += 1;
                                        }
                                        if temp < 0.25f64 * rho * rho {
                                            ibrnch = 1 as libc::c_int as libc::c_long;
                                            current_block = 12414752556692412193;
                                            continue;
                                        }
                                    }
                                    resnew = zero;
                                    *con.offset((mp - 1 as libc::c_int as libc::c_long) as isize) =
                                        zero;
                                    sum = zero;
                                    k = 1 as libc::c_int as libc::c_long;
                                    while k <= mp {
                                        sum = *con.offset(
                                            (k - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                        i = 1 as libc::c_int as libc::c_long;
                                        while i <= n {
                                            sum -= *a.offset(
                                                (i - 1 as libc::c_int as libc::c_long
                                                    + n * (k - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            ) * *dx.offset(
                                                (i - 1 as libc::c_int as libc::c_long) as isize,
                                            );
                                            i += 1;
                                        }
                                        if k < mp {
                                            resnew = if resnew >= sum { resnew } else { sum };
                                        }
                                        k += 1;
                                    }
                                    prerec = *datmat.offset(
                                        (mpp - 1 as libc::c_int as libc::c_long
                                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) - resnew;
                                    barmu = if prerec > zero { sum / prerec } else { zero };
                                    if !(parmu < 1.5f64 * barmu) {
                                        break;
                                    }
                                    parmu = 2.0f64 * barmu;
                                    if iprint >= 2 as libc::c_int as libc::c_long {
                                        println!("\n   Increase in PARMU to {}\n", parmu);
                                    }
                                    phi = *datmat.offset(
                                        (mp - 1 as libc::c_int as libc::c_long
                                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) + parmu
                                        * *datmat.offset(
                                            (mpp - 1 as libc::c_int as libc::c_long
                                                + mpp * (np - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                    j = 1 as libc::c_int as libc::c_long;
                                    loop {
                                        if !(j <= n) {
                                            break 'c_12399;
                                        }
                                        temp = *datmat.offset(
                                            (mp - 1 as libc::c_int as libc::c_long
                                                + mpp * (j - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) + parmu
                                            * *datmat.offset(
                                                (mpp - 1 as libc::c_int as libc::c_long
                                                    + mpp * (j - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            );
                                        if temp < phi {
                                            current_block = 10213206415157612706;
                                            break;
                                        }
                                        if temp == phi && parmu == zero {
                                            if *datmat.offset(
                                                (mpp - 1 as libc::c_int as libc::c_long
                                                    + mpp * (j - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            ) < *datmat.offset(
                                                (mpp - 1 as libc::c_int as libc::c_long
                                                    + mpp * (np - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            ) {
                                                current_block = 10213206415157612706;
                                                break;
                                            }
                                        }
                                        j += 1;
                                    }
                                } else {
                                    jdrop = 0 as libc::c_int as libc::c_long;
                                    temp = pareta;
                                    j = 1 as libc::c_int as libc::c_long;
                                    while j <= n {
                                        if *veta
                                            .offset((j - 1 as libc::c_int as libc::c_long) as isize)
                                            > temp
                                        {
                                            jdrop = j;
                                            temp = *veta.offset(
                                                (j - 1 as libc::c_int as libc::c_long) as isize,
                                            );
                                        }
                                        j += 1;
                                    }
                                    if jdrop == 0 as libc::c_int as libc::c_long {
                                        j = 1 as libc::c_int as libc::c_long;
                                        while j <= n {
                                            if *vsig.offset(
                                                (j - 1 as libc::c_int as libc::c_long) as isize,
                                            ) < temp
                                            {
                                                jdrop = j;
                                                temp = *vsig.offset(
                                                    (j - 1 as libc::c_int as libc::c_long) as isize,
                                                );
                                            }
                                            j += 1;
                                        }
                                    }
                                    temp = gamma
                                        * rho
                                        * *vsig.offset(
                                            (jdrop - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                    i = 1 as libc::c_int as libc::c_long;
                                    while i <= n {
                                        *dx.offset(
                                            (i - 1 as libc::c_int as libc::c_long) as isize,
                                        ) = temp
                                            * *simi.offset(
                                                (jdrop - 1 as libc::c_int as libc::c_long
                                                    + n * (i - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            );
                                        i += 1;
                                    }
                                    cvmaxp = zero;
                                    cvmaxm = zero;
                                    sum = zero;
                                    k = 1 as libc::c_int as libc::c_long;
                                    while k <= mp {
                                        sum = zero;
                                        i = 1 as libc::c_int as libc::c_long;
                                        while i <= n {
                                            sum = sum
                                                + *a.offset(
                                                    (i - 1 as libc::c_int as libc::c_long
                                                        + n * (k - 1 as libc::c_int
                                                            as libc::c_long))
                                                        as isize,
                                                ) * *dx.offset(
                                                    (i - 1 as libc::c_int as libc::c_long) as isize,
                                                );
                                            i += 1;
                                        }
                                        if k < mp {
                                            temp = *datmat.offset(
                                                (k - 1 as libc::c_int as libc::c_long
                                                    + mpp * (np - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            );
                                            cvmaxp = if cvmaxp >= -sum - temp {
                                                cvmaxp
                                            } else {
                                                -sum - temp
                                            };
                                            cvmaxm = if cvmaxm >= sum - temp {
                                                cvmaxm
                                            } else {
                                                sum - temp
                                            };
                                        }
                                        k += 1;
                                    }
                                    if parmu * (cvmaxp - cvmaxm) > sum + sum {
                                        dxsign = -one;
                                    } else {
                                        dxsign = one;
                                    }
                                    temp = zero;
                                    i = 1 as libc::c_int as libc::c_long;
                                    while i <= n {
                                        *dx.offset(
                                            (i - 1 as libc::c_int as libc::c_long) as isize,
                                        ) *= dxsign;
                                        *sim.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (jdrop - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) = *dx.offset(
                                            (i - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                        temp += *simi.offset(
                                            (jdrop - 1 as libc::c_int as libc::c_long
                                                + n * (i - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) * *dx.offset(
                                            (i - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                        i += 1;
                                    }
                                    i = 1 as libc::c_int as libc::c_long;
                                    while i <= n {
                                        *simi.offset(
                                            (jdrop - 1 as libc::c_int as libc::c_long
                                                + n * (i - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) /= temp;
                                        i += 1;
                                    }
                                    j = 1 as libc::c_int as libc::c_long;
                                    while j <= n {
                                        if j != jdrop {
                                            temp = zero;
                                            i = 1 as libc::c_int as libc::c_long;
                                            while i <= n {
                                                temp += *simi.offset(
                                                    (j - 1 as libc::c_int as libc::c_long
                                                        + n * (i - 1 as libc::c_int
                                                            as libc::c_long))
                                                        as isize,
                                                ) * *dx.offset(
                                                    (i - 1 as libc::c_int as libc::c_long) as isize,
                                                );
                                                i += 1;
                                            }
                                            i = 1 as libc::c_int as libc::c_long;
                                            while i <= n {
                                                *simi.offset(
                                                    (j - 1 as libc::c_int as libc::c_long
                                                        + n * (i - 1 as libc::c_int
                                                            as libc::c_long))
                                                        as isize,
                                                ) -= temp
                                                    * *simi.offset(
                                                        (jdrop - 1 as libc::c_int as libc::c_long
                                                            + n * (i - 1 as libc::c_int
                                                                as libc::c_long))
                                                            as isize,
                                                    );
                                                i += 1;
                                            }
                                        }
                                        *x.offset(
                                            (j - 1 as libc::c_int as libc::c_long) as isize,
                                        ) = *sim.offset(
                                            (j - 1 as libc::c_int as libc::c_long
                                                + n * (np - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) + *dx.offset(
                                            (j - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                        j += 1;
                                    }
                                    current_block = 14453151562619017203;
                                    continue 'c_12387;
                                }
                            }
                        }
                    }
                }
                prerem = parmu * prerec - sum;
                i = 1 as libc::c_int as libc::c_long;
                while i <= n {
                    *x.offset((i - 1 as libc::c_int as libc::c_long) as isize) = *sim.offset(
                        (i - 1 as libc::c_int as libc::c_long
                            + n * (np - 1 as libc::c_int as libc::c_long))
                            as isize,
                    ) + *dx
                        .offset((i - 1 as libc::c_int as libc::c_long) as isize);
                    i += 1;
                }
                ibrnch = 1 as libc::c_int as libc::c_long;
                current_block = 14453151562619017203;
            }
        }
    }
    match current_block {
        2880604979707412946 => {
            i = 1 as libc::c_int as libc::c_long;
            while i <= n {
                *x.offset((i - 1 as libc::c_int as libc::c_long) as isize) = *sim.offset(
                    (i - 1 as libc::c_int as libc::c_long
                        + n * (np - 1 as libc::c_int as libc::c_long)) as isize,
                );
                i += 1;
            }
            f = *datmat.offset(
                (mp - 1 as libc::c_int as libc::c_long
                    + mpp * (np - 1 as libc::c_int as libc::c_long)) as isize,
            );
            resmax = *datmat.offset(
                (mpp - 1 as libc::c_int as libc::c_long
                    + mpp * (np - 1 as libc::c_int as libc::c_long)) as isize,
            );
            current_block = 18000475564953935197;
        }
        _ => {}
    }
    match current_block {
        18000475564953935197 => {
            if iprint >= 1 as libc::c_int as libc::c_long {
                print_calcfc(
                    n,
                    nfvals,
                    f,
                    resmax,
                    &mut *x.offset(0 as isize) as *mut libc::c_double as *const libc::c_double,
                );
            }
        }
        _ => {}
    }
    (*ctx).nfvals = nfvals;
    (*ctx).parmu = parmu;
    (*ctx).parsig = parsig;
    (*ctx).prerec = prerec;
    (*ctx).prerem = prerem;
    (*ctx).rho = rho;
    (*ctx).f = f;
    (*ctx).ibrnch = ibrnch;
    (*ctx).iflag = iflag;
    (*ctx).ifull = ifull;
    (*ctx).jdrop = jdrop;
    (*ctx).status = status;
    return status;
}

#[allow(clippy::too_many_arguments)]
unsafe fn cobylb(
    n: libc::c_long,
    m: libc::c_long,
    mut calcfc: Option<cobyla_calcfc>,
    mut calcfc_data: *mut libc::c_void,
    mut x: *mut libc::c_double,
    rhobeg: libc::c_double,
    rhoend: libc::c_double,
    iprint: libc::c_long,
    mut _maxfun: *mut libc::c_long,
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
    mut iact: *mut libc::c_long,
) -> libc::c_int {
    let mut current_block: u64;
    let zero: libc::c_double = 0.0f64;
    let one: libc::c_double = 1.0f64;
    let alpha: libc::c_double = 0.25f64;
    let beta: libc::c_double = 2.1f64;
    let gamma: libc::c_double = 0.5f64;
    let delta: libc::c_double = 1.1f64;
    let mut parmu: libc::c_double = 0.;
    let mut parsig: libc::c_double = 0.;
    let mut prerec: libc::c_double = 0.;
    let mut prerem: libc::c_double = 0.;
    let mut rho: libc::c_double = 0.;
    let mut barmu: libc::c_double = 0.;
    let mut cvmaxm: libc::c_double = 0.;
    let mut cvmaxp: libc::c_double = 0.;
    let mut dxsign: libc::c_double = 0.;
    let mut edgmax: libc::c_double = 0.;
    let mut error: libc::c_double = 0.;
    let mut pareta: libc::c_double = 0.;
    let mut phi: libc::c_double = 0.;
    let mut phimin: libc::c_double = 0.;
    let mut ratio: libc::c_double = 0.;
    let mut resmax: libc::c_double = 0.;
    let mut resnew: libc::c_double = 0.;
    let mut sum: libc::c_double = 0.;
    let mut temp: libc::c_double = 0.;
    let mut tempa: libc::c_double = 0.;
    let mut trured: libc::c_double = 0.;
    let mut vmnew: libc::c_double = 0.;
    let mut vmold: libc::c_double = 0.;
    let mut ibrnch: libc::c_long = 0;
    let mut iflag: libc::c_long = 0;
    let mut ifull: libc::c_long = 0;
    let mut jdrop: libc::c_long = 0;
    let mut nfvals: libc::c_long = 0;
    let mut maxfun: libc::c_long = 0;
    let mut i: libc::c_long = 0;
    let mut j: libc::c_long = 0;
    let mut k: libc::c_long = 0;
    let mut l: libc::c_long = 0;
    let mut mp: libc::c_long = 0;
    let mut mpp: libc::c_long = 0;
    let mut np: libc::c_long = 0;
    let mut nbest: libc::c_long = 0;
    let mut f: libc::c_double = 0.;
    let mut status: libc::c_int = 0;
    np = n + 1 as libc::c_int as libc::c_long;
    mp = m + 1 as libc::c_int as libc::c_long;
    mpp = m + 2 as libc::c_int as libc::c_long;
    maxfun = *_maxfun;
    nfvals = 0 as libc::c_int as libc::c_long;
    rho = rhobeg;
    parmu = zero;
    jdrop = np;
    ibrnch = 0 as libc::c_int as libc::c_long;
    temp = one / rho;
    i = 1 as libc::c_int as libc::c_long;
    while i <= n {
        *sim.offset(
            (i - 1 as libc::c_int as libc::c_long + n * (np - 1 as libc::c_int as libc::c_long))
                as isize,
        ) = *x.offset((i - 1 as libc::c_int as libc::c_long) as isize);
        j = 1 as libc::c_int as libc::c_long;
        while j <= n {
            *sim.offset(
                (i - 1 as libc::c_int as libc::c_long + n * (j - 1 as libc::c_int as libc::c_long))
                    as isize,
            ) = zero;
            *simi.offset(
                (i - 1 as libc::c_int as libc::c_long + n * (j - 1 as libc::c_int as libc::c_long))
                    as isize,
            ) = zero;
            j += 1;
        }
        *sim.offset(
            (i - 1 as libc::c_int as libc::c_long + n * (i - 1 as libc::c_int as libc::c_long))
                as isize,
        ) = rho;
        *simi.offset(
            (i - 1 as libc::c_int as libc::c_long + n * (i - 1 as libc::c_int as libc::c_long))
                as isize,
        ) = temp;
        i += 1;
    }
    if iprint >= 2 as libc::c_int as libc::c_long {
        println!(
            "\n   The initial value of RHO is {}  and PARMU is set to zero.\n",
            rho,
        );
    }
    'c_3166: loop {
        if nfvals >= maxfun
            && nfvals > 0 as libc::c_int as libc::c_long
            && iprint > 0 as libc::c_int as libc::c_long
        {
            status = -(2 as libc::c_int);
            println!(
                "Return from subroutine COBYLA because {:?}.\n",
                cobyla_reason(status),
            );
            current_block = 18034400475116118263;
            break;
        } else {
            f = calcfc.expect("non-null function pointer")(
                n,
                m,
                x as *const libc::c_double,
                con,
                calcfc_data,
            );
            nfvals += 1;
            resmax = zero;
            k = 0 as libc::c_int as libc::c_long;
            while k < m {
                temp = -*con.offset(k as isize);
                if resmax < temp {
                    resmax = temp;
                }
                k += 1;
            }
            if nfvals == iprint - 1 as libc::c_int as libc::c_long
                || iprint == 3 as libc::c_int as libc::c_long
            {
                print_calcfc(
                    n,
                    nfvals,
                    f,
                    resmax,
                    &mut *x.offset(0) as *mut libc::c_double as *const libc::c_double,
                );
            }
            *con.offset((mp - 1 as libc::c_int as libc::c_long) as isize) = f;
            *con.offset((mpp - 1 as libc::c_int as libc::c_long) as isize) = resmax;
            if ibrnch == 1 as libc::c_int as libc::c_long {
                vmold = *datmat.offset(
                    (mp - 1 as libc::c_int as libc::c_long
                        + mpp * (np - 1 as libc::c_int as libc::c_long))
                        as isize,
                ) + parmu
                    * *datmat.offset(
                        (mpp - 1 as libc::c_int as libc::c_long
                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                            as isize,
                    );
                vmnew = f + parmu * resmax;
                trured = vmold - vmnew;
                if parmu == zero
                    && f == *datmat.offset(
                        (mp - 1 as libc::c_int as libc::c_long
                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                            as isize,
                    )
                {
                    prerem = prerec;
                    trured = *datmat.offset(
                        (mpp - 1 as libc::c_int as libc::c_long
                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                            as isize,
                    ) - resmax;
                }
                ratio = if trured <= zero { one } else { zero };
                jdrop = 0 as libc::c_int as libc::c_long;
                j = 1 as libc::c_int as libc::c_long;
                while j <= n {
                    temp = zero;
                    i = 1 as libc::c_int as libc::c_long;
                    while i <= n {
                        temp += *simi.offset(
                            (j - 1 as libc::c_int as libc::c_long
                                + n * (i - 1 as libc::c_int as libc::c_long))
                                as isize,
                        ) * *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                        i += 1;
                    }
                    temp = temp.abs();
                    if temp > ratio {
                        jdrop = j;
                        ratio = temp;
                    }
                    *sigbar.offset((j - 1 as libc::c_int as libc::c_long) as isize) =
                        temp * *vsig.offset((j - 1 as libc::c_int as libc::c_long) as isize);
                    j += 1;
                }
                edgmax = delta * rho;
                l = 0 as libc::c_int as libc::c_long;
                j = 1 as libc::c_int as libc::c_long;
                while j <= n {
                    if *sigbar.offset((j - 1 as libc::c_int as libc::c_long) as isize) >= parsig
                        || *sigbar.offset((j - 1 as libc::c_int as libc::c_long) as isize)
                            >= *vsig.offset((j - 1 as libc::c_int as libc::c_long) as isize)
                    {
                        temp = *veta.offset((j - 1 as libc::c_int as libc::c_long) as isize);
                        if trured > zero {
                            temp = zero;
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                let mut tempb: libc::c_double = *dx
                                    .offset((i - 1 as libc::c_int as libc::c_long) as isize)
                                    - *sim.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (j - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                temp += tempb * tempb;
                                i += 1;
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
                if l > 0 as libc::c_int as libc::c_long {
                    jdrop = l;
                }
                if jdrop == 0 as libc::c_int as libc::c_long {
                    current_block = 12206176237338959610;
                } else {
                    temp = zero;
                    i = 1 as libc::c_int as libc::c_long;
                    while i <= n {
                        *sim.offset(
                            (i - 1 as libc::c_int as libc::c_long
                                + n * (jdrop - 1 as libc::c_int as libc::c_long))
                                as isize,
                        ) = *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                        temp += *simi.offset(
                            (jdrop - 1 as libc::c_int as libc::c_long
                                + n * (i - 1 as libc::c_int as libc::c_long))
                                as isize,
                        ) * *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                        i += 1;
                    }
                    i = 1 as libc::c_int as libc::c_long;
                    while i <= n {
                        *simi.offset(
                            (jdrop - 1 as libc::c_int as libc::c_long
                                + n * (i - 1 as libc::c_int as libc::c_long))
                                as isize,
                        ) = *simi.offset(
                            (jdrop - 1 as libc::c_int as libc::c_long
                                + n * (i - 1 as libc::c_int as libc::c_long))
                                as isize,
                        ) / temp;
                        i += 1;
                    }
                    j = 1 as libc::c_int as libc::c_long;
                    while j <= n {
                        if j != jdrop {
                            temp = zero;
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                temp += *simi.offset(
                                    (j - 1 as libc::c_int as libc::c_long
                                        + n * (i - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) * *dx
                                    .offset((i - 1 as libc::c_int as libc::c_long) as isize);
                                i += 1;
                            }
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                *simi.offset(
                                    (j - 1 as libc::c_int as libc::c_long
                                        + n * (i - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) -= temp
                                    * *simi.offset(
                                        (jdrop - 1 as libc::c_int as libc::c_long
                                            + n * (i - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                i += 1;
                            }
                        }
                        j += 1;
                    }
                    k = 1 as libc::c_int as libc::c_long;
                    while k <= mpp {
                        *datmat.offset(
                            (k - 1 as libc::c_int as libc::c_long
                                + mpp * (jdrop - 1 as libc::c_int as libc::c_long))
                                as isize,
                        ) = *con.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                        k += 1;
                    }
                    if trured > zero && trured >= 0.1f64 * prerem {
                        current_block = 1854612888966593411;
                    } else {
                        current_block = 12206176237338959610;
                    }
                }
            } else {
                k = 1 as libc::c_int as libc::c_long;
                while k <= mpp {
                    *datmat.offset(
                        (k - 1 as libc::c_int as libc::c_long
                            + mpp * (jdrop - 1 as libc::c_int as libc::c_long))
                            as isize,
                    ) = *con.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                    k += 1;
                }
                if !(nfvals > np) {
                    if jdrop <= n {
                        if *datmat.offset(
                            (mp - 1 as libc::c_int as libc::c_long
                                + mpp * (np - 1 as libc::c_int as libc::c_long))
                                as isize,
                        ) <= f
                        {
                            *x.offset((jdrop - 1 as libc::c_int as libc::c_long) as isize) = *sim
                                .offset(
                                    (jdrop - 1 as libc::c_int as libc::c_long
                                        + n * (np - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                );
                        } else {
                            *sim.offset(
                                (jdrop - 1 as libc::c_int as libc::c_long
                                    + n * (np - 1 as libc::c_int as libc::c_long))
                                    as isize,
                            ) = *x.offset((jdrop - 1 as libc::c_int as libc::c_long) as isize);
                            k = 1 as libc::c_int as libc::c_long;
                            while k <= mpp {
                                *datmat.offset(
                                    (k - 1 as libc::c_int as libc::c_long
                                        + mpp * (jdrop - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) = *datmat.offset(
                                    (k - 1 as libc::c_int as libc::c_long
                                        + mpp * (np - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                );
                                *datmat.offset(
                                    (k - 1 as libc::c_int as libc::c_long
                                        + mpp * (np - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) = *con.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                                k += 1;
                            }
                            k = 1 as libc::c_int as libc::c_long;
                            while k <= jdrop {
                                temp = zero;
                                *sim.offset(
                                    (jdrop - 1 as libc::c_int as libc::c_long
                                        + n * (k - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) = -rho;
                                i = k;
                                while i <= jdrop {
                                    temp -= *simi.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (k - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                    i += 1;
                                }
                                *simi.offset(
                                    (jdrop - 1 as libc::c_int as libc::c_long
                                        + n * (k - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) = temp;
                                k += 1;
                            }
                        }
                    }
                    if nfvals <= n {
                        jdrop = nfvals;
                        *x.offset((jdrop - 1 as libc::c_int as libc::c_long) as isize) += rho;
                        continue;
                    }
                }
                ibrnch = 1 as libc::c_int as libc::c_long;
                current_block = 1854612888966593411;
            }
            'c_3180: loop {
                match current_block {
                    12206176237338959610 => {
                        if iflag == 0 as libc::c_int as libc::c_long {
                            ibrnch = 0 as libc::c_int as libc::c_long;
                            current_block = 1854612888966593411;
                        } else if rho > rhoend {
                            rho = 0.5f64 * rho;
                            if rho <= 1.5f64 * rhoend {
                                rho = rhoend;
                            }
                            if parmu > zero {
                                let mut cmin: libc::c_double = 0.;
                                let mut cmax: libc::c_double = 0.;
                                let mut denom: libc::c_double = 0.;
                                denom = zero;
                                k = 1 as libc::c_int as libc::c_long;
                                while k <= mp {
                                    cmin = *datmat.offset(
                                        (k - 1 as libc::c_int as libc::c_long
                                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                    cmax = cmin;
                                    i = 1 as libc::c_int as libc::c_long;
                                    while i <= n {
                                        temp = *datmat.offset(
                                            (k - 1 as libc::c_int as libc::c_long
                                                + mpp * (i - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                        cmin = if cmin <= temp { cmin } else { temp };
                                        cmax = if cmax >= temp { cmax } else { temp };
                                        i += 1;
                                    }
                                    if k <= m && cmin < 0.5f64 * cmax {
                                        temp = (if cmax >= zero { cmax } else { zero }) - cmin;
                                        if denom <= zero {
                                            denom = temp;
                                        } else {
                                            denom = if denom <= temp { denom } else { temp };
                                        }
                                    }
                                    k += 1;
                                }
                                if denom == zero {
                                    parmu = zero;
                                } else if cmax - cmin < parmu * denom {
                                    parmu = (cmax - cmin) / denom;
                                }
                            }
                            if iprint >= 2 as libc::c_int as libc::c_long {
                                println!(
                                    "\n   Reduction in RHO to {}  and PARMU = {}\n",
                                    rho, parmu,
                                );
                                if iprint == 2 as libc::c_int as libc::c_long {
                                    print_calcfc(
                                        n,
                                        nfvals,
                                        *datmat.offset(
                                            (mp - 1 as libc::c_int as libc::c_long
                                                + mpp * (np - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ),
                                        *datmat.offset(
                                            (mpp - 1 as libc::c_int as libc::c_long
                                                + mpp * (np - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ),
                                        &mut *sim.offset(
                                            (n * (np - 1 as libc::c_int as libc::c_long)) as isize,
                                        )
                                            as *mut libc::c_double
                                            as *const libc::c_double,
                                    );
                                }
                            }
                            current_block = 1854612888966593411;
                        } else {
                            if iprint >= 1 as libc::c_int as libc::c_long {
                                println!("\n   Normal return from subroutine COBYLA\n");
                            }
                            status = 0 as libc::c_int;
                            if ifull == 1 as libc::c_int as libc::c_long {
                                current_block = 1854573226567202969;
                                break 'c_3166;
                            } else {
                                current_block = 18034400475116118263;
                                break 'c_3166;
                            }
                        }
                    }
                    _ => {
                        phimin = *datmat.offset(
                            (mp - 1 as libc::c_int as libc::c_long
                                + mpp * (np - 1 as libc::c_int as libc::c_long))
                                as isize,
                        ) + parmu
                            * *datmat.offset(
                                (mpp - 1 as libc::c_int as libc::c_long
                                    + mpp * (np - 1 as libc::c_int as libc::c_long))
                                    as isize,
                            );
                        nbest = np;
                        j = 1 as libc::c_int as libc::c_long;
                        while j <= n {
                            temp = *datmat.offset(
                                (mp - 1 as libc::c_int as libc::c_long
                                    + mpp * (j - 1 as libc::c_int as libc::c_long))
                                    as isize,
                            ) + parmu
                                * *datmat.offset(
                                    (mpp - 1 as libc::c_int as libc::c_long
                                        + mpp * (j - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                );
                            if temp < phimin {
                                nbest = j;
                                phimin = temp;
                            } else if temp == phimin && parmu == zero {
                                if *datmat.offset(
                                    (mpp - 1 as libc::c_int as libc::c_long
                                        + mpp * (j - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) < *datmat.offset(
                                    (mpp - 1 as libc::c_int as libc::c_long
                                        + mpp * (nbest - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) {
                                    nbest = j;
                                }
                            }
                            j += 1;
                        }
                        if nbest <= n {
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= mpp {
                                temp = *datmat.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + mpp * (np - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                );
                                *datmat.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + mpp * (np - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) = *datmat.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + mpp * (nbest - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                );
                                *datmat.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + mpp * (nbest - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) = temp;
                                i += 1;
                            }
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                temp = *sim.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (nbest - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                );
                                *sim.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (nbest - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) = zero;
                                *sim.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (np - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) += temp;
                                tempa = zero;
                                k = 1 as libc::c_int as libc::c_long;
                                while k <= n {
                                    *sim.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (k - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) -= temp;
                                    tempa -= *simi.offset(
                                        (k - 1 as libc::c_int as libc::c_long
                                            + n * (i - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                    k += 1;
                                }
                                *simi.offset(
                                    (nbest - 1 as libc::c_int as libc::c_long
                                        + n * (i - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) = tempa;
                                i += 1;
                            }
                        }
                        error = zero;
                        i = 1 as libc::c_int as libc::c_long;
                        while i <= n {
                            j = 1 as libc::c_int as libc::c_long;
                            while j <= n {
                                temp = if i == j { -one } else { zero };
                                k = 1 as libc::c_int as libc::c_long;
                                while k <= n {
                                    temp += *simi.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (k - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) * *sim.offset(
                                        (k - 1 as libc::c_int as libc::c_long
                                            + n * (j - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                    k += 1;
                                }
                                temp = temp.abs();
                                error = if error >= temp { error } else { temp };
                                j += 1;
                            }
                            i += 1;
                        }
                        if error > 0.1f64 {
                            status = -(1 as libc::c_int);
                            if iprint >= 1 as libc::c_int as libc::c_long {
                                println!(
                                    "Return from subroutine COBYLA because {:?}.\n",
                                    cobyla_reason(status),
                                );
                            }
                            current_block = 18034400475116118263;
                            break 'c_3166;
                        } else {
                            k = 1 as libc::c_int as libc::c_long;
                            while k <= mp {
                                *con.offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                    -*datmat.offset(
                                        (k - 1 as libc::c_int as libc::c_long
                                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                j = 1 as libc::c_int as libc::c_long;
                                while j <= n {
                                    *w.offset((j - 1 as libc::c_int as libc::c_long) as isize) =
                                        *datmat.offset(
                                            (k - 1 as libc::c_int as libc::c_long
                                                + mpp * (j - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) + *con.offset(
                                            (k - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                    j += 1;
                                }
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    temp = zero;
                                    j = 1 as libc::c_int as libc::c_long;
                                    while j <= n {
                                        temp += *w.offset(
                                            (j - 1 as libc::c_int as libc::c_long) as isize,
                                        ) * *simi.offset(
                                            (j - 1 as libc::c_int as libc::c_long
                                                + n * (i - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                        j += 1;
                                    }
                                    *a.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (k - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = if k == mp { -temp } else { temp };
                                    i += 1;
                                }
                                k += 1;
                            }
                            iflag = 1 as libc::c_int as libc::c_long;
                            parsig = alpha * rho;
                            pareta = beta * rho;
                            j = 1 as libc::c_int as libc::c_long;
                            while j <= n {
                                let mut wsig: libc::c_double = zero;
                                let mut weta: libc::c_double = zero;
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    wsig += *simi.offset(
                                        (j - 1 as libc::c_int as libc::c_long
                                            + n * (i - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) * *simi.offset(
                                        (j - 1 as libc::c_int as libc::c_long
                                            + n * (i - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                    weta += *sim.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (j - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) * *sim.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (j - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                    i += 1;
                                }
                                *vsig.offset((j - 1 as libc::c_int as libc::c_long) as isize) =
                                    one / wsig.sqrt();
                                *veta.offset((j - 1 as libc::c_int as libc::c_long) as isize) =
                                    weta.sqrt();
                                if *vsig.offset((j - 1 as libc::c_int as libc::c_long) as isize)
                                    < parsig
                                    || *veta.offset((j - 1 as libc::c_int as libc::c_long) as isize)
                                        > pareta
                                {
                                    iflag = 0 as libc::c_int as libc::c_long;
                                }
                                j += 1;
                            }
                            if ibrnch == 1 as libc::c_int as libc::c_long
                                || iflag == 1 as libc::c_int as libc::c_long
                            {
                                let mut z: *mut libc::c_double = w;
                                let mut zdota: *mut libc::c_double = z.offset((n * n) as isize);
                                let mut vmc: *mut libc::c_double = zdota.offset(n as isize);
                                let mut sdirn: *mut libc::c_double = vmc.offset(mp as isize);
                                let mut dxnew: *mut libc::c_double = sdirn.offset(n as isize);
                                let mut vmd: *mut libc::c_double = dxnew.offset(n as isize);
                                trstlp(
                                    n,
                                    m,
                                    a as *const libc::c_double,
                                    con as *const libc::c_double,
                                    rho,
                                    dx,
                                    &mut ifull,
                                    iact,
                                    z,
                                    zdota,
                                    vmc,
                                    sdirn,
                                    dxnew,
                                    vmd,
                                );
                                if ifull == 0 as libc::c_int as libc::c_long {
                                    temp = zero;
                                    i = 1 as libc::c_int as libc::c_long;
                                    while i <= n {
                                        temp += *dx.offset(
                                            (i - 1 as libc::c_int as libc::c_long) as isize,
                                        ) * *dx.offset(
                                            (i - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                        i += 1;
                                    }
                                    if temp < 0.25f64 * rho * rho {
                                        ibrnch = 1 as libc::c_int as libc::c_long;
                                        current_block = 12206176237338959610;
                                        continue;
                                    }
                                }
                                resnew = zero;
                                *con.offset((mp - 1 as libc::c_int as libc::c_long) as isize) =
                                    zero;
                                sum = zero;
                                k = 1 as libc::c_int as libc::c_long;
                                while k <= mp {
                                    sum = *con
                                        .offset((k - 1 as libc::c_int as libc::c_long) as isize);
                                    i = 1 as libc::c_int as libc::c_long;
                                    while i <= n {
                                        sum -= *a.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (k - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) * *dx.offset(
                                            (i - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                        i += 1;
                                    }
                                    if k < mp {
                                        resnew = if resnew >= sum { resnew } else { sum };
                                    }
                                    k += 1;
                                }
                                prerec = *datmat.offset(
                                    (mpp - 1 as libc::c_int as libc::c_long
                                        + mpp * (np - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) - resnew;
                                barmu = if prerec > zero { sum / prerec } else { zero };
                                if !(parmu < 1.5f64 * barmu) {
                                    break;
                                }
                                parmu = 2.0f64 * barmu;
                                if iprint >= 2 as libc::c_int as libc::c_long {
                                    println!("\n   Increase in PARMU to {}\n", parmu);
                                }
                                phi = *datmat.offset(
                                    (mp - 1 as libc::c_int as libc::c_long
                                        + mpp * (np - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) + parmu
                                    * *datmat.offset(
                                        (mpp - 1 as libc::c_int as libc::c_long
                                            + mpp * (np - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                j = 1 as libc::c_int as libc::c_long;
                                loop {
                                    if !(j <= n) {
                                        break 'c_3180;
                                    }
                                    temp = *datmat.offset(
                                        (mp - 1 as libc::c_int as libc::c_long
                                            + mpp * (j - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) + parmu
                                        * *datmat.offset(
                                            (mpp - 1 as libc::c_int as libc::c_long
                                                + mpp * (j - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                    if temp < phi {
                                        current_block = 1854612888966593411;
                                        break;
                                    }
                                    if temp == phi && parmu == zero {
                                        if *datmat.offset(
                                            (mpp - 1 as libc::c_int as libc::c_long
                                                + mpp * (j - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) < *datmat.offset(
                                            (mpp - 1 as libc::c_int as libc::c_long
                                                + mpp * (np - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) {
                                            current_block = 1854612888966593411;
                                            break;
                                        }
                                    }
                                    j += 1;
                                }
                            } else {
                                jdrop = 0 as libc::c_int as libc::c_long;
                                temp = pareta;
                                j = 1 as libc::c_int as libc::c_long;
                                while j <= n {
                                    if *veta.offset((j - 1 as libc::c_int as libc::c_long) as isize)
                                        > temp
                                    {
                                        jdrop = j;
                                        temp = *veta.offset(
                                            (j - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                    }
                                    j += 1;
                                }
                                if jdrop == 0 as libc::c_int as libc::c_long {
                                    j = 1 as libc::c_int as libc::c_long;
                                    while j <= n {
                                        if *vsig
                                            .offset((j - 1 as libc::c_int as libc::c_long) as isize)
                                            < temp
                                        {
                                            jdrop = j;
                                            temp = *vsig.offset(
                                                (j - 1 as libc::c_int as libc::c_long) as isize,
                                            );
                                        }
                                        j += 1;
                                    }
                                }
                                temp = gamma
                                    * rho
                                    * *vsig.offset(
                                        (jdrop - 1 as libc::c_int as libc::c_long) as isize,
                                    );
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize) =
                                        temp * *simi.offset(
                                            (jdrop - 1 as libc::c_int as libc::c_long
                                                + n * (i - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                    i += 1;
                                }
                                cvmaxp = zero;
                                cvmaxm = zero;
                                sum = zero;
                                k = 1 as libc::c_int as libc::c_long;
                                while k <= mp {
                                    sum = zero;
                                    i = 1 as libc::c_int as libc::c_long;
                                    while i <= n {
                                        sum = sum
                                            + *a.offset(
                                                (i - 1 as libc::c_int as libc::c_long
                                                    + n * (k - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            ) * *dx.offset(
                                                (i - 1 as libc::c_int as libc::c_long) as isize,
                                            );
                                        i += 1;
                                    }
                                    if k < mp {
                                        temp = *datmat.offset(
                                            (k - 1 as libc::c_int as libc::c_long
                                                + mpp * (np - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                        cvmaxp = if cvmaxp >= -sum - temp {
                                            cvmaxp
                                        } else {
                                            -sum - temp
                                        };
                                        cvmaxm = if cvmaxm >= sum - temp {
                                            cvmaxm
                                        } else {
                                            sum - temp
                                        };
                                    }
                                    k += 1;
                                }
                                if parmu * (cvmaxp - cvmaxm) > sum + sum {
                                    dxsign = -one;
                                } else {
                                    dxsign = one;
                                }
                                temp = zero;
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize) *=
                                        dxsign;
                                    *sim.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (jdrop - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                                    temp += *simi.offset(
                                        (jdrop - 1 as libc::c_int as libc::c_long
                                            + n * (i - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) * *dx
                                        .offset((i - 1 as libc::c_int as libc::c_long) as isize);
                                    i += 1;
                                }
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    *simi.offset(
                                        (jdrop - 1 as libc::c_int as libc::c_long
                                            + n * (i - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) /= temp;
                                    i += 1;
                                }
                                j = 1 as libc::c_int as libc::c_long;
                                while j <= n {
                                    if j != jdrop {
                                        temp = zero;
                                        i = 1 as libc::c_int as libc::c_long;
                                        while i <= n {
                                            temp += *simi.offset(
                                                (j - 1 as libc::c_int as libc::c_long
                                                    + n * (i - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            ) * *dx.offset(
                                                (i - 1 as libc::c_int as libc::c_long) as isize,
                                            );
                                            i += 1;
                                        }
                                        i = 1 as libc::c_int as libc::c_long;
                                        while i <= n {
                                            *simi.offset(
                                                (j - 1 as libc::c_int as libc::c_long
                                                    + n * (i - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            ) -= temp
                                                * *simi.offset(
                                                    (jdrop - 1 as libc::c_int as libc::c_long
                                                        + n * (i - 1 as libc::c_int
                                                            as libc::c_long))
                                                        as isize,
                                                );
                                            i += 1;
                                        }
                                    }
                                    *x.offset((j - 1 as libc::c_int as libc::c_long) as isize) =
                                        *sim.offset(
                                            (j - 1 as libc::c_int as libc::c_long
                                                + n * (np - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) + *dx.offset(
                                            (j - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                    j += 1;
                                }
                                continue 'c_3166;
                            }
                        }
                    }
                }
            }
            prerem = parmu * prerec - sum;
            i = 1 as libc::c_int as libc::c_long;
            while i <= n {
                *x.offset((i - 1 as libc::c_int as libc::c_long) as isize) = *sim.offset(
                    (i - 1 as libc::c_int as libc::c_long
                        + n * (np - 1 as libc::c_int as libc::c_long)) as isize,
                ) + *dx
                    .offset((i - 1 as libc::c_int as libc::c_long) as isize);
                i += 1;
            }
            ibrnch = 1 as libc::c_int as libc::c_long;
        }
    }
    match current_block {
        18034400475116118263 => {
            i = 1 as libc::c_int as libc::c_long;
            while i <= n {
                *x.offset((i - 1 as libc::c_int as libc::c_long) as isize) = *sim.offset(
                    (i - 1 as libc::c_int as libc::c_long
                        + n * (np - 1 as libc::c_int as libc::c_long)) as isize,
                );
                i += 1;
            }
            f = *datmat.offset(
                (mp - 1 as libc::c_int as libc::c_long
                    + mpp * (np - 1 as libc::c_int as libc::c_long)) as isize,
            );
            resmax = *datmat.offset(
                (mpp - 1 as libc::c_int as libc::c_long
                    + mpp * (np - 1 as libc::c_int as libc::c_long)) as isize,
            );
        }
        _ => {}
    }
    if iprint >= 1 as libc::c_int as libc::c_long {
        print_calcfc(
            n,
            nfvals,
            f,
            resmax,
            &mut *x.offset(0) as *mut libc::c_double as *const libc::c_double,
        );
    }
    *_maxfun = nfvals;
    return status;
}

#[allow(clippy::too_many_arguments)]
unsafe fn trstlp(
    n: libc::c_long,
    m: libc::c_long,
    mut a: *const libc::c_double,
    mut b: *const libc::c_double,
    rho: libc::c_double,
    mut dx: *mut libc::c_double,
    mut ifull: *mut libc::c_long,
    mut iact: *mut libc::c_long,
    mut z: *mut libc::c_double,
    mut zdota: *mut libc::c_double,
    mut vmultc: *mut libc::c_double,
    mut sdirn: *mut libc::c_double,
    mut dxnew: *mut libc::c_double,
    mut vmultd: *mut libc::c_double,
) {
    let mut current_block: u64;
    let zero: libc::c_double = 0.0f64;
    let one: libc::c_double = 1.0f64;
    let tiny: libc::c_double = 1.0e-6f64;
    let Op1: libc::c_double = 0.1f64;
    let Op2: libc::c_double = 0.2f64;
    let mut acca: libc::c_double = 0.;
    let mut accb: libc::c_double = 0.;
    let mut alpha: libc::c_double = 0.;
    let mut beta: libc::c_double = 0.;
    let mut dd: libc::c_double = 0.;
    let mut optnew: libc::c_double = 0.;
    let mut optold: libc::c_double = 0.;
    let mut ratio: libc::c_double = 0.;
    let mut resmax: libc::c_double = 0.;
    let mut resold: libc::c_double = 0.;
    let mut sd: libc::c_double = 0.;
    let mut sp: libc::c_double = 0.;
    let mut spabs: libc::c_double = 0.;
    let mut ss: libc::c_double = 0.;
    let mut step: libc::c_double = 0.;
    let mut stpful: libc::c_double = 0.;
    let mut sum: libc::c_double = 0.;
    let mut sumabs: libc::c_double = 0.;
    let mut temp: libc::c_double = 0.;
    let mut tempa: libc::c_double = 0.;
    let mut tempb: libc::c_double = 0.;
    let mut tot: libc::c_double = 0.;
    let mut vsave: libc::c_double = 0.;
    let mut zdotv: libc::c_double = 0.;
    let mut zdotw: libc::c_double = 0.;
    let mut zdvabs: libc::c_double = 0.;
    let mut zdwabs: libc::c_double = 0.;
    let mut i: libc::c_long = 0;
    let mut icon: libc::c_long = 0;
    let mut icount: libc::c_long = 0;
    let mut isave: libc::c_long = 0;
    let mut j: libc::c_long = 0;
    let mut k: libc::c_long = 0;
    let mut kk: libc::c_long = 0;
    let mut kl: libc::c_long = 0;
    let mut kp: libc::c_long = 0;
    let mut kw: libc::c_long = 0;
    let mut mcon: libc::c_long = 0;
    let mut nact: libc::c_long = 0;
    let mut nactx: libc::c_long = 0;
    *ifull = 1 as libc::c_int as libc::c_long;
    icon = 0 as libc::c_int as libc::c_long;
    mcon = m;
    nact = 0 as libc::c_int as libc::c_long;
    resmax = zero;
    resold = zero;
    i = 1 as libc::c_int as libc::c_long;
    while i <= n {
        j = 1 as libc::c_int as libc::c_long;
        while j <= n {
            *z.offset(
                (i - 1 as libc::c_int as libc::c_long + n * (j - 1 as libc::c_int as libc::c_long))
                    as isize,
            ) = zero;
            j += 1;
        }
        *z.offset(
            (i - 1 as libc::c_int as libc::c_long + n * (i - 1 as libc::c_int as libc::c_long))
                as isize,
        ) = one;
        *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize) = zero;
        i += 1;
    }
    if m >= 1 as libc::c_int as libc::c_long {
        k = 1 as libc::c_int as libc::c_long;
        while k <= m {
            if *b.offset((k - 1 as libc::c_int as libc::c_long) as isize) > resmax {
                resmax = *b.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                icon = k;
            }
            k += 1;
        }
        k = 1 as libc::c_int as libc::c_long;
        while k <= m {
            *iact.offset((k - 1 as libc::c_int as libc::c_long) as isize) = k;
            *vmultc.offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                resmax - *b.offset((k - 1 as libc::c_int as libc::c_long) as isize);
            k += 1;
        }
    }
    if resmax == zero {
        current_block = 16746262731645592041;
    } else {
        i = 1 as libc::c_int as libc::c_long;
        while i <= n {
            *sdirn.offset((i - 1 as libc::c_int as libc::c_long) as isize) = zero;
            i += 1;
        }
        current_block = 10785318058620693244;
    }
    'c_5472: loop {
        match current_block {
            16746262731645592041 => {
                mcon = m + 1 as libc::c_int as libc::c_long;
                icon = mcon;
                *iact.offset((mcon - 1 as libc::c_int as libc::c_long) as isize) = mcon;
                *vmultc.offset((mcon - 1 as libc::c_int as libc::c_long) as isize) = zero;
                current_block = 10785318058620693244;
            }
            _ => {
                optold = zero;
                icount = 0 as libc::c_int as libc::c_long;
                loop {
                    if mcon == m {
                        optnew = resmax;
                    } else {
                        optnew = zero;
                        i = 1 as libc::c_int as libc::c_long;
                        while i <= n {
                            optnew -= *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize)
                                * *a.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (mcon - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                );
                            i += 1;
                        }
                    }
                    if icount == 0 as libc::c_int as libc::c_long || optnew < optold {
                        optold = optnew;
                        nactx = nact;
                        icount = 3 as libc::c_int as libc::c_long;
                    } else if nact > nactx {
                        nactx = nact;
                        icount = 3 as libc::c_int as libc::c_long;
                    } else {
                        icount -= 1;
                        if icount == 0 as libc::c_int as libc::c_long {
                            break;
                        }
                    }
                    if icon <= nact {
                        if icon < nact {
                            isave =
                                *iact.offset((icon - 1 as libc::c_int as libc::c_long) as isize);
                            vsave =
                                *vmultc.offset((icon - 1 as libc::c_int as libc::c_long) as isize);
                            k = icon;
                            loop {
                                kp = k + 1 as libc::c_int as libc::c_long;
                                kk = *iact.offset((kp - 1 as libc::c_int as libc::c_long) as isize);
                                sp = zero;
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    sp += *z.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (k - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) * *a.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (kk - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                    i += 1;
                                }
                                temp = (sp * sp
                                    + *zdota
                                        .offset((kp - 1 as libc::c_int as libc::c_long) as isize)
                                        * *zdota.offset(
                                            (kp - 1 as libc::c_int as libc::c_long) as isize,
                                        ))
                                .sqrt();
                                alpha = *zdota
                                    .offset((kp - 1 as libc::c_int as libc::c_long) as isize)
                                    / temp;
                                beta = sp / temp;
                                *zdota.offset((kp - 1 as libc::c_int as libc::c_long) as isize) =
                                    alpha
                                        * *zdota.offset(
                                            (k - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                *zdota.offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                    temp;
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    temp = alpha
                                        * *z.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (kp - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        )
                                        + beta
                                            * *z.offset(
                                                (i - 1 as libc::c_int as libc::c_long
                                                    + n * (k - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            );
                                    *z.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (kp - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = alpha
                                        * *z.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (k - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        )
                                        - beta
                                            * *z.offset(
                                                (i - 1 as libc::c_int as libc::c_long
                                                    + n * (kp - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            );
                                    *z.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (k - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = temp;
                                    i += 1;
                                }
                                *iact.offset((k - 1 as libc::c_int as libc::c_long) as isize) = kk;
                                *vmultc.offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                    *vmultc
                                        .offset((kp - 1 as libc::c_int as libc::c_long) as isize);
                                k = kp;
                                if !(k < nact) {
                                    break;
                                }
                            }
                            *iact.offset((k - 1 as libc::c_int as libc::c_long) as isize) = isave;
                            *vmultc.offset((k - 1 as libc::c_int as libc::c_long) as isize) = vsave;
                        }
                        nact -= 1;
                        if mcon > m {
                            current_block = 14831642685178214089;
                        } else {
                            temp = zero;
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                temp += *sdirn
                                    .offset((i - 1 as libc::c_int as libc::c_long) as isize)
                                    * *z.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (nact + 1 as libc::c_int as libc::c_long
                                                - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                i += 1;
                            }
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                *sdirn.offset((i - 1 as libc::c_int as libc::c_long) as isize) -=
                                    temp * *z.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (nact + 1 as libc::c_int as libc::c_long
                                                - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                i += 1;
                            }
                            current_block = 18040478258510813106;
                        }
                    } else {
                        kk = *iact.offset((icon - 1 as libc::c_int as libc::c_long) as isize);
                        i = 1 as libc::c_int as libc::c_long;
                        while i <= n {
                            *dxnew.offset((i - 1 as libc::c_int as libc::c_long) as isize) = *a
                                .offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (kk - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                );
                            i += 1;
                        }
                        tot = zero;
                        k = n;
                        while k > nact {
                            sp = zero;
                            spabs = zero;
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                temp = *z.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (k - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) * *dxnew
                                    .offset((i - 1 as libc::c_int as libc::c_long) as isize);
                                sp += temp;
                                spabs += temp.abs();
                                i += 1;
                            }
                            acca = spabs + Op1 * sp.abs();
                            accb = spabs + Op2 * sp.abs();
                            if spabs >= acca || acca >= accb {
                                sp = zero;
                            }
                            if tot == zero {
                                tot = sp;
                            } else {
                                kp = k + 1 as libc::c_int as libc::c_long;
                                temp = (sp * sp + tot * tot).sqrt();
                                alpha = sp / temp;
                                beta = tot / temp;
                                tot = temp;
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    temp = alpha
                                        * *z.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (k - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        )
                                        + beta
                                            * *z.offset(
                                                (i - 1 as libc::c_int as libc::c_long
                                                    + n * (kp - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            );
                                    *z.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (kp - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = alpha
                                        * *z.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (kp - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        )
                                        - beta
                                            * *z.offset(
                                                (i - 1 as libc::c_int as libc::c_long
                                                    + n * (k - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            );
                                    *z.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (k - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) = temp;
                                    i += 1;
                                }
                            }
                            k -= 1;
                        }
                        if tot != zero {
                            nact += 1;
                            *zdota.offset((nact - 1 as libc::c_int as libc::c_long) as isize) = tot;
                            *vmultc.offset((icon - 1 as libc::c_int as libc::c_long) as isize) =
                                *vmultc.offset((nact - 1 as libc::c_int as libc::c_long) as isize);
                            *vmultc.offset((nact - 1 as libc::c_int as libc::c_long) as isize) =
                                zero;
                        } else {
                            ratio = -one;
                            k = nact;
                            loop {
                                zdotv = zero;
                                zdvabs = zero;
                                i = 1 as libc::c_int as libc::c_long;
                                while i <= n {
                                    temp = *z.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (k - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    ) * *dxnew
                                        .offset((i - 1 as libc::c_int as libc::c_long) as isize);
                                    zdotv += temp;
                                    zdvabs += temp.abs();
                                    i += 1;
                                }
                                acca = zdvabs + Op1 * zdotv.abs();
                                accb = zdvabs + Op2 * zdotv.abs();
                                if zdvabs < acca && acca < accb {
                                    temp = zdotv
                                        / *zdota.offset(
                                            (k - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                    if temp > zero
                                        && *iact
                                            .offset((k - 1 as libc::c_int as libc::c_long) as isize)
                                            <= m
                                    {
                                        tempa = *vmultc.offset(
                                            (k - 1 as libc::c_int as libc::c_long) as isize,
                                        ) / temp;
                                        if ratio < zero || tempa < ratio {
                                            ratio = tempa;
                                        }
                                    }
                                    if k >= 2 as libc::c_int as libc::c_long {
                                        kw = *iact.offset(
                                            (k - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                        i = 1 as libc::c_int as libc::c_long;
                                        while i <= n {
                                            *dxnew.offset(
                                                (i - 1 as libc::c_int as libc::c_long) as isize,
                                            ) -= temp
                                                * *a.offset(
                                                    (i - 1 as libc::c_int as libc::c_long
                                                        + n * (kw
                                                            - 1 as libc::c_int as libc::c_long))
                                                        as isize,
                                                );
                                            i += 1;
                                        }
                                    }
                                    *vmultd
                                        .offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                        temp;
                                } else {
                                    *vmultd
                                        .offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                        zero;
                                }
                                k -= 1;
                                if !(k > 0 as libc::c_int as libc::c_long) {
                                    break;
                                }
                            }
                            if ratio < zero {
                                break;
                            }
                            k = 1 as libc::c_int as libc::c_long;
                            while k <= nact {
                                tempb = *vmultc
                                    .offset((k - 1 as libc::c_int as libc::c_long) as isize)
                                    - ratio
                                        * *vmultd.offset(
                                            (k - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                *vmultc.offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                    if zero >= tempb { zero } else { tempb };
                                k += 1;
                            }
                            if icon < nact {
                                isave = *iact
                                    .offset((icon - 1 as libc::c_int as libc::c_long) as isize);
                                vsave = *vmultc
                                    .offset((icon - 1 as libc::c_int as libc::c_long) as isize);
                                k = icon;
                                loop {
                                    kp = k + 1 as libc::c_int as libc::c_long;
                                    kw = *iact
                                        .offset((kp - 1 as libc::c_int as libc::c_long) as isize);
                                    sp = zero;
                                    i = 1 as libc::c_int as libc::c_long;
                                    while i <= n {
                                        sp += *z.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (k - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) * *a.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (kw - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                        i += 1;
                                    }
                                    temp = (sp * sp
                                        + *zdota.offset(
                                            (kp - 1 as libc::c_int as libc::c_long) as isize,
                                        ) * *zdota.offset(
                                            (kp - 1 as libc::c_int as libc::c_long) as isize,
                                        ))
                                    .sqrt();
                                    alpha = *zdota
                                        .offset((kp - 1 as libc::c_int as libc::c_long) as isize)
                                        / temp;
                                    beta = sp / temp;
                                    *zdota
                                        .offset((kp - 1 as libc::c_int as libc::c_long) as isize) =
                                        alpha
                                            * *zdota.offset(
                                                (k - 1 as libc::c_int as libc::c_long) as isize,
                                            );
                                    *zdota
                                        .offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                        temp;
                                    i = 1 as libc::c_int as libc::c_long;
                                    while i <= n {
                                        temp = alpha
                                            * *z.offset(
                                                (i - 1 as libc::c_int as libc::c_long
                                                    + n * (kp - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            )
                                            + beta
                                                * *z.offset(
                                                    (i - 1 as libc::c_int as libc::c_long
                                                        + n * (k - 1 as libc::c_int
                                                            as libc::c_long))
                                                        as isize,
                                                );
                                        *z.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (kp - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) = alpha
                                            * *z.offset(
                                                (i - 1 as libc::c_int as libc::c_long
                                                    + n * (k - 1 as libc::c_int as libc::c_long))
                                                    as isize,
                                            )
                                            - beta
                                                * *z.offset(
                                                    (i - 1 as libc::c_int as libc::c_long
                                                        + n * (kp
                                                            - 1 as libc::c_int as libc::c_long))
                                                        as isize,
                                                );
                                        *z.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (k - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        ) = temp;
                                        i += 1;
                                    }
                                    *iact.offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                        kw;
                                    *vmultc
                                        .offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                        *vmultc.offset(
                                            (kp - 1 as libc::c_int as libc::c_long) as isize,
                                        );
                                    k = kp;
                                    if !(k < nact) {
                                        break;
                                    }
                                }
                                *iact.offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                    isave;
                                *vmultc.offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                    vsave;
                            }
                            temp = zero;
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                temp += *z.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (nact - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) * *a.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (kk - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                );
                                i += 1;
                            }
                            if temp == zero {
                                break;
                            }
                            *zdota.offset((nact - 1 as libc::c_int as libc::c_long) as isize) =
                                temp;
                            *vmultc.offset((icon - 1 as libc::c_int as libc::c_long) as isize) =
                                zero;
                            *vmultc.offset((nact - 1 as libc::c_int as libc::c_long) as isize) =
                                ratio;
                        }
                        *iact.offset((icon - 1 as libc::c_int as libc::c_long) as isize) =
                            *iact.offset((nact - 1 as libc::c_int as libc::c_long) as isize);
                        *iact.offset((nact - 1 as libc::c_int as libc::c_long) as isize) = kk;
                        if mcon > m && kk != mcon {
                            k = nact - 1 as libc::c_int as libc::c_long;
                            sp = zero;
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                sp += *z.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (k - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) * *a.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (kk - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                );
                                i += 1;
                            }
                            temp = (sp * sp
                                + *zdota
                                    .offset((nact - 1 as libc::c_int as libc::c_long) as isize)
                                    * *zdota.offset(
                                        (nact - 1 as libc::c_int as libc::c_long) as isize,
                                    ))
                            .sqrt();
                            alpha = *zdota
                                .offset((nact - 1 as libc::c_int as libc::c_long) as isize)
                                / temp;
                            beta = sp / temp;
                            *zdota.offset((nact - 1 as libc::c_int as libc::c_long) as isize) =
                                alpha
                                    * *zdota
                                        .offset((k - 1 as libc::c_int as libc::c_long) as isize);
                            *zdota.offset((k - 1 as libc::c_int as libc::c_long) as isize) = temp;
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                temp = alpha
                                    * *z.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (nact - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    )
                                    + beta
                                        * *z.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (k - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                *z.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (nact - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) = alpha
                                    * *z.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (k - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    )
                                    - beta
                                        * *z.offset(
                                            (i - 1 as libc::c_int as libc::c_long
                                                + n * (nact - 1 as libc::c_int as libc::c_long))
                                                as isize,
                                        );
                                *z.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (k - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) = temp;
                                i += 1;
                            }
                            *iact.offset((nact - 1 as libc::c_int as libc::c_long) as isize) =
                                *iact.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                            *iact.offset((k - 1 as libc::c_int as libc::c_long) as isize) = kk;
                            temp = *vmultc.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                            *vmultc.offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                                *vmultc.offset((nact - 1 as libc::c_int as libc::c_long) as isize);
                            *vmultc.offset((nact - 1 as libc::c_int as libc::c_long) as isize) =
                                temp;
                        }
                        if mcon > m {
                            current_block = 14831642685178214089;
                        } else {
                            kk = *iact.offset((nact - 1 as libc::c_int as libc::c_long) as isize);
                            temp = zero;
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                temp += *sdirn
                                    .offset((i - 1 as libc::c_int as libc::c_long) as isize)
                                    * *a.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (kk - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                i += 1;
                            }
                            temp = (temp - one)
                                / *zdota.offset((nact - 1 as libc::c_int as libc::c_long) as isize);
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                *sdirn.offset((i - 1 as libc::c_int as libc::c_long) as isize) -=
                                    temp * *z.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (nact - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                i += 1;
                            }
                            current_block = 18040478258510813106;
                        }
                    }
                    match current_block {
                        14831642685178214089 => {
                            temp = one
                                / *zdota.offset((nact - 1 as libc::c_int as libc::c_long) as isize);
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                *sdirn.offset((i - 1 as libc::c_int as libc::c_long) as isize) =
                                    temp * *z.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (nact - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                                i += 1;
                            }
                        }
                        _ => {}
                    }
                    dd = rho * rho;
                    sd = zero;
                    ss = zero;
                    i = 1 as libc::c_int as libc::c_long;
                    while i <= n {
                        if (*dx.offset((i - 1 as libc::c_int as libc::c_long) as isize)).abs()
                            >= tiny * rho
                        {
                            dd -= *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize)
                                * *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                        }
                        sd += *sdirn.offset((i - 1 as libc::c_int as libc::c_long) as isize)
                            * *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                        ss += *sdirn.offset((i - 1 as libc::c_int as libc::c_long) as isize)
                            * *sdirn.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                        i += 1;
                    }
                    if dd <= zero {
                        break;
                    }
                    temp = (ss * dd).sqrt();
                    if sd.abs() >= tiny * temp {
                        temp = (ss * dd + sd * sd).sqrt();
                    }
                    stpful = dd / (temp + sd);
                    step = stpful;
                    if mcon == m {
                        acca = step + Op1 * resmax;
                        accb = step + Op2 * resmax;
                        if step >= acca || acca >= accb {
                            current_block = 16746262731645592041;
                            continue 'c_5472;
                        }
                        step = if step <= resmax { step } else { resmax };
                    }
                    i = 1 as libc::c_int as libc::c_long;
                    while i <= n {
                        *dxnew.offset((i - 1 as libc::c_int as libc::c_long) as isize) = *dx
                            .offset((i - 1 as libc::c_int as libc::c_long) as isize)
                            + step * *sdirn.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                        i += 1;
                    }
                    if mcon == m {
                        resold = resmax;
                        resmax = zero;
                        k = 1 as libc::c_int as libc::c_long;
                        while k <= nact {
                            kk = *iact.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                            temp = *b.offset((kk - 1 as libc::c_int as libc::c_long) as isize);
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                temp -= *a.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (kk - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) * *dxnew
                                    .offset((i - 1 as libc::c_int as libc::c_long) as isize);
                                i += 1;
                            }
                            resmax = if resmax >= temp { resmax } else { temp };
                            k += 1;
                        }
                    }
                    k = nact;
                    loop {
                        zdotw = zero;
                        zdwabs = zero;
                        i = 1 as libc::c_int as libc::c_long;
                        while i <= n {
                            temp = *z.offset(
                                (i - 1 as libc::c_int as libc::c_long
                                    + n * (k - 1 as libc::c_int as libc::c_long))
                                    as isize,
                            ) * *dxnew
                                .offset((i - 1 as libc::c_int as libc::c_long) as isize);
                            zdotw += temp;
                            zdwabs += temp.abs();
                            i += 1;
                        }
                        acca = zdwabs + Op1 * zdotw.abs();
                        accb = zdwabs + Op2 * zdotw.abs();
                        if zdwabs >= acca || acca >= accb {
                            zdotw = zero;
                        }
                        *vmultd.offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                            zdotw / *zdota.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                        if k < 2 as libc::c_int as libc::c_long {
                            break;
                        }
                        kk = *iact.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                        i = 1 as libc::c_int as libc::c_long;
                        while i <= n {
                            *dxnew.offset((i - 1 as libc::c_int as libc::c_long) as isize) -=
                                *vmultd.offset((k - 1 as libc::c_int as libc::c_long) as isize)
                                    * *a.offset(
                                        (i - 1 as libc::c_int as libc::c_long
                                            + n * (kk - 1 as libc::c_int as libc::c_long))
                                            as isize,
                                    );
                            i += 1;
                        }
                        k -= 1;
                    }
                    if mcon > m
                        && *vmultd.offset((nact - 1 as libc::c_int as libc::c_long) as isize) < zero
                    {
                        *vmultd.offset((nact - 1 as libc::c_int as libc::c_long) as isize) = zero;
                    }
                    i = 1 as libc::c_int as libc::c_long;
                    while i <= n {
                        *dxnew.offset((i - 1 as libc::c_int as libc::c_long) as isize) = *dx
                            .offset((i - 1 as libc::c_int as libc::c_long) as isize)
                            + step * *sdirn.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                        i += 1;
                    }
                    if mcon > nact {
                        kl = nact + 1 as libc::c_int as libc::c_long;
                        k = kl;
                        while k <= mcon {
                            kk = *iact.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                            sum = resmax
                                - *b.offset((kk - 1 as libc::c_int as libc::c_long) as isize);
                            sumabs = resmax
                                + (*b.offset((kk - 1 as libc::c_int as libc::c_long) as isize))
                                    .abs();
                            i = 1 as libc::c_int as libc::c_long;
                            while i <= n {
                                temp = *a.offset(
                                    (i - 1 as libc::c_int as libc::c_long
                                        + n * (kk - 1 as libc::c_int as libc::c_long))
                                        as isize,
                                ) * *dxnew
                                    .offset((i - 1 as libc::c_int as libc::c_long) as isize);
                                sum += temp;
                                sumabs += temp.abs();
                                i += 1;
                            }
                            acca = sumabs + Op1 * sum.abs();
                            accb = sumabs + Op2 * sum.abs();
                            if sumabs >= acca || acca >= accb {
                                sum = zero;
                            }
                            *vmultd.offset((k - 1 as libc::c_int as libc::c_long) as isize) = sum;
                            k += 1;
                        }
                    }
                    ratio = one;
                    icon = 0 as libc::c_int as libc::c_long;
                    k = 1 as libc::c_int as libc::c_long;
                    while k <= mcon {
                        if *vmultd.offset((k - 1 as libc::c_int as libc::c_long) as isize) < zero {
                            temp = *vmultc.offset((k - 1 as libc::c_int as libc::c_long) as isize)
                                / (*vmultc.offset((k - 1 as libc::c_int as libc::c_long) as isize)
                                    - *vmultd
                                        .offset((k - 1 as libc::c_int as libc::c_long) as isize));
                            if temp < ratio {
                                ratio = temp;
                                icon = k;
                            }
                        }
                        k += 1;
                    }
                    temp = one - ratio;
                    i = 1 as libc::c_int as libc::c_long;
                    while i <= n {
                        *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize) = temp
                            * *dx.offset((i - 1 as libc::c_int as libc::c_long) as isize)
                            + ratio
                                * *dxnew.offset((i - 1 as libc::c_int as libc::c_long) as isize);
                        i += 1;
                    }
                    k = 1 as libc::c_int as libc::c_long;
                    while k <= mcon {
                        tempb = temp
                            * *vmultc.offset((k - 1 as libc::c_int as libc::c_long) as isize)
                            + ratio
                                * *vmultd.offset((k - 1 as libc::c_int as libc::c_long) as isize);
                        *vmultc.offset((k - 1 as libc::c_int as libc::c_long) as isize) =
                            if zero >= tempb { zero } else { tempb };
                        k += 1;
                    }
                    if mcon == m {
                        resmax = resold + ratio * (resmax - resold);
                    }
                    if icon > 0 as libc::c_int as libc::c_long {
                        continue;
                    }
                    if step == stpful {
                        break 'c_5472;
                    } else {
                        current_block = 16746262731645592041;
                        continue 'c_5472;
                    }
                }
                if mcon == m {
                    current_block = 16746262731645592041;
                    continue;
                }
                *ifull = 0 as libc::c_int as libc::c_long;
                break;
            }
        }
    }
}
#[no_mangle]
pub unsafe fn cobyla_reason(mut status: libc::c_int) -> *const libc::c_char {
    match status {
        1 => {
            return b"user requested to compute F(X) and C(X)" as *const u8 as *const libc::c_char;
        }
        0 => return b"algorithm was successful" as *const u8 as *const libc::c_char,
        -1 => {
            return b"rounding errors are becoming damaging" as *const u8 as *const libc::c_char;
        }
        -2 => {
            return b"MAXFUN limit has been reached" as *const u8 as *const libc::c_char;
        }
        -3 => return b"illegal NULL address" as *const u8 as *const libc::c_char,
        -4 => {
            return b"unexpected parameter or corrupted workspace" as *const u8
                as *const libc::c_char;
        }
        _ => return b"unknown status" as *const u8 as *const libc::c_char,
    };
}
unsafe fn print_calcfc(
    mut n: libc::c_long,
    mut nfvals: libc::c_long,
    mut f: libc::c_double,
    mut maxcv: libc::c_double,
    mut x: *const libc::c_double,
) {
    let mut i: libc::c_long = 0;
    println!(
        "\n   NFVALS ={}   F ={}    MAXCV ={}\n   X ={}",
        nfvals as libc::c_int,
        f,
        maxcv,
        *x.offset(0 as libc::c_int as isize),
    );
    i = 1 as libc::c_int as libc::c_long;
    while i < n {
        // let fmt = if i % 5 as libc::c_int as libc::c_long == 0 as libc::c_int as libc::c_long {
        //     "\n{}"
        // } else {
        //     "{}"
        // };
        println!("{}", *x.offset(i as isize));
        i += 1;
    }
    println!("\n");
}
