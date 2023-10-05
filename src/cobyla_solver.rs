use crate::cobyla::{
    cobyla_context_t, cobyla_create, cobyla_delete, cobyla_get_status, cobyla_iterate,
    cobyla_reason, CobylaStatus,
};
use crate::cobyla_state::*;
use std::mem::ManuallyDrop;

use argmin::core::{CostFunction, Problem, Solver, State, TerminationStatus, KV};
use serde::{Deserialize, Serialize};

/// [Argmin Solver](https://www.argmin-rs.org/book/index.html) which implements COBYLA method.
///
/// ```
/// use argmin::core::{CostFunction, Error, Executor};
/// use argmin::core::observers::{ObserverMode, SlogLogger};
/// use cobyla::CobylaSolver;
///
/// struct ParaboloidProblem;
/// impl CostFunction for ParaboloidProblem {
///     type Param = Vec<f64>;
///     type Output = Vec<f64>;
///
///     // Minimize 10*(x0+1)^2 + x1^2 subject to x0 >= 0
///     fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> {
///         Ok(vec![10. * (x[0] + 1.).powf(2.) + x[1].powf(2.), x[0]])
///     }
/// }
///
/// let pb = ParaboloidProblem;
/// let solver = CobylaSolver::new(vec![1., 1.]);
///
/// let res = Executor::new(pb, solver)
///             .configure(|state| state.max_iters(100))
///             .add_observer(SlogLogger::term(), ObserverMode::Always)
///             .run()
///             .unwrap();
///
/// // Wait a second (lets the logger flush everything before printing again)
/// std::thread::sleep(std::time::Duration::from_secs(1));
///
/// println!("Result of COBYLA:\n{}", res);
/// ```
#[derive(Clone, Serialize, Deserialize)]
pub struct CobylaSolver {
    /// Initial guess for x value
    x0: Vec<f64>,
}

impl CobylaSolver {
    pub fn new(x0: Vec<f64>) -> Self {
        CobylaSolver { x0 }
    }
}

impl<O> Solver<O, CobylaState> for CobylaSolver
where
    O: CostFunction<Param = Vec<f64>, Output = Vec<f64>>,
{
    const NAME: &'static str = "COBYLA";

    /// Initializes the algorithm.
    ///
    /// Executed before any iterations are performed and has access to the optimization problem
    /// definition and the internal state of the solver.
    /// Returns an updated `state` and optionally a `KV` which holds key-value pairs used in
    /// [Observers](`argmin::core::observers::Observe`).
    /// The default implementation returns the unaltered `state` and no `KV`.
    #[allow(clippy::useless_conversion)]
    fn init(
        &mut self,
        problem: &mut Problem<O>,
        state: CobylaState,
    ) -> std::result::Result<(CobylaState, Option<KV>), argmin::core::Error> {
        let n = self.x0.len() as i32;
        let fx0 = problem.cost(&self.x0)?;
        let m = (fx0.len() - 1) as i32;
        let rhobeg = state.rhobeg();
        let rhoend = state.get_rhoend();
        let iprint = state.get_iprint();
        let maxfun = state.get_maxfun();
        let mut initial_state = state;
        let ptr = unsafe {
            cobyla_create(
                n.into(),
                m.into(),
                rhobeg,
                rhoend,
                iprint.into(),
                maxfun.into(),
            )
        };
        initial_state.cobyla_context = Some(ManuallyDrop::new(ptr));

        let initial_state = initial_state.param(self.x0.clone()).cost(fx0);
        Ok((initial_state, None))
    }

    /// Computes a single iteration of the algorithm and has access to the optimization problem
    /// definition and the internal state of the solver.
    /// Returns an updated `state` and optionally a `KV` which holds key-value pairs used in
    /// [Observers](`argmin::core::observers::Observe`).
    fn next_iter(
        &mut self,
        problem: &mut Problem<O>,
        state: CobylaState,
    ) -> std::result::Result<(CobylaState, Option<KV>), argmin::core::Error> {
        let mut x = state.get_param().unwrap().clone();
        if let Some(ctx) = state.cobyla_context.as_ref() {
            let cost = problem.cost(&x)?;
            let f = cost[0];
            let mut c = Box::new(cost[1..].to_vec());

            let _status = unsafe {
                cobyla_iterate(
                    **ctx as *mut cobyla_context_t,
                    f,
                    x.as_mut_ptr(),
                    c.as_mut_ptr(),
                )
            };
            let fx = problem.cost(&x)?;
            let state = state.param(x).cost(fx);
            return Ok((state, None));
        }

        Ok((state, None))
    }

    /// Used to implement stopping criteria, in particular criteria which are not covered by
    /// ([`terminate_internal`](`Solver::terminate_internal`).
    ///
    /// This method has access to the internal state and returns an `TerminationReason`.
    fn terminate(&mut self, state: &CobylaState) -> TerminationStatus {
        if let Some(ctx) = state.cobyla_context.as_ref() {
            let status = unsafe {
                let ctx_ptr = **ctx;
                cobyla_get_status(ctx_ptr)
            };
            if status == CobylaStatus::COBYLA_ITERATE as i32 {
                return TerminationStatus::NotTerminated;
            } else {
                let cstr = unsafe { std::ffi::CStr::from_ptr(cobyla_reason(status)) };
                let reason = cstr.to_str().unwrap().to_string();
                unsafe { cobyla_delete(**ctx as *mut cobyla_context_t) }
                if reason == "algorithm was successful" {
                    return TerminationStatus::Terminated(
                        argmin::core::TerminationReason::SolverConverged,
                    );
                }
                return TerminationStatus::Terminated(argmin::core::TerminationReason::SolverExit(
                    reason,
                ));
            }
        }
        TerminationStatus::Terminated(argmin::core::TerminationReason::SolverExit(
            "Unknown".to_string(),
        ))
    }
}
