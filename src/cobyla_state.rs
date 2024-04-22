use crate::cobyla::cobyla_context_t;
/// Implementation of `argmin::IterState` for Cobyla optimizer
use argmin::core::{Problem, State, TerminationReason, TerminationStatus};
#[cfg(feature = "serde1")]
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::mem::ManuallyDrop;

/// Maintains the state from iteration to iteration of the [crate::CobylaSolver].
///
/// This struct is passed from one iteration of an algorithm to the next.
///
/// Keeps track of
///
/// * parameter vector of current and previous iteration
/// * best parameter vector of current and previous iteration
/// * cost function value (objective and constraint functions values) of current and previous iteration
/// * current and previous best cost function value
/// * target cost function value
/// * current iteration number
/// * iteration number where the last best parameter vector was found
/// * maximum number of iterations that will be executed
/// * problem function evaluation counts
/// * elapsed time
/// * termination status
/// * COBYLA specific parameters: rhobeg, rhoend, iprint, maxfun
///
#[derive(Clone, Debug, Default)]
#[cfg_attr(feature = "serde1", derive(Serialize, Deserialize))]
pub struct CobylaState {
    /// Current parameter vector
    pub param: Option<Vec<f64>>,
    /// Previous parameter vector
    pub prev_param: Option<Vec<f64>>,
    /// Current best parameter vector
    pub best_param: Option<Vec<f64>>,
    /// Previous best parameter vector
    pub prev_best_param: Option<Vec<f64>>,

    /// Current cost function value
    pub cost: Option<Vec<f64>>,
    /// Previous cost function value
    pub prev_cost: Option<Vec<f64>>,
    /// Current best cost function value
    pub best_cost: Option<Vec<f64>>,
    /// Previous best cost function value
    pub prev_best_cost: Option<Vec<f64>>,
    /// Target cost function value
    pub target_cost: f64,

    /// Current iteration
    pub iter: u64,
    /// Iteration number of last best cost
    pub last_best_iter: u64,
    /// Maximum number of iterations
    pub max_iters: u64,
    /// Evaluation counts
    pub counts: HashMap<String, u64>,
    /// Time required so far
    pub time: Option<instant::Duration>,
    /// Status of optimization execution
    pub termination_status: TerminationStatus,

    /// Rho start value
    pub rhobeg: f64,
    /// Rho end value
    pub rhoend: f64,
    /// Control of traces
    pub iprint: i32,
    /// Cost function calls budget
    pub maxfun: i32,

    #[cfg_attr(feature = "serde1", serde(skip))]
    pub cobyla_context: Option<ManuallyDrop<*mut cobyla_context_t>>,
}

impl CobylaState
where
    Self: State<Float = f64>,
{
    /// Set parameter vector. This shifts the stored parameter vector to the previous parameter
    /// vector.
    ///
    /// # Example
    ///
    /// ```
    /// # use argmin::core::{IterState, State};
    /// # use cobyla::CobylaState;
    /// # let state = CobylaState::new();
    /// # let param_old = vec![1.0f64, 2.0f64];
    /// # let state = state.param(param_old);
    /// # assert!(state.prev_param.is_none());
    /// # assert_eq!(state.param.as_ref().unwrap()[0].to_ne_bytes(), 1.0f64.to_ne_bytes());
    /// # assert_eq!(state.param.as_ref().unwrap()[1].to_ne_bytes(), 2.0f64.to_ne_bytes());
    /// # let param = vec![0.0f64, 3.0f64];
    /// let state = state.param(param);
    /// # assert_eq!(state.prev_param.as_ref().unwrap()[0].to_ne_bytes(), 1.0f64.to_ne_bytes());
    /// # assert_eq!(state.prev_param.as_ref().unwrap()[1].to_ne_bytes(), 2.0f64.to_ne_bytes());
    /// # assert_eq!(state.param.as_ref().unwrap()[0].to_ne_bytes(), 0.0f64.to_ne_bytes());
    /// # assert_eq!(state.param.as_ref().unwrap()[1].to_ne_bytes(), 3.0f64.to_ne_bytes());
    /// ```
    #[must_use]
    pub fn param(mut self, param: Vec<f64>) -> Self {
        std::mem::swap(&mut self.prev_param, &mut self.param);
        self.param = Some(param);
        self
    }

    /// Set target cost.
    ///
    /// When this cost is reached, the algorithm will stop. The default is
    /// `Self::Float::NEG_INFINITY`.
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let state: CobylaState = CobylaState::new();
    /// # assert_eq!(state.target_cost.to_ne_bytes(), f64::NEG_INFINITY.to_ne_bytes());
    /// let state = state.target_cost(0.0);
    /// # assert_eq!(state.target_cost.to_ne_bytes(), 0.0f64.to_ne_bytes());
    /// ```
    #[must_use]
    pub fn target_cost(mut self, target_cost: f64) -> Self {
        self.target_cost = target_cost;
        self
    }

    /// Set maximum number of iterations
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let state: CobylaState = CobylaState::new();
    /// # assert_eq!(state.max_iters, std::u64::MAX);
    /// let state = state.max_iters(1000);
    /// # assert_eq!(state.max_iters, 1000);
    /// ```
    #[must_use]
    pub fn max_iters(mut self, iters: u64) -> Self {
        self.max_iters = iters;
        self
    }
    /// Alias for max_iters using historic cobyla terminology
    #[must_use]
    pub fn maxfun(mut self, maxfun: u64) -> Self {
        self.max_iters = maxfun;
        self
    }

    /// Set maximum number of iterations
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let state: CobylaState = CobylaState::new();
    /// # assert_eq!(state.iprint, 1);
    /// let state = state.iprint(0);
    /// # assert_eq!(state.iprint, 0);
    /// ```
    #[must_use]
    pub fn iprint(mut self, iprint: i32) -> Self {
        self.iprint = iprint;
        self
    }

    /// Set the current cost function value. This shifts the stored cost function value to the
    /// previous cost function value.
    ///
    /// # Example
    ///
    /// ```
    /// # use argmin::core::State;
    /// # use cobyla::CobylaState;
    /// # let state: CobylaState = CobylaState::new();
    /// # let cost_old = 1.0f64;
    /// # let state = state.cost(vec![cost_old]);
    /// # assert!(state.prev_cost.is_none());
    /// # assert_eq!(state.cost.as_ref().unwrap()[0].to_ne_bytes(), 1.0f64.to_ne_bytes());
    /// # let cost = 0.0f64;
    /// let state = state.cost(vec![cost]);
    /// # assert_eq!(state.prev_cost.as_ref().unwrap()[0].to_ne_bytes(), 1.0f64.to_ne_bytes());
    /// # assert_eq!(state.cost.as_ref().unwrap()[0].to_ne_bytes(), 0.0f64.to_ne_bytes());
    /// ```
    #[must_use]
    pub fn cost(mut self, cost: Vec<f64>) -> Self {
        std::mem::swap(&mut self.prev_cost, &mut self.cost);
        self.cost = Some(cost);
        self
    }

    /// Returns current cost (ie objective) function and constraint values.
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # state.cost = Some(vec![12.0, 0.1]);
    /// let cost = state.get_full_cost();
    /// # assert_eq!(cost.unwrap()[0].to_ne_bytes(), 12.0f64.to_ne_bytes());
    /// # assert_eq!(cost.unwrap()[1].to_ne_bytes(), 0.1f64.to_ne_bytes());
    /// ```
    pub fn get_full_cost(&self) -> Option<&Vec<f64>> {
        self.cost.as_ref()
    }

    /// Returns current cost (ie objective) function and constraint values.
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # state.best_cost = Some(vec![12.0, 0.1]);
    /// let cost = state.get_full_best_cost();
    /// # assert_eq!(cost.unwrap()[0].to_ne_bytes(), 12.0f64.to_ne_bytes());
    /// # assert_eq!(cost.unwrap()[1].to_ne_bytes(), 0.1f64.to_ne_bytes());
    /// ```
    pub fn get_full_best_cost(&self) -> Option<&Vec<f64>> {
        self.best_cost.as_ref()
    }

    /// Returns the rho start value
    pub fn rhobeg(&self) -> f64 {
        self.rhobeg
    }

    /// Returns the rho end value
    pub fn get_rhoend(&self) -> f64 {
        self.rhoend
    }

    /// Returns the level of printing
    pub fn get_iprint(&self) -> i32 {
        self.iprint
    }

    /// Returns cost function calls budget
    pub fn get_maxfun(&self) -> i32 {
        self.max_iters as i32
    }
}

impl State for CobylaState {
    /// Type of parameter vector
    type Param = Vec<f64>;
    /// Floating point precision
    type Float = f64;

    /// Create new `CobylaState` instance
    ///
    /// # Example
    ///
    /// ```
    /// # extern crate instant;
    /// # use instant;
    /// # use std::collections::HashMap;
    /// # use argmin::core::{State, TerminationStatus};
    /// use cobyla::CobylaState;
    /// let state: CobylaState = CobylaState::new();
    ///
    /// # assert!(state.param.is_none());
    /// # assert!(state.prev_param.is_none());
    /// # assert!(state.best_param.is_none());
    /// # assert!(state.prev_best_param.is_none());
    /// # assert!(state.cost.is_none());
    /// # assert!(state.prev_cost.is_none());
    /// # assert!(state.best_cost.is_none());
    /// # assert!(state.prev_best_cost.is_none());
    /// # assert_eq!(state.target_cost, f64::NEG_INFINITY);
    /// # assert_eq!(state.iter, 0);
    /// # assert_eq!(state.last_best_iter, 0);
    /// # assert_eq!(state.max_iters, std::u64::MAX);
    /// # assert_eq!(state.counts, HashMap::new());
    /// # assert_eq!(state.time.unwrap(), instant::Duration::new(0, 0));
    /// # assert_eq!(state.termination_status, TerminationStatus::NotTerminated);
    /// ```
    fn new() -> Self {
        CobylaState {
            param: None,
            prev_param: None,
            best_param: None,
            prev_best_param: None,

            cost: None,
            prev_cost: None,
            best_cost: None,
            prev_best_cost: None,
            target_cost: f64::NEG_INFINITY,

            iter: 0,
            last_best_iter: 0,
            max_iters: std::u64::MAX,
            counts: HashMap::new(),
            time: Some(instant::Duration::new(0, 0)),
            termination_status: TerminationStatus::NotTerminated,

            rhobeg: 0.5,
            rhoend: 1e-4,
            iprint: 1,
            maxfun: 2000,

            cobyla_context: None,
        }
    }

    /// Checks if the current parameter vector is better than the previous best parameter value. If
    /// a new best parameter vector was found, the state is updated accordingly.
    ///
    /// # Example
    ///
    /// ```
    /// # use argmin::core::{State, ArgminFloat};
    /// # use cobyla::CobylaState;
    ///
    /// let mut state: CobylaState = CobylaState::new();
    ///
    /// // Simulating a new parameter vector
    /// state.param = Some(vec![2.0f64]);
    /// state.cost = Some(vec![5.0]);
    ///
    /// // Calling update
    /// state.update();
    ///
    /// // Check if update was successful
    /// assert_eq!(state.best_param.as_ref().unwrap()[0], 2.0f64);
    /// assert_eq!(state.best_cost.as_ref().unwrap()[0], 5.0);
    /// assert!(state.is_best());
    /// ```
    fn update(&mut self) {
        if let Some(cost) = self.cost.as_ref() {
            if let Some(param) = self.param.as_ref().cloned() {
                std::mem::swap(&mut self.prev_best_param, &mut self.best_param);
                self.best_param = Some(param);
            }
            std::mem::swap(&mut self.prev_best_cost, &mut self.best_cost);
            self.best_cost = Some(cost.clone());
            self.last_best_iter = self.iter;
        }
    }

    /// Returns a reference to the current parameter vector
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # assert!(state.param.is_none());
    /// # state.param = Some(vec![1.0, 2.0]);
    /// # assert_eq!(state.param.as_ref().unwrap()[0].to_ne_bytes(), 1.0f64.to_ne_bytes());
    /// # assert_eq!(state.param.as_ref().unwrap()[1].to_ne_bytes(), 2.0f64.to_ne_bytes());
    /// let param = state.get_param();  // Option<&P>
    /// # assert_eq!(param.as_ref().unwrap()[0].to_ne_bytes(), 1.0f64.to_ne_bytes());
    /// # assert_eq!(param.as_ref().unwrap()[1].to_ne_bytes(), 2.0f64.to_ne_bytes());
    /// ```
    fn get_param(&self) -> Option<&Vec<f64>> {
        self.param.as_ref()
    }

    /// Returns a reference to the current best parameter vector
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    ///
    /// # let mut state: CobylaState = CobylaState::new();

    /// # assert!(state.best_param.is_none());
    /// # state.best_param = Some(vec![1.0, 2.0]);
    /// # assert_eq!(state.best_param.as_ref().unwrap()[0].to_ne_bytes(), 1.0f64.to_ne_bytes());
    /// # assert_eq!(state.best_param.as_ref().unwrap()[1].to_ne_bytes(), 2.0f64.to_ne_bytes());
    /// let best_param = state.get_best_param();  // Option<&P>
    /// # assert_eq!(best_param.as_ref().unwrap()[0].to_ne_bytes(), 1.0f64.to_ne_bytes());
    /// # assert_eq!(best_param.as_ref().unwrap()[1].to_ne_bytes(), 2.0f64.to_ne_bytes());
    /// ```
    fn get_best_param(&self) -> Option<&Vec<f64>> {
        self.best_param.as_ref()
    }

    /// Sets the termination reason
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat, TerminationReason, TerminationStatus};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # assert_eq!(state.termination_status, TerminationStatus::NotTerminated);
    /// let state = state.terminate_with(TerminationReason::SolverConverged);
    /// # assert_eq!(state.termination_status, TerminationStatus::Terminated(TerminationReason::SolverConverged));
    /// ```
    fn terminate_with(mut self, reason: TerminationReason) -> Self {
        self.termination_status = TerminationStatus::Terminated(reason);
        self
    }

    /// Sets the time required so far.
    ///
    /// # Example
    ///
    /// ```
    /// # extern crate instant;
    /// # use instant;
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat, TerminationReason};
    /// # let mut state: CobylaState = CobylaState::new();
    /// let state = state.time(Some(instant::Duration::new(0, 12)));
    /// # assert_eq!(state.time.unwrap(), instant::Duration::new(0, 12));
    /// ```
    fn time(&mut self, time: Option<instant::Duration>) -> &mut Self {
        self.time = time;
        self
    }

    /// Returns current cost function value.
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # state.cost = Some(vec![12.0]);
    /// let cost = state.get_cost();
    /// # assert_eq!(cost.to_ne_bytes(), 12.0f64.to_ne_bytes());
    /// ```
    fn get_cost(&self) -> Self::Float {
        match self.cost.as_ref() {
            Some(c) => *(c.first().unwrap_or(&f64::INFINITY)),
            None => f64::INFINITY,
        }
    }

    /// Returns current best cost function value.
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # state.best_cost = Some(vec![12.0]);
    /// let best_cost = state.get_best_cost();
    /// # assert_eq!(best_cost.to_ne_bytes(), 12.0f64.to_ne_bytes());
    /// ```
    fn get_best_cost(&self) -> Self::Float {
        match self.best_cost.as_ref() {
            Some(c) => *(c.first().unwrap_or(&f64::INFINITY)),
            None => f64::INFINITY,
        }
    }

    /// Returns target cost function value.
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # state.target_cost = 12.0;
    /// let target_cost = state.get_target_cost();
    /// # assert_eq!(target_cost.to_ne_bytes(), 12.0f64.to_ne_bytes());
    /// ```
    fn get_target_cost(&self) -> Self::Float {
        self.target_cost
    }

    /// Returns current number of iterations.
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # state.iter = 12;
    /// let iter = state.get_iter();
    /// # assert_eq!(iter, 12);
    /// ```
    fn get_iter(&self) -> u64 {
        self.iter
    }

    /// Returns iteration number of last best parameter vector.
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # state.last_best_iter = 12;
    /// let last_best_iter = state.get_last_best_iter();
    /// # assert_eq!(last_best_iter, 12);
    /// ```
    fn get_last_best_iter(&self) -> u64 {
        self.last_best_iter
    }

    /// Returns the maximum number of iterations.
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # state.max_iters = 12;
    /// let max_iters = state.get_max_iters();
    /// # assert_eq!(max_iters, 12);
    /// ```
    fn get_max_iters(&self) -> u64 {
        self.max_iters
    }

    /// Returns the termination status.
    ///
    /// # Example
    ///
    /// ```
    /// # use argmin::core::{State, ArgminFloat, TerminationStatus};
    /// # use cobyla::CobylaState;
    /// # let mut state = CobylaState::new();
    /// let termination_status = state.get_termination_status();
    /// # assert_eq!(*termination_status, TerminationStatus::NotTerminated);
    /// ```
    fn get_termination_status(&self) -> &TerminationStatus {
        &self.termination_status
    }

    /// Returns the termination reason if terminated, otherwise None.
    ///
    /// # Example
    ///
    /// ```
    /// # use argmin::core::{State, ArgminFloat, TerminationReason};
    /// # use cobyla::CobylaState;
    /// # let mut state = CobylaState::new();
    /// let termination_reason = state.get_termination_reason();
    /// # assert_eq!(termination_reason, None);
    /// ```
    fn get_termination_reason(&self) -> Option<&TerminationReason> {
        match &self.termination_status {
            TerminationStatus::Terminated(reason) => Some(reason),
            TerminationStatus::NotTerminated => None,
        }
    }

    /// Returns the time elapsed since the start of the optimization.
    ///
    /// # Example
    ///
    /// ```
    /// # extern crate instant;
    /// # use instant;
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// let time = state.get_time();
    /// # assert_eq!(time.unwrap(), instant::Duration::new(0, 0));
    /// ```
    fn get_time(&self) -> Option<instant::Duration> {
        self.time
    }

    /// Increments the number of iterations by one
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # assert_eq!(state.iter, 0);
    /// state.increment_iter();
    /// # assert_eq!(state.iter, 1);
    /// ```
    fn increment_iter(&mut self) {
        self.iter += 1;
    }

    /// Set all function evaluation counts to the evaluation counts of another `Problem`.
    ///
    /// ```
    /// # use std::collections::HashMap;
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{Problem, State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # assert_eq!(state.counts, HashMap::new());
    /// # state.counts.insert("test2".to_string(), 10u64);
    /// #
    /// # #[derive(Eq, PartialEq, Debug)]
    /// # struct UserDefinedProblem {};
    /// #
    /// # let mut problem = Problem::new(UserDefinedProblem {});
    /// # problem.counts.insert("test1", 10u64);
    /// # problem.counts.insert("test2", 2);
    /// state.func_counts(&problem);
    /// # let mut hm = HashMap::new();
    /// # hm.insert("test1".to_string(), 10u64);
    /// # hm.insert("test2".to_string(), 2u64);
    /// # assert_eq!(state.counts, hm);
    /// ```
    fn func_counts<O>(&mut self, problem: &Problem<O>) {
        for (k, &v) in problem.counts.iter() {
            let count = self.counts.entry(k.to_string()).or_insert(0);
            *count = v
        }
    }

    /// Returns function evaluation counts
    ///
    /// # Example
    ///
    /// ```
    /// # use std::collections::HashMap;
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # assert_eq!(state.counts, HashMap::new());
    /// # state.counts.insert("test2".to_string(), 10u64);
    /// let counts = state.get_func_counts();
    /// # let mut hm = HashMap::new();
    /// # hm.insert("test2".to_string(), 10u64);
    /// # assert_eq!(*counts, hm);
    /// ```
    fn get_func_counts(&self) -> &HashMap<String, u64> {
        &self.counts
    }

    /// Returns whether the current parameter vector is also the best parameter vector found so
    /// far.
    ///
    /// # Example
    ///
    /// ```
    /// # use cobyla::CobylaState;
    /// # use argmin::core::{State, ArgminFloat};
    /// # let mut state: CobylaState = CobylaState::new();
    /// # state.last_best_iter = 12;
    /// # state.iter = 12;
    /// let is_best = state.is_best();
    /// # assert!(is_best);
    /// # state.last_best_iter = 12;
    /// # state.iter = 21;
    /// # let is_best = state.is_best();
    /// # assert!(!is_best);
    /// ```
    fn is_best(&self) -> bool {
        self.last_best_iter == self.iter
    }
}
