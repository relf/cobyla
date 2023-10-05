# cobyla

[![tests](https://github.com/relf/cobyla/workflows/tests/badge.svg)](https://github.com/relf/cobyla/actions?query=workflow%3Atests)
[![crates.io](https://img.shields.io/crates/v/cobyla)](https://crates.io/crates/cobyla)
[![docs](https://docs.rs/cobyla/badge.svg)](https://docs.rs/cobyla)

COBYLA an algorithm for minimizing a function of many variables. The method is derivatives free (only the function values are needed) 
and take into account constraints on the variables. The algorithm is described in:

  > M.J.D. Powell, "A direct search optimization method that models the objective and constraint functions by linear interpolation," in 
  > Advances in Optimization and Numerical Analysis Mathematics and Its Applications, vol. 275 (eds. Susana Gomez and Jean-Pierre Hennart), 
  > Kluwer Academic Publishers, pp. 51-67 (1994).

## cobyla 0.4.x

COBYLA implementation coming from the NLopt project is now also available (as `nlopt_cobyla`) allowing to compare 
with the initial implementation (`fmin_cobyla`). This new implementation went through the same process of using c2rust
transpilation from initial C implementation.

```bash
cargo run --example paraboloid
```

## cobyla 0.3.x

COBYLA is now also implemented as an argmin::Solver to benefit from [argmin framework](https://github.com/argmin-rs) tooling. See [example](./examples/paraboloid.rs)

```bash
cargo run --example paraboloid
```

## cobyla 0.2.x

COBYLA C code has been translated to Rust using [c2rust](https://github.com/immunant/c2rust) then manually edited.

## cobyla 0.1.x

Rust wrapper for COBYLA optimizer (COBYLA stands for Constrained Optimization BY Linear Approximations). 
COBYLA C code was copied from [here](https://github.com/emmt/Algorithms/tree/master/cobyla) and wrapped 
using the callback trick implemented in [nlopt-rust](https://github.com/adwhit/rust-nlopt) project.

