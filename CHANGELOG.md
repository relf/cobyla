# Changelog

## [Unreleased]

## [0.8.0] - 2025-10-20

* Gate COBYLA argmin COBYLA solver behind `argmin` feature.
  Migration guide: If you use `CobylaSolver` change your Cargo.toml as follows:
  `cobyla = { version = "0.7", features = ["argmin"] }`

## [0.7.0] - 2025-10-04

* Upgrade to argmin 0.11.0
* Bump to edition 2024
* Use `web-time` instead of deprecated `instant`

## [0.6.0] - 2024-04-22

* Upgrade to argmin 0.10.0:
  * new `serde1` feature to enable serialization/deserialization (required for checkpointing)
  * add `argmin-observer-slog` dependency

## [0.5.1] - 2023-10-19

* Fix default constraint tolerance (set to zero, was 2e-4)

## [0.5.0] - 2023-10-13

* Remove `fmin_cobyla` implementation as `nlopt_cobyla` now renamed `minimize` based on NLopt implementation is more powerful
Nevertheless Cobyla argmin solver is still based on initial implementation, not on NLopt one.
* Remove gradient from `ObjFn` trait which is now renamed `Func`
* `minimize` API rework to make it more Rusty.
* Not all options of NLopt version are managed. Currently `cobyla::minimize` handles `xinit`, `bounds`, `max_eval`, `ftol_rel`,
`ftol_abs`, `xtol_rel`, `xtol_abs` and only inequality contraints (like c >= 0). See documentation for further details.

## [0.4.0] - 2023-10-06

COBYLA implementation coming from the NLopt project is now also available (as `nlopt_cobyla`) allowing to compare
with the initial implementation (`fmin_cobyla`). This new implementation went through the same process of using c2rust
transpilation from initial C implementation.

## [0.3.0] - 2023-01-28

COBYLA is now also implemented as an argmin::Solver to benefit from [argmin framework](https://github.com/argmin-rs) tooling. See [example](./examples/paraboloid.rs)

## [0.2.0] - 2023-01-09

COBYLA C code has been translated to Rust using [c2rust](https://github.com/immunant/c2rust) then manually edited.

## [0.1.0] - 2022-05-06

Rust wrapper for COBYLA optimizer (COBYLA stands for Constrained Optimization BY Linear Approximations).
COBYLA C code was copied from [here](https://github.com/emmt/Algorithms/tree/master/cobyla) and wrapped
using the callback trick implemented in [nlopt-rust](https://github.com/adwhit/rust-nlopt) project.
