# cobyla

[![tests](https://github.com/relf/cobyla/workflows/tests/badge.svg)](https://github.com/relf/cobyla/actions?query=workflow%3Atests)
[![crates.io](https://img.shields.io/crates/v/cobyla)](https://crates.io/crates/cobyla)
[![docs](https://docs.rs/cobyla/badge.svg)](https://docs.rs/cobyla)

COBYLA an algorithm for minimizing a function of many variables. The method is derivatives free (only the function values are needed) and take into account constraints on the variables. The algorithm is described in:

  > M.J.D. Powell, "A direct search optimization method that models the objective and constraint functions by linear interpolation," in 
  > Advances in Optimization and Numerical Analysis Mathematics and Its Applications, vol. 275 (eds. Susana Gomez and Jean-Pierre Hennart), 
  > Kluwer Academic Publishers, pp. 51-67 (1994).

## Example

```bash
cargo run --example paraboloid
```