# cobyla - a pure Rust implementation

[![tests](https://github.com/relf/cobyla/workflows/tests/badge.svg)](https://github.com/relf/cobyla/actions?query=workflow%3Atests)
[![crates.io](https://img.shields.io/crates/v/cobyla)](https://crates.io/crates/cobyla)
[![docs](https://docs.rs/cobyla/badge.svg)](https://docs.rs/cobyla)

COBYLA is an algorithm for minimizing a function of many variables. The method is derivatives-free (only the function values are needed) 
and take into account constraints on the variables. The algorithm is described in:

  > M.J.D. Powell, "A direct search optimization method that models the objective and constraint functions by linear interpolation," in 
  > Advances in Optimization and Numerical Analysis Mathematics and Its Applications, vol. 275 (eds. Susana Gomez and Jean-Pierre Hennart), 
  > Kluwer Academic Publishers, pp. 51-67 (1994).

The algorithm comes into two flavours :
* As an [argmin solver](), the Rust code was generated from the C code from [here](https://github.com/emmt/Algorithms/tree/master/cobyla) 
* As a function `minimize`, the Rust code was generated from the C code of the [NLopt](https://github.com/stevengj/nlopt) project (version 2.7.1)  

In both cases, an initial transpilation was done with [c2rust](https://github.com/immunant/c2rust) then the code was manually edited to make it work. 
The callback mechanismn is inspired from the Rust binding of NLopt, namely [rust-nlopt](https://github.com/adwhit/rust-nlopt)

## Example

```bash
cargo run --example paraboloid
```

## Related projects

* [rust-nlopt](https://github.com/adwhit/rust-nlopt): the Rust binding of the [NLopt project](https://nlopt.readthedocs.io)
* [argmin](https://github.com/argmin-rs/argmin): the pure-Rust optimization framework
* [slsqp](https://github.com/relf/slsqp): a pure Rust implementation of the SLSQP algorithm. 

## License

The project is released under MIT License.
