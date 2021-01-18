use cobyla::{fmin_cobyla, CstrFn};

fn paraboloid(x: &[f64], _data: &mut ()) -> f64 {
    10. * (x[0] + 1.).powf(2.) + x[1].powf(2.)
}

fn main() {
    let mut x = vec![1., 1.];

    #[allow(bare_trait_objects)]
    let mut cons: Vec<&CstrFn> = vec![];
    let cstr1 = |x: &[f64]| x[0];
    cons.push(&cstr1);

    let (status, x_opt) = fmin_cobyla(paraboloid, &mut x, &cons, (), 0.5, 1e-4, 200, 1);
    println!("status = {}", status);
    println!("x = {:?}", x_opt);
}
