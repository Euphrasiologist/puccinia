// okay, we start with the simplest model
// hopefully we can make the machinery more generic
// with a nice API at the end.

/// A model trait with methods to be implemented
trait Model {
    /// Check the parameters
    fn check(&self);
    /// The differential equations
    fn diff(&self, pop: Pop, dpop: &mut DPop);
    /// The Runge-Kutta routine
    fn runge_kutta()
}

/// The basic input parameters
#[derive(Clone, Copy)]
struct BasicParameters {
    /// Beta - is the transmission rate and incorporates the
    /// encounter rate between susceptible and infectious
    /// individuals together with the probability of transmission.
    beta: f64,
    /// Gamma - is called the removal or recovery rate, though often
    /// we are more interested in its reciprocal (1/γ) which determines
    /// the average infectious period.
    gamma: f64,
    /// s - is the initial proportion of the population that are susceptible.
    s: f64,
    /// i - is the initial proportion of the population that are infectious.
    i: f64,
}

impl BasicParameters {
    fn check(&self) {
        let b_pos = self.beta.is_sign_positive();
        todo!()
    }
}

/// The maximum time the algorithm should run for
const MAXTIME: f64 = 70.0;

/// A population structure with SIR parameters
#[derive(Copy, Clone, Debug, Default)]
struct Pop {
    s: f64,
    i: f64,
    r: f64,
}

/// A (derived?) population structure?
/// Seems redundant, just use the above...
#[derive(Clone, Copy, Debug, Default)]
struct DPop {
    s: f64,
    i: f64,
    r: f64,
}

// diff is where the differential equations are defined.
fn diff(bp: BasicParameters, pop: Pop, dpop: &mut DPop) {
    let temp_s = pop.s;
    let temp_i = pop.i;
    let _temp_r = pop.r;

    let BasicParameters {
        beta,
        gamma,
        s: _,
        i: _,
    } = bp;

    dpop.s = -beta * temp_s * temp_i;
    dpop.i = beta * temp_s * temp_i - gamma * temp_i;
    dpop.r = gamma * temp_i;
}

fn runge_kutta(
    step: f64,
    bp: BasicParameters,
    dpop: &mut DPop,
    s: &mut f64,
    i: &mut f64,
    r: &mut f64,
) {
    let mut dpop1 = DPop::default();
    let mut dpop2 = DPop::default();
    let mut dpop3 = DPop::default();
    let mut dpop4 = DPop::default();
    let mut tmp_pop = Pop::default();
    let initial_pop = Pop {
        s: *s,
        i: *i,
        r: *r,
    };

    diff(bp, initial_pop, dpop);

    fn rk_inner(dpopn: &mut DPop, dpop: &mut DPop, tmp_pop: &mut Pop, initial_pop: Pop, step: f64) {
        dpopn.s = dpop.s;
        dpopn.i = dpop.i;
        dpopn.r = dpop.r;

        tmp_pop.s = initial_pop.s + step * dpopn.s / 2.0;
        tmp_pop.i = initial_pop.i + step * dpopn.i / 2.0;
        tmp_pop.r = initial_pop.r + step * dpopn.r / 2.0;
    }

    rk_inner(&mut dpop1, dpop, &mut tmp_pop, initial_pop, step);
    diff(bp, tmp_pop, dpop);

    rk_inner(&mut dpop2, dpop, &mut tmp_pop, initial_pop, step);
    diff(bp, tmp_pop, dpop);

    rk_inner(&mut dpop3, dpop, &mut tmp_pop, initial_pop, step);
    diff(bp, tmp_pop, dpop);

    dpop4.s = dpop.s;
    dpop4.i = dpop.i;
    dpop4.r = dpop.r;

    tmp_pop.s =
        initial_pop.s + (dpop1.s / 6.0 + dpop2.s / 3.0 + dpop3.s / 3.0 + dpop4.s / 4.0) * step;
    tmp_pop.i =
        initial_pop.i + (dpop1.i / 6.0 + dpop2.i / 3.0 + dpop3.i / 3.0 + dpop4.i / 4.0) * step;
    tmp_pop.r =
        initial_pop.r + (dpop1.r / 6.0 + dpop2.r / 3.0 + dpop3.r / 3.0 + dpop4.r / 4.0) * step;

    *s = tmp_pop.s;
    *i = tmp_pop.i;
    *r = tmp_pop.r;
}

fn main() {
    // let's do the bulk of the calc here
    let bp = BasicParameters {
        beta: 520.0 / 365.0,
        gamma: 1.0 / 7.0,
        s: 1.0 - 1e-6,
        i: 1e-6,
    };

    let mut s = bp.s;
    let mut i = bp.i;
    let mut r = 1.0 - s - i;

    let step = 0.01 * ((bp.beta + bp.gamma) * s);

    let mut every = (1.0 / ((bp.beta + bp.gamma) * s)).log10().floor().powi(10);

    while MAXTIME / every > 10000.0 {
        every *= 10.0;
    }

    let mut t = 0.0;

    let mut dpop = DPop::default();

    while t < MAXTIME {
        runge_kutta(step, bp, &mut dpop, &mut s, &mut i, &mut r);
        t += step;

        if (t / every).floor() > ((t - step) / every).floor() {
            println!("{t}\t{s}\t{i}\t{r}");
        }
    }
}
