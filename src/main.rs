// okay, we start with the simplest model
// hopefully we can make the machinery more generic
// with a nice API at the end.

/// A model trait with methods to be implemented
/// P are the population parameters
/// D are the differential population parameters
trait Model<P, D> {
    /// The maximum time the algorithm should run for.
    const MAXTIME: f64;
    /// Check the parameters
    fn check(&self);
    /// The differential equations
    fn diff(&mut self, pop: P, dpop: &mut D);
    /// The Runge-Kutta routine
    fn runge_kutta(&mut self, step: f64, dpop: &mut D);
}

/// The basic input parameters
#[derive(Clone, Copy)]
struct SIR {
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
    /// r - the recovered cohort (1 - s - i)
    r: f64,
}

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

impl Model<Pop, DPop> for SIR {
    const MAXTIME: f64 = 70.0;

    fn check(&self) {
        todo!()
    }

    fn diff(&mut self, pop: Pop, dpop: &mut DPop) {
        let temp_s = pop.s;
        let temp_i = pop.i;
        let _temp_r = pop.r;

        dpop.s = -self.beta * temp_s * temp_i;
        dpop.i = self.beta * temp_s * temp_i - self.gamma * temp_i;
        dpop.r = self.gamma * temp_i;
    }

    fn runge_kutta(&mut self, step: f64, dpop: &mut DPop) {
        let mut dpop1 = DPop::default();
        let mut dpop2 = DPop::default();
        let mut dpop3 = DPop::default();
        let mut dpop4 = DPop::default();
        let mut tmp_pop = Pop::default();
        let initial_pop = Pop {
            s: self.s,
            i: self.i,
            r: self.r,
        };

        self.diff(initial_pop, dpop);

        fn rk_inner(
            dpopn: &mut DPop,
            dpop: &mut DPop,
            tmp_pop: &mut Pop,
            initial_pop: Pop,
            step: f64,
        ) {
            dpopn.s = dpop.s;
            dpopn.i = dpop.i;
            dpopn.r = dpop.r;

            tmp_pop.s = initial_pop.s + step * dpopn.s / 2.0;
            tmp_pop.i = initial_pop.i + step * dpopn.i / 2.0;
            tmp_pop.r = initial_pop.r + step * dpopn.r / 2.0;
        }

        rk_inner(&mut dpop1, dpop, &mut tmp_pop, initial_pop, step);
        self.diff(tmp_pop, dpop);

        rk_inner(&mut dpop2, dpop, &mut tmp_pop, initial_pop, step);
        self.diff(tmp_pop, dpop);

        rk_inner(&mut dpop3, dpop, &mut tmp_pop, initial_pop, step);
        self.diff(tmp_pop, dpop);

        dpop4.s = dpop.s;
        dpop4.i = dpop.i;
        dpop4.r = dpop.r;

        tmp_pop.s =
            initial_pop.s + (dpop1.s / 6.0 + dpop2.s / 3.0 + dpop3.s / 3.0 + dpop4.s / 4.0) * step;
        tmp_pop.i =
            initial_pop.i + (dpop1.i / 6.0 + dpop2.i / 3.0 + dpop3.i / 3.0 + dpop4.i / 4.0) * step;
        tmp_pop.r =
            initial_pop.r + (dpop1.r / 6.0 + dpop2.r / 3.0 + dpop3.r / 3.0 + dpop4.r / 4.0) * step;

        // TODO: does this work?
        self.s = tmp_pop.s;
        self.i = tmp_pop.i;
        self.r = tmp_pop.r;
    }
}

fn main() {
    // let's do the bulk of the calc here
    let mut bp = SIR {
        beta: 520.0 / 365.0,
        gamma: 1.0 / 7.0,
        s: 1.0 - 1e-6,
        i: 1e-6,
        r: 1.0 - (1.0 - 1e-6) - 1e-6,
    };

    let step = 0.01 * ((bp.beta + bp.gamma) * bp.s);

    let mut every = (1.0 / ((bp.beta + bp.gamma) * bp.s))
        .log10()
        .floor()
        .powi(10);

    while SIR::MAXTIME / every > 10000.0 {
        every *= 10.0;
    }

    let mut t = 0.0;

    let mut dpop = DPop::default();

    while t < SIR::MAXTIME {
        bp.runge_kutta(step, &mut dpop);

        t += step;

        if (t / every).floor() > ((t - step) / every).floor() {
            let s = bp.s;
            let i = bp.i;
            let r = bp.r;
            println!("{t}\t{s}\t{i}\t{r}");
        }
    }
}
