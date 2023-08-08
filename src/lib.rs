/// One of the most basic models we can have.
pub mod sir;

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
