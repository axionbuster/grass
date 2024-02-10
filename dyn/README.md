# dyn: Primitives of Dynamics.

## Overview.

This directory lists mostly independent building blocks used to calculate the gravitational potential between pairs of
circular particles of uniform mass distribution.

The main goals are to:

- Compute the acceleration of a test particle due to the presence of another particle in space while considering the
  case where the two particles are too close, at which moment the inherent singularity in Newtonian gravity becomes
  significant.
- Integrate particles' motion in a numerically stable way while preserving the kinetic and potential energies.

## Acceleration

The Gravity class computes the gravitational interaction between two circular particles.

When the two circles are disjoint, acceleration is computed using the usual Newtonian gravitational formula. When one
circle fully contains the other, neither particle will feel any acceleration due to the other due to Newton's shell
theorem, which says that the gravitational potential is zero inside a radially symmetrical body in the vacuum. Outside
it, the gravitational potential is equal to the case of a point particle.

When two circles touch each other, numerical integration is used to compute the net force on the test particle due to
the source particle. The test particle is the one that feels the net force. So, the test particle is divided into small
regions of identical areas (and thus having identical masses).

Then, the shell theorem is applied: Each region is imagined to be so small that it can be approximated as a single
particle, say p.

If p is outside the source particle's circle, then, by the shell theorem, the case of point particles apply (the usual
gravity formula is used to calculate the force that p itself feels due to the source particle). If p is inside it, then,
again, by the shell theorem, p will feel zero force.

After looping over many such points (p), the small net forces that each "p" feels (dF) are added up to produce the net
force (F) that the test particle itself feels.

**Reminder**: Such "p" points make up the body of the test particle.

## Integration

Yoshida's fourth-order symplectic integrator is provided.

This integrator is "symplectic," meaning it tends to preserve the energy of a system of particles. This
energy-preserving property is in contrast to, say, a high-order Runge-Kutta method, which generally does not preserve
energy and will quickly lead to global catastrophes (for example, all particles abnormally clump together in a single
blob or all particles spread out to infinity).

Yoshida's integrator takes in a function that relates a particle's position to its acceleration. So, a
velocity-dependent function cannot be used.

The integrator will call this function three times during each step of the algorithm and "average" the results together
in a certain way to update the position and velocity values. Time is supposed to "have been frozen" across the three
calls, meaning that the three calls are parallel, contingent possibilities from a single common past; there is no
supposition of any cause-and-effect relationship between the three calls. At least one of the three calls will use a
negative multiple of the step size. Regarding the model time, it steps both "backward" and "forward" in time.

**Reminder**: There are precisely three acceleration calls per step.

## Low-Discrepancy Sequence (Halton)

Suppose a set of points has been uniformly distributed inside a circle. Then, in the limit of an infinite number of
points in the circle, any subset of the circle's interior ("disk") will contain approximately the same ratio of points
to the total number of points as the ratio of its area to the total area. However, this guarantee is only made at the
limit of infinite samples.

Large "nicely-shaped" areas that are either extraordinarily dense or sparse for any finite number of points could
exist (and almost always do). This fact is detrimental to the geometric integration used in the Newtonian gravity class
in the case of two intersecting circular particles, which assumes that, given any "nice" region of a circle, the number
of particles in the region, area spanned by the region, and the mass contained in the region make identical ratios to
the same quantities for the whole circle, even in the case of as few as thirty regions.

Low-discrepancy sequences enforce this assumption, unlike random sequences. Low-discrepancy sequences are designed to
distribute finite samples of points as evenly as possible while maintaining a degree of unpredictability.

This quality of having a uniform distribution of finite points saves the number of samples drastically. Whereas about
500 regions are needed for the usual uniform random distribution, using Halton low-discrepancy sequences cuts the
requirement to 50. In practice, 30 regions are enough for reasonable behavior in the few-second scale. Thirty is a good
approximation since particles rarely stay intersecting for that long.

## Header Files

In newton.h (Gravity class):

- The public field "G" is the universal gravitational constant. Its default value is 1, which may be set to any other
  constant during a Gravity object's lifetime.
- The "field" public member function calculates the gravitational acceleration between a pair of particles. It
  internally uses a pre-generated "random" (see "Halton" above) set of points for integration, but that random sequence
  frequently has to be regenerated to avoid accumulation effects.
- The "refresh_disk" public member function regenerates the internal set of points.

In kahan.h (Kahan class):

- Kahan's compensated summation: Floating-point summation can accumulate rounding errors. A compensated summation keeps
  an extra variable that keeps track of the rounding errors during a summation and continuously compensates for it.
  Kahan's summation is fast but presumably not vectorizable and requires subnormal numbers to exist. Hence, I recommend
  that compensated summation be used immediately outside the "hot inner loops" unless numerical analysis says otherwise.

Other files:

- circle.h (Circle class)
- halton.h (Halton class)
- yoshida.h (Yoshida class)
