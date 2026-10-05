# Airway discretization in lumped simulation models of pulmonary gas exchange: how many compartments are enough?

Code and Supplemental Material for the paper of the same title by Emre Dincel
(Istanbul Technical University). The Supplemental Material S1 to S6 is
`Supplemental_Material.pdf` in this repository.

The repository contains a self-contained implementation of the dead-space chain
analyzed in the paper, together with two drivers that regenerate every table the
isolated chain produces. It depends on nothing but base MATLAB. It reproduces
Tables 1 to 3 of the manuscript, the figures quoted in the text of Sections 4.1
to 4.3, and Supplements S1.1 to S1.5, S2, S4, S4.1 and the residues of S6, and it
computes every quantity plotted in Figures 1 to 3, though it does not draw them.

Tested on MATLAB R2024b. No toolboxes required.

## Quick start

```matlab
reproduce_tables         % Tables 1-2, the Section 4.1-4.3 items, about 70 s
reproduce_tables('all')  % adds Table 3, the cost search; about three minutes more
reproduce_tables(0)      % just the exactness checks, instant
reproduce_tables(1)      % one table; also 'K', 'courant', 'convergence'
supplement_checks        % Supplements S1.4, S1.5, S4 and S4.1, about 70 s
```

## What is here

Seven functions, one per job.

| File | What it is |
|---|---|
| `bohr_closed.m` | Eqs. (8) and (9). The dead space a chain delivers, in closed form. Called with a Courant number it gives the discrete result for the scheme that is actually run; called without one it gives the Co → 0 Erlang limit, which is the same expression in the limit. Its second output is the Courant number a run would actually realize. |
| `bohr_chain_sim.m` | The chain itself. Explicit upwind, periodic steady state, both flow profiles, constant or swinging alveolar boundary. Its second output carries the chain contents at the end of each phase, on which Proposition 3 rests, and the step count of each phase can be set directly. |
| `bohr_parcel_sim.m` | The transport-exact reference. Ordered parcel list, no mixing between neighbours. |
| `bohr_matched_sim.m` | The volume-matched schedule of Remark 1. Each step carries exactly one compartment volume, so the Courant number is one whatever the flow does. Exact above the dead space; below it, the partial step ending each phase leaves the residue Supplement S6 reports. |
| `criterion_nstar.m` | Eq. (15). The compartment count a stated accuracy requires. |
| `reproduce_tables.m` | Driver for the main text. Regenerates the tables and prints them beside their assumptions. |
| `supplement_checks.m` | Driver for the supplement: the effective Courant number of a half-sine (S1.4), Proposition 3 against simulation (S1.5), the criterion read in reverse (S4), and the two `K` tables that show why two is a bound and 2.86 is not (S4.1). |

## What each item needs

| Item | Content | Produced by |
|---|---|---|
| Table 1 | Compartment count `N*` from the criterion | `reproduce_tables(1)` |
| Table 2 | Alveolar PCO₂ with V̇A taken from the chain | `reproduce_tables(2)` |
| Table 3 | Cheapest configuration reaching one percent | `reproduce_tables(3)`, slow |
| Tables 4, 5, 6 | Reference model, open and closed loop | **not in this package** |
| Section 4.1, S1.1–S1.3 | Exactness of Eqs. (8) and (9), requested against realized Courant number | `reproduce_tables(0)` |
| Section 4.1 | Periodic factor `K` under square flow, 161 configurations, and its invariance to rate, dead space, PA and I:E | `reproduce_tables('K')` |
| Section 4.1, S1.5 | The residue identity of Proposition 3, and where its hypothesis fails | `supplement_checks(4)` |
| Section 4.2, Figure 2 | Dead-space error against Courant number | `reproduce_tables('courant')` |
| Section 4.3 | Convergence at r = 0.50, and what to reference it to | `reproduce_tables('convergence')` |
| S1.4 | Effective Courant number of a half-sine | `supplement_checks(1)` |
| S2 | Cost of the parcel scheme | `bohr_parcel_sim(0.150, 0.165, struct('swing', s))` with s = 0, 0.5, 2 and 4; the second output holds the buffer and the updates per breath |
| S4 | The criterion read in reverse, `r_max` for a given `N` | `supplement_checks(2)` |
| S4.1 | `K` at the peak and at half the peak | `supplement_checks(3)` |
| S6 | Residues of the volume-matched schedule | `bohr_matched_sim`, as its help text shows |
| S3, S5 | Reference model | **not in this package** |

Tables 4, 5 and 6 come from a closed-loop Simulink model of respiratory mechanics,
gas exchange, blood chemistry and chemoreflex control. That model is not
distributed here. Everything the paper concludes about the chain itself, and
every theoretical result, is reproducible from this repository alone; the
reference model shows what the error does once the chain sits inside a working
physiological model, and its role is stated as such in the paper.

## Three things worth knowing before you use this

**The step schedule is not a detail, and setting a step is not reporting it.**
Under square flow the flow is constant within each phase, so a fixed step holds
the Courant number exactly at the requested value. Under a half-sine the step is
set from the Courant number at *peak* flow and then held fixed, so the realized
value tracks the instantaneous flow and lies below that peak everywhere but the
peak itself. A half-sine run is therefore *identified* by its peak, and results
are not comparable across profiles without saying which was used. Every table
above prints its own profile.

But when such a run is compared against `bohr_closed` or fed to
`criterion_nstar`, the Courant number to pass is **half the peak**, not the peak:
that is the value the chain effectively runs at, and `supplement_checks(1)`
measures the equivalent at 0.49 to 0.58 of the peak over peaks from 0.10 to 0.99.
Reading the closed form at the peak instead inflates the discrepancy without
bound as `N` grows, which is the whole content of `supplement_checks(3)`. At the
peak the ratio reaches 2.86 and is still climbing at `N = 80`, while at half the
peak it stays below 1.94 and the bound of Proposition 3 holds.

**Requested and realized Courant numbers differ.** A phase must contain a whole
number of steps, so the value a run realizes is not in general the value asked
for. Evaluating `bohr_closed` at the requested value leaves residuals of order
1e-6 against a single-pass simulation, and at the realized value they return to
roundoff. That is why the realized value is the function's second output:

```matlab
[~, Cr] = bohr_closed(5, 0.150, 0.300, 0.70);   % asked 0.70, gets 0.7143
v       = bohr_closed(5, 0.150, 0.300, Cr);     % now exact
```

`reproduce_tables(0)` demonstrates it across ten configurations. The same care is
needed when forming `K`. Against the requested value it reads as high as 2.0143,
above the bound the paper proves, and that excess is the rounding rather than the
chain.

**`min(VD,VT)` is only sometimes the right reference.** It is the
transport-exact answer when the alveolar partial pressure is held constant across
the breath. Give the alveolar compartment an intra-breath swing, as every
closed-loop model has, and the transport-exact answer moves below it by an amount
that owes nothing to discretization. Measuring against `min(VD,VT)` then charges
every scheme for an error none of them can remove, and produces a spurious
optimum near ten compartments. `reproduce_tables('convergence')` shows this directly. Use
`bohr_parcel_sim` as the reference in that case.

## Reuse

The two functions most likely to be useful outside this paper are `bohr_closed`,
which gives the dead space any upwind chain will deliver without running it and
which carries the realized-Courant correction that cost the most time here, and
`criterion_nstar`, which inverts it to a compartment count.
