function [Nstar, predicted] = criterion_nstar(r, tol, Co, VD, K, Nmax)
%CRITERION_NSTAR  Eq. (15): the compartment count a stated accuracy requires.
%
%   NSTAR = CRITERION_NSTAR(R, TOL, CO) returns the smallest N satisfying
%
%       K * e_disc(N, Co, r) <= tol,        e_disc = [min(VD,VT) - VD_sp]/min(VD,VT)
%
%   with e_disc from the discrete closed form and K = 2 by default.  R is the
%   dead-space fraction VD/VT and TOL is the tolerance on the dead space, not on
%   alveolar ventilation; dividing through by the amplification r/(1-r) converts
%   between the two, and the difference is not small.
%
%   [NSTAR, PREDICTED] = ... also returns K*e_disc at the returned count.
%
%   WHY K = 2.  Under square flow at a constant Courant number, Proposition 3
%   bounds the ratio of the periodic deficit to the single-pass one by 2 and it
%   is approached from below, so taking the limit makes the criterion
%   conservative by construction.  It is conservative in fact as well: at all
%   eighteen points of Table 1 the achieved error is below target, by 3 to 25
%   per cent of it at seventeen of them and by 83 per cent at r = 0.30 under the
%   five per cent tolerance, where the count returned is two.
%
%   CO UNDER A HALF-SINE.  The bound holds at the Courant number the chain
%   effectively runs at, not at the peak its step was built from.  A model
%   driving a half-sine should therefore pass HALF its peak Courant number to
%   this function; entering the peak itself makes the criterion
%   anti-conservative by an amount that grows with N.  Section 4.1 of the paper
%   and Supplement S4.1 set this out.
%
%   TWO CAVEATS carried from Section 2.6.  For a FIXED time step Co is
%   proportional to N, so the two cannot be chosen independently and the search
%   must be run along that line rather than at fixed Courant number.  And e_disc
%   is a deficit relative to min(VD,VT), which is the transport-exact answer only
%   when the alveolar partial pressure is held constant across the breath; see
%   Section 4.3 and REPRODUCE_TABLES('convergence').
%
%   The left-hand side is non-increasing in N (Section 2.6 gives the
%   convex-order argument), so the first admissible count is the smallest.  The
%   search below is a plain upward scan returning that count.
%
%   Defaults: CO = 0.10, VD = 0.150 l, K = 2, NMAX = 4000.
%
%   See also BOHR_CLOSED, BOHR_CHAIN_SIM.

if nargin < 3 || isempty(Co),   Co   = 0.10;  end
if nargin < 4 || isempty(VD),   VD   = 0.150; end
if nargin < 5 || isempty(K),    K    = 2.0;   end
if nargin < 6 || isempty(Nmax), Nmax = 4000;  end

VT  = VD/r;
lim = min(VD, VT);

Nstar = NaN;  predicted = NaN;
for N = 1:Nmax
    e = (lim - bohr_closed(N, VD, VT, Co))/lim;
    if K*e <= tol
        Nstar = N;  predicted = 100*K*e;  return
    end
end
warning('criterion_nstar:notFound', ...
        'No N <= %d meets tol = %g at r = %g, Co = %g.', Nmax, tol, r, Co);
end
