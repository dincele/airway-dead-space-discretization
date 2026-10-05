function [VDsp, Co_realised] = bohr_closed(N, VD, VT, Co)
%BOHR_CLOSED  Dead space an N-tank chain delivers on a single pass, in closed form.
% 
%   VDSP = BOHR_CLOSED(N, VD, VT) evaluates Eq. (8), the continuous limit. The
%   transit volume is Erlang with shape N and rate N/VD, and
%
%       VD_sp = VT [1 - F_N(VT)] + VD F_{N+1}(VT)
%
%   VDSP = BOHR_CLOSED(N, VD, VT, CO) evaluates Eq. (9), which is what the scheme
%   actually delivers at Courant number CO. Each step advances a given portion of
%   gas by one compartment with probability CO, so crossing N compartments takes
%   a negative-binomial number of steps and the transit volume X has
%
%       E[X]   = VD                       exactly
%       Var[X] = VD^2 (1-Co) / N          exactly
%
%   and the delivered dead space is the limited expectation E[min(X,VT)]. Eq. (8)
%   is the CO -> 0 limit of this, so passing CO = 0 or omitting it are the same
%   thing.
%
%   [VDSP, CO_REALISED] = BOHR_CLOSED(N, VD, VT, CO) also returns the Courant
%   number a run would actually realise. A phase must contain a whole number of
%   steps, so a run asked for CO delivers VT/(n*Vc) with n = round(VT/(CO*Vc)).
%   Evaluating this function at the requested value leaves residuals of order
%   1e-6 against a single-pass simulation; evaluating it at CO_REALISED returns
%   them to roundoff:
%
%       [~, Cr] = bohr_closed(5, 0.150, 0.300, 0.70);
%       v       = bohr_closed(5, 0.150, 0.300, Cr);     % now exact
%
%   The same care is needed when forming the periodic factor K. Against the
%   requested value it reads as high as 2.0143, above the bound Proposition 3
%   proves; that excess is the rounding, not the chain.
%
%   EXACTNESS. Eq. (9) is exact for the scheme when VT is a whole number of step
%   volumes CO*VD/N, which is what CO_REALISED enforces. It also assumes the
%   single-pass premise, a Courant number constant within the phase, and a
%   constant alveolar partial pressure across the expiration.
%
%   TRUNCATION. The sum over step counts is finite in effect: min(k*d,VT)
%   saturates at VT once k*d >= VT, so KMAX below covers the saturation point
%   plus 60 standard deviations of the step count. Whatever mass is left over is
%   assigned to VT, the value min() takes there anyway, so the truncation cannot
%   bias the result upward.
%
%   N, VD, VT in litres. CO in [0,1].
%
%   Examples reproducing figures quoted in Section 4.1:
%       bohr_closed(100, 0.150, 0.165)          % -> 0.148633 l   (Eq. 8)
%       bohr_closed(5,   0.150, 0.300, 0.10)    % Eq. (9)
%
%   See also BOHR_CHAIN_SIM, BOHR_PARCEL_SIM, CRITERION_NSTAR.

if nargin < 4 || isempty(Co), Co = 0; end

if nargout > 1
    Vc = VD/N;
    n  = max(1, round(VT/(max(Co,eps)*Vc)));
    Co_realised = VT/(n*Vc);
end

% ------------------------------------------------ Eq. (8), the Co -> 0 limit
if Co <= 0
    x    = N .* VT ./ VD;
    VDsp = VT .* (1 - gammainc(x, N)) + VD .* gammainc(x, N+1);
    return
end

% ------------------------------------------------ Eq. (9), the discrete form
d    = Co * VD / N;                       % volume carried by one step
kmax = max(N, ceil(VT/d) + 40 + ceil(60*sqrt(N*(1-Co))/Co));
k    = (N:kmax)';

% log pmf of the negative binomial, formed in logs to stay stable at large N
logp = gammaln(k) - gammaln(N) - gammaln(k-N+1) ...
       + N*log(Co) + (k-N).*log1p(-Co);
p    = exp(logp);

VDsp = sum(p .* min(k*d, VT)) + max(0, 1 - sum(p)) * VT;
end
