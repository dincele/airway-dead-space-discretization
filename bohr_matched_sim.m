function [VDb, info] = bohr_matched_sim(N, VD, VT, opt)
%BOHR_MATCHED_SIM  The volume-matched step schedule of Remark 1.
%
%   VDB = BOHR_MATCHED_SIM(N, VD, VT) advances the chain in steps that each
%   carry exactly one compartment volume VD/N, so the Courant number is one
%   whatever the flow waveform does and the update is a shift register.  Where
%   the phase volume is not a whole number of compartment volumes, the phase
%   ends with one partial step at the leftover Courant number.
%
%   This is the scheme of Kooi and Huygen, reached there by integration along
%   characteristics.  It is exact above the dead space, where the Bohr
%   measurement, Eq. (6), returns min(VD,VT) to machine precision.  Below the
%   dead space it is exact only when VT is a
%   whole number of compartment volumes; otherwise the partial step ending
%   each phase leaves a residue, which is what this function is used to
%   measure in Supplement S6.
%
%   Because Eq. (2) makes the chain a function of passed volume alone, no flow
%   waveform, respiratory rate or I:E ratio enters.  OPT accepts only
%     .PA       alveolar partial pressure, mmHg                 (default 40)
%     .PI       inspired partial pressure, mmHg                 (default 0)
%     .warmup   breaths before the measured one                 (default 400)
%
%   INFO returns the step counts and the Courant number of the partial step.
%
%   Reproduces the alveolar volumes quoted in Supplement S6:
%       for VT = [0.030 0.045 0.100 0.140]
%           1000*(VT - bohr_matched_sim(5, 0.150, VT))
%       end
%
%   See also BOHR_CHAIN_SIM, BOHR_PARCEL_SIM.

if nargin < 4, opt = struct; end
if ~isfield(opt,'PA'),     opt.PA     = 40;  end
if ~isfield(opt,'PI'),     opt.PI     = 0;   end
if ~isfield(opt,'warmup'), opt.warmup = 400; end

Vc  = VD/N;
nfull = floor(VT/Vc + 1e-12);          % whole compartment volumes in a phase
frac  = VT/Vc - nfull;                 % the partial step that ends it
if frac < 1e-12, frac = 0; end
a = [ones(1,nfull), frac];             % the phase's Courant sequence
a = a(a > 0);

P = linspace(opt.PI, opt.PA, N);
for br = 1:(opt.warmup+1)
    measure = (br == opt.warmup+1);
    if measure, num = 0; den = 0; end

    for k = 1:numel(a)                              % inspiration, mouth inward
        P = advect(P, a(k), opt.PI, +1);
    end

    for k = 1:numel(a)                              % expiration, alveolus inward
        if measure                                  % P_1 before the update
            num = num + P(1)*a(k)*Vc;  den = den + a(k)*Vc;
        end
        P = advect(P, a(k), opt.PA, -1);
    end
end

VDb = VT*(opt.PA - num/den)/(opt.PA - opt.PI);

if nargout > 1
    info.steps_per_phase = numel(a);
    info.full_steps      = nfull;
    info.partial_Co      = frac;
end
end

% -------------------------------------------------------------------------
function P = advect(P, a, Pboundary, direction)
if direction < 0, P = fliplr(P); end
n = numel(P);
F = [Pboundary, P];
P = P - a*(F(2:n+1) - F(1:n));
if direction < 0, P = fliplr(P); end
end
