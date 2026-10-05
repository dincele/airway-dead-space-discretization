function [VDb, info] = bohr_chain_sim(N, VD, VT, Co, opt)
%BOHR_CHAIN_SIM  Periodic-steady-state Bohr dead space of an N-compartment chain.
% 
%   VDB = BOHR_CHAIN_SIM(N, VD, VT, CO) integrates Eq. (1) with explicit Euler
%   until the chain is in periodic steady state, then measures the Bohr dead
%   space on the following breath from the flow-weighted mixed expired partial
%   pressure:
%
%       Pbar_E   = int_exp P_1 |Vdot| dt  /  int_exp |Vdot| dt
%       VD_Bohr  = VT (P_A - Pbar_E) / (P_A - P_I)
%
%   [VDB, INFO] = ... also returns the step counts and the realised Courant
%   number, which differs from the requested one whenever the phase does not
%   contain a whole number of steps.
%
%   OPT is an optional struct:
%     .profile  'square' (default) or 'halfsine'
%     .swing    peak-to-peak/2 amplitude of the alveolar boundary within
%               expiration, mmHg; 0 (default) holds it constant
%     .RR       respiratory rate, min^-1                        (default 15)
%     .IE       inspiratory:expiratory ratio as I/E             (default 0.5)
%     .PA       alveolar partial pressure, mmHg                 (default 40)
%     .PI       inspired partial pressure, mmHg                 (default 0)
%     .warmup   breaths before the measured one                 (default 400)
%     .nI, .nE  steps in inspiration and in expiration; override the counts
%               derived from CO, so that the two phases can be given the same
%               Courant sequence or different ones on purpose (Supplement S1.5)
%
%   STEP SCHEDULE.  Under square flow the flow is constant within each phase, so
%   a fixed step holds Co exactly at its requested value and the step count is
%   rounded to the nearest integer.  Under a half-sine the step is set from Co at
%   PEAK flow and then held fixed, so the realised Courant number tracks the
%   instantaneous flow and lies below that peak everywhere but the peak itself.
%   CO is therefore a peak value here, and a half-sine run is identified by it.
%
%   Setting the step and reporting it are separate matters.  When such a run is
%   compared against BOHR_CLOSED or fed to CRITERION_NSTAR, the Courant number to
%   use is HALF the peak, which is the value the chain effectively runs at; the
%   measured equivalent is 0.49 to 0.58 of the peak (Supplement S1.4).  Reading
%   the closed form at the peak instead inflates the discrepancy without bound as
%   N grows.  Section 4.1 of the paper and Supplement S4.1 set this out.
%
%   Volumes in litres, pressures in mmHg.
%
%   INFO.CO_REALISED is the Courant number the run actually realises under square
%   flow, which differs from the requested one because the step count must be a
%   whole number. It is NaN under a half-sine, where no single value describes
%   the phase. See BOHR_CLOSED for why this matters.  With .nI and .nE set it is
%   the inspiratory value.
%
%   INFO.CE and INFO.R are the chain contents at the end of expiration and at the
%   end of inspiration of the measured breath, each as an alveolar-gas volume,
%   sum((P - P_I)/(P_A - P_I)) * VD/N.  Proposition 3 rests on VD_Bohr = Ce - R,
%   which is conservation, and VD - Ce = R, which needs matched Courant
%   sequences; together they give VD - VD_Bohr = 2R.  SUPPLEMENT_CHECKS(4)
%   tests both.
%
%   See also BOHR_CLOSED, BOHR_PARCEL_SIM.

if nargin < 5, opt = struct; end
if ~isfield(opt,'profile'), opt.profile = 'square'; end
if ~isfield(opt,'swing'),   opt.swing   = 0;        end
if ~isfield(opt,'RR'),      opt.RR      = 15;       end
if ~isfield(opt,'IE'),      opt.IE      = 0.5;      end
if ~isfield(opt,'PA'),      opt.PA      = 40;       end
if ~isfield(opt,'PI'),      opt.PI      = 0;        end
if ~isfield(opt,'warmup'),  opt.warmup  = 400;      end
if ~isfield(opt,'nI'),      opt.nI      = [];       end
if ~isfield(opt,'nE'),      opt.nE      = [];       end

T  = 60/opt.RR;
Ti = T*opt.IE/(1+opt.IE);
Te = T - Ti;
Vc = VD/N;
halfsine = strcmpi(opt.profile,'halfsine');

if halfsine
    % peak flows of a half-sine whose integral over the phase is VT
    qI = pi*VT/(2*Ti);  qE = pi*VT/(2*Te);
    dt = Co*Vc/max(qI,qE);              % step from the peak Courant number
    nI = max(2, ceil(Ti/dt));  nE = max(2, ceil(Te/dt));
else
    qI = VT/Ti;  qE = VT/Te;            % constant flow within each phase
    nI = max(1, round(Ti/(Co*Vc/qI)));
    nE = max(1, round(Te/(Co*Vc/qE)));
end
if ~isempty(opt.nI), nI = opt.nI; end
if ~isempty(opt.nE), nE = opt.nE; end
dtI = Ti/nI;  dtE = Te/nE;

P = linspace(opt.PI, opt.PA, N);        % any initial state in [PI,PA] will do
for br = 1:(opt.warmup+1)
    measure = (br == opt.warmup+1);
    if measure, num = 0; den = 0; end

    for k = 1:nI                                    % inspiration, mouth inward
        q = qI;  if halfsine, q = qI*sin(pi*((k-0.5)*dtI)/Ti); end
        P = advect(P, q*dtI/Vc, opt.PI, +1);
    end
    if measure, R = sum((P - opt.PI)/(opt.PA - opt.PI))*Vc; end

    for k = 1:nE                                    % expiration, alveolus inward
        q = qE;  if halfsine, q = qE*sin(pi*((k-0.5)*dtE)/Te); end
        PA = opt.PA;
        if opt.swing ~= 0
            PA = opt.PA + opt.swing*sin(2*pi*((k-0.5)*dtE)/Te);
        end
        if measure, num = num + P(1)*q*dtE;  den = den + q*dtE; end
        P = advect(P, q*dtE/Vc, PA, -1);
    end
    if measure, Ce = sum((P - opt.PI)/(opt.PA - opt.PI))*Vc; end
end

VDb = VT*(opt.PA - num/den)/(opt.PA - opt.PI);

if nargout > 1
    info.nI = nI;  info.nE = nE;
    if halfsine
        % the flow varies within the phase, so no single realised value exists
        info.Co_realised = NaN;
    else
        info.Co_realised = VT/(nI*Vc);
    end
    info.Ce = Ce;  info.R = R;
end
end

% -------------------------------------------------------------------------
function P = advect(P, a, Pboundary, direction)
%ADVECT  One explicit-upwind step at Courant number a, boundary value Pboundary.
%   direction +1 advects from compartment 1 (mouth) inward; -1 from N inward.
if direction < 0, P = fliplr(P); end
n = numel(P);
F = [Pboundary, P];
P = P - a*(F(2:n+1) - F(1:n));
if direction < 0, P = fliplr(P); end
end
