function [VDb, info] = bohr_parcel_sim(VD, VT, opt)
%BOHR_PARCEL_SIM  Transport-exact reference: the airway as an ordered parcel list.
%
%   VDB = BOHR_PARCEL_SIM(VD, VT) carries the dead space as an ordered list of
%   parcels, each holding a volume and a partial pressure.  Fresh gas is pushed
%   in at one end and an equal volume displaced from the other, and neighbouring
%   parcels of equal composition are merged.  No parcel is ever mixed with its
%   neighbour, so transport is plug flow to machine precision whatever the step.
%   What the step controls is how finely the alveolar boundary is sampled, not
%   how faithfully gas is carried.
%
%   This is the reference the scheme comparison of Section 4.5 is measured
%   against.  It has to be, because once the alveolar boundary swings within the
%   breath, min(VD,VT) is no longer the exact answer and measuring a tolerance
%   against it charges every scheme for an error none of them can remove.
%
%   OPT accepts the same fields as BOHR_CHAIN_SIM, plus
%     .steps   steps per phase (default 400)
%
%   COST CAVEAT.  Against a CONSTANT boundary every arriving parcel merges on
%   contact and the buffer collapses to about 1.8 parcels, which is an artefact
%   of the test problem and not a property of the scheme.  Any swing at all,
%   even half a millimetre of mercury, raises the buffer to a few hundred.  Cost
%   figures taken against a constant boundary are meaningless.
%
%   See also BOHR_CHAIN_SIM.

if nargin < 3, opt = struct; end
if ~isfield(opt,'swing'),  opt.swing  = 0;   end
if ~isfield(opt,'RR'),     opt.RR     = 15;  end
if ~isfield(opt,'IE'),     opt.IE     = 0.5; end
if ~isfield(opt,'PA'),     opt.PA     = 40;  end
if ~isfield(opt,'PI'),     opt.PI     = 0;   end
if ~isfield(opt,'warmup'), opt.warmup = 200; end
if ~isfield(opt,'steps'),  opt.steps  = 400; end

T  = 60/opt.RR;
Ti = T*opt.IE/(1+opt.IE);
Te = T - Ti;
qI = pi*VT/(2*Ti);  qE = pi*VT/(2*Te);      % half-sine peak flows

dt  = min(Ti,Te)/opt.steps;
nI  = ceil(Ti/dt);  nE = ceil(Te/dt);
dtI = Ti/nI;  dtE = Te/nE;

v = VD;  p = opt.PA;                         % one parcel, wholly alveolar
maxbuf = 1;  totbuf = 0;  nsamp = 0;

for br = 1:(opt.warmup+1)
    measure = (br == opt.warmup+1);
    if measure, num = 0; den = 0; end

    for k = 1:nI
        dV = qI*sin(pi*((k-0.5)*dtI)/Ti)*dtI;
        [v,p] = push_front(v, p, dV, opt.PI);
        [v,p] = drop_back(v, p, dV);
    end

    for k = 1:nE
        dV = qE*sin(pi*((k-0.5)*dtE)/Te)*dtE;
        PA = opt.PA;
        if opt.swing ~= 0
            PA = opt.PA + opt.swing*sin(2*pi*((k-0.5)*dtE)/Te);
        end
        [v,p] = push_back(v, p, dV, PA);
        [v,p,mass,vol] = drop_front(v, p, dV);
        if measure, num = num + mass;  den = den + vol; end
        maxbuf = max(maxbuf, numel(v));  totbuf = totbuf + numel(v);  nsamp = nsamp + 1;
    end
end

VDb = VT*(opt.PA - num/den)/(opt.PA - opt.PI);

if nargout > 1
    info.steps_per_breath = nI + nE;
    info.buffer_mean      = totbuf/max(1,nsamp);
    info.buffer_max       = maxbuf;
    info.updates_per_breath = (nI+nE)*info.buffer_mean;
end
end

% ------------------------------------------------------------------ helpers
function [v,p] = push_front(v, p, dV, val)
if ~isempty(v) && abs(p(1)-val) < 1e-12, v(1) = v(1) + dV;
else, v = [dV v];  p = [val p]; end
end

function [v,p] = push_back(v, p, dV, val)
if ~isempty(v) && abs(p(end)-val) < 1e-12, v(end) = v(end) + dV;
else, v = [v dV];  p = [p val]; end
end

function [v,p,mass,vol] = drop_front(v, p, dV)
mass = 0;  vol = 0;  rem = dV;
while rem > 1e-15
    t = min(rem, v(1));
    mass = mass + t*p(1);  vol = vol + t;
    v(1) = v(1) - t;  rem = rem - t;
    if v(1) <= 1e-15 && numel(v) > 1, v(1) = [];  p(1) = [];
    elseif v(1) <= 1e-15, break; end
end
end

function [v,p] = drop_back(v, p, dV)
rem = dV;
while rem > 1e-15
    t = min(rem, v(end));
    v(end) = v(end) - t;  rem = rem - t;
    if v(end) <= 1e-15 && numel(v) > 1, v(end) = [];  p(end) = [];
    elseif v(end) <= 1e-15, break; end
end
end
