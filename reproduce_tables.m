function reproduce_tables(which)
%REPRODUCE_TABLES  Regenerate every table in the paper that the isolated chain produces.
%
%   REPRODUCE_TABLES         runs everything except the scheme-cost search
%   REPRODUCE_TABLES('all')  runs everything, including the search (slow)
%   REPRODUCE_TABLES(n)      runs one table, n = 1..3, or 0 for the checks
%   REPRODUCE_TABLES(key)    runs one item the paper reports in its text:
%       'K'            Section 4.1, the periodic factor under square flow
%       'courant'      Section 4.2 and Figure 2, error against Courant number
%       'convergence'  Section 4.3, what a convergence study is referenced to
%
%   Covers Tables 1 to 3, the verification checks of Section 4.1 (Supplement
%   S1.1 to S1.3) and the three text items above.  Tables 4, 5 and 6 come from
%   the closed-loop reference model, which is not part of this package; see
%   README.md.
%
%   Baseline throughout: VD = 150 ml, PA = 40 mmHg, PI = 0, rate 15 min^-1,
%   I:E 1:2.  Square flow at a constant Courant number unless stated.

if nargin < 1, which = 'default'; end
run_all = ischar(which) && strcmpi(which,'all');
sel = @(k) strcmpi(which,'default') || run_all || isequal(which,k) || ...
           (ischar(which) && ischar(k) && strcmpi(which,k));

VD = 0.150;  t0 = tic;

% ==================================================================== checks
if sel(0)
fprintf('\n===== Section 4.1 checks (Supplement S1.1 to S1.3) =====\n\n');

fprintf('(a) Eq. (8) reaches the plug-flow limit.  N = 100:\n');
fprintf('    %8s %12s %12s %10s\n','VT (ml)','min(VD,VT)','Eq.(8) ml','error (%)');
for VT = [0.500 0.200 0.165 0.120]
    v = bohr_closed(100, VD, VT);  m = min(VD,VT);
    fprintf('    %8.0f %12.2f %12.4f %9.3f\n', VT*1e3, m*1e3, v*1e3, 100*(v-m)/m);
end
fprintf('    The worst case is the tidal volume closest to VD, not the smallest.\n\n');

fprintf('(b) Eq. (9) is exact for the scheme, at the REALISED Courant number.\n');
fprintf('    %4s %8s %10s %12s | %13s %13s\n', ...
        'N','VT(ml)','requested','realised','resid @ req','resid @ real');
for N = [5 40]
    if N==5, VT = 0.300; else, VT = 0.165; end
    for Co = [0.10 0.30 0.50 0.70 0.90]
        [~, Cr] = bohr_closed(N, VD, VT, Co);
        sp = single_pass(N, VD, VT, Cr);
        r1 = abs(bohr_closed(N,VD,VT,Co) - sp)/sp;
        r2 = abs(bohr_closed(N,VD,VT,Cr) - sp)/sp;
        fprintf('    %4d %8.0f %10.2f %12.4f | %13.2e %13.2e\n', ...
                N, VT*1e3, Co, Cr, r1, r2);
    end
end
fprintf('    Last column at roundoff => the residuals belong to the run,\n');
fprintf('    that is to the integer step count, and not to the formula.\n');
end

% ======================================================== Section 4.1, K
if sel('K')
fprintf('\n===== Section 4.1: the periodic factor K, square flow =====\n\n');
Cos = [0.05 0.10 0.20 0.30 0.50];
fprintf('    %5s %9s %9s %9s %16s\n','N','r=0.50','r=0.75','r=0.86','spread over Co');
for N = [5 10 20]
    fprintf('    %5d', N);  sp = [];
    for r = [0.50 0.75 0.86]
        VT = VD/r;  kk = [];
        for Co = Cos
            [~, Cr] = bohr_closed(N, VD, VT, Co);
            ed = (VD - bohr_closed(N,VD,VT,Cr))/VD;
            if ed < 1e-6, continue; end
            es = (VD - bohr_chain_sim(N,VD,VT,Co))/VD;
            kk(end+1) = es/ed; %#ok<AGROW>
        end
        fprintf('%9.2f', mean(kk));  sp(end+1) = max(kk)-min(kk); %#ok<AGROW>
    end
    fprintf('%12.2f--%.2f\n', min(sp), max(sp));
end

fprintf('\n    Range over all 161 configurations examined\n');
fprintf('    (N = 5..80, r = 0.50..0.95, Co = 0.05..0.50):\n');
K = [];  tag = {};
for N = [5 10 20 40 80]
    for r = [0.50 0.68 0.75 0.83 0.86 0.91 0.95]
        VT = VD/r;
        for Co = Cos
            [~, Cr] = bohr_closed(N, VD, VT, Co);
            ed = (VD - bohr_closed(N,VD,VT,Cr))/VD;
            if ed < 1e-6, continue; end
            K(end+1) = ((VD - bohr_chain_sim(N,VD,VT,Co))/VD)/ed; %#ok<AGROW>
            tag{end+1} = sprintf('N=%d r=%.2f Co=%.2f',N,r,Co);    %#ok<AGROW>
        end
    end
end
[mn,i1] = min(K);  [mx,i2] = max(K);
fprintf('      points   : %d\n', numel(K));
fprintf('      smallest : %.4f  (%s)\n', mn, tag{i1});
fprintf('      largest  : %.4f  (%s)\n', mx, tag{i2});
fprintf('      mean     : %.4f\n', mean(K));
fprintf('      above 2  : %d   <- Proposition 3 bounds this at zero\n', sum(K > 2+1e-9));

fprintf('\n    Invariance of K at N = 5, r = 0.75, Co = 0.10, one setting changed\n');
fprintf('    at a time from the baseline:\n');
vary = {'baseline',          VD,    struct(); ...
        'rate 30 min^-1',    VD,    struct('RR',30); ...
        'dead space 100 ml', 0.100, struct(); ...
        'dead space 220 ml', 0.220, struct(); ...
        'PA 60 mmHg',        VD,    struct('PA',60); ...
        'I:E 1:1',           VD,    struct('IE',1); ...
        'I:E 1:4',           VD,    struct('IE',0.25)};
for i = 1:size(vary,1)
    VDx = vary{i,2};  VTx = VDx/0.75;
    [~, Cr] = bohr_closed(5, VDx, VTx, 0.10);
    ed = (VDx - bohr_closed(5, VDx, VTx, Cr))/VDx;
    es = (VDx - bohr_chain_sim(5, VDx, VTx, 0.10, vary{i,3}))/VDx;
    fprintf('      %-18s K = %.6f\n', vary{i,1}, es/ed);
end
end

% ========================================== Section 4.2 and Figure 2, Co
if sel('courant')
fprintf('\n===== Section 4.2, Figure 2: dead-space error (%%) against Courant number, r = 0.50 =====\n\n');
VT = 0.300;
fprintf('    %6s | %10s %10s %10s | %10s %10s %10s\n', ...
        'Co','sq N=5','N=10','N=40','sine N=5','N=10','N=40');
for Co = [0.10 0.30 0.50 0.75 0.99]
    fprintf('    %6.2f |', Co);
    for prof = {'square','halfsine'}
        for N = [5 10 40]
            o.profile = prof{1};
            b = bohr_chain_sim(N, VD, VT, Co, o);
            fprintf(' %9.3f', 100*(min(VD,VT)-b)/min(VD,VT));
        end
        if strcmp(prof{1},'square'), fprintf(' |'); end
    end
    fprintf('\n');
end
fprintf('\n    Square flow reaches Co = 1 exactly and the error vanishes.  A half-sine\n');
fprintf('    under a fixed step cannot.  The constant Courant number reproducing\n');
fprintf('    the same error is 0.49 to 0.58 of the peak over peaks from 0.10 to\n');
fprintf('    0.99, so a peak of 0.99 behaves like a square wave at about 0.49.\n');
fprintf('    Neither the time average (2/pi) nor the flow-weighted average (pi/4)\n');
fprintf('    predicts this.  A half-sine run is identified here by its PEAK, but\n');
fprintf('    compared against the closed form, or entered into the criterion, at\n');
fprintf('    HALF that peak.  See Section 4.1 and Supplement S4.1.\n');
end

% =================================================================== Table 1
if sel(1)
fprintf('\n===== Table 1: compartment count from Eq. (15), square flow, Co = 0.10 =====\n\n');
rs = [0.30 0.50 0.60 0.68 0.75 0.81 0.86 0.91 0.95];
for tol = [0.05 0.01]
    fprintf('    %.0f%% tolerance on dead space\n', 100*tol);
    fprintf('      %-22s', 'r');       fprintf('%8.2f', rs);      fprintf('\n');
    Ns = zeros(size(rs));  em = Ns;
    for i = 1:numel(rs)
        Ns(i) = criterion_nstar(rs(i), tol, 0.10, VD);
        VT = VD/rs(i);  m = min(VD,VT);
        em(i) = 100*(m - bohr_chain_sim(Ns(i), VD, VT, 0.10))/m;
    end
    fprintf('      %-22s', 'N*');              fprintf('%8d', Ns);  fprintf('\n');
    fprintf('      %-22s', 'error obtained (%)'); fprintf('%8.2f', em); fprintf('\n');
    fprintf('      %-22s', 'error in VA (%)');
    fprintf('%8.1f', em.*rs./(1-rs));  fprintf('\n\n');
end
end

% ======================================================== Section 4.3
if sel('convergence')
fprintf('\n===== Section 4.3: what a convergence study should be referenced to =====\n');
fprintf('       r = 0.50, half-sine, 2 mmHg alveolar swing, peak Co held at 0.10\n\n');
r = 0.50;  VT = VD/r;
o = struct('profile','halfsine','swing',2,'warmup',200);
vp = bohr_parcel_sim(VD, VT, o);
fprintf('    transport-exact reference (parcel run) : %.2f ml\n', vp*1e3);
fprintf('    min(VD,VT)                             : %.2f ml\n\n', min(VD,VT)*1e3);
fprintf('    %5s %14s %16s %18s\n','N','chain (ml)','vs min(VD,VT) %','vs exact transport %');
for N = [5 8 10 20 40 160]
    vc = bohr_chain_sim(N, VD, VT, 0.10, o);
    fprintf('    %5d %14.2f %16.2f %18.2f\n', ...
            N, vc*1e3, 100*(min(VD,VT)-vc)/min(VD,VT), 100*(vp-vc)/vp);
end
fprintf('\n    Against min(VD,VT) the error is NOT monotone in N: it falls to a\n');
fprintf('    minimum near ten compartments and rises again.  That apparent\n');
fprintf('    optimum is a cancellation between the discretisation error and the\n');
fprintf('    boundary term, which have opposite signs there.  It is not accuracy.\n');
end

% =================================================================== Table 2
if sel(2)
fprintf('\n===== Table 2: alveolar PCO2 with VA taken from the chain =====\n');
fprintf('       minute ventilation 7.5 l/min, VCO2 195 ml/min, half-sine flow\n\n');
VE = 7.5;  VCO2 = 0.195;
fs = [15 25 35 45 50];          % stopped at the asymptote, as the paper is
o6 = struct('profile','halfsine','warmup',80);
fprintf('    %-12s','f (min^-1)');  fprintf('%9d', fs);  fprintf('\n');
fprintf('    %-12s','VT (ml)');     fprintf('%9.1f', 1e3*VE./fs);  fprintf('\n');
for N = [5 10 20 40]
    fprintf('    %-12s', sprintf('N = %d', N));
    for f = fs
        VT = VE/f;
        b  = bohr_chain_sim(N, VD, VT, 0.10, o6);
        VA = f*(VT - b);
        fprintf('%9.1f', 863*VCO2/VA);
    end
    fprintf('\n');
end
fprintf('    %-12s','plug flow');
for f = fs
    VT = VE/f;  VA = f*(VT - min(VD,VT));
    if VA <= 1e-9, fprintf('%9s','---'); else, fprintf('%9.1f', 863*VCO2/VA); end
end
fprintf('\n\n    The correct behaviour has a vertical asymptote at 50 min^-1,\n');
fprintf('    where tidal volume equals dead space.  No finite chain has one.\n');
fprintf('    The sweep stops there: beyond it the relation keeps rising in every\n');
fprintf('    row, but towards an alveolar pressure that has no upper bound where\n');
fprintf('    the gas does, 760 - 47 = 713 mmHg once saturated with water vapour.\n');
fprintf('    By the alveolar gas equation an alveolar PCO2 of 125 mmHg exhausts\n');
fprintf('    the alveolar oxygen on room air, so entries above that describe a gas\n');
fprintf('    a patient cannot hold; the paper italicises them for that reason.\n');
fprintf('    They show the ordering, not attainable states.\n');
end

% =================================================================== Table 3
if run_all || isequal(which,3)      % never in the default run: it is slow
fprintf('\n===== Table 3: cheapest configuration reaching one per cent =====\n');
fprintf('       (joint search over N and Co; this takes several minutes)\n\n');
o = struct('profile','halfsine','swing',2,'warmup',200);
fprintf('    %6s %6s %8s %12s %22s\n', ...
        'r','N','Co','upwind cost','parcel cost (see below)');
for r = [0.50 0.83 0.91 0.95]
    VT = VD/r;
    [vp, ip] = bohr_parcel_sim(VD, VT, o);
    best = Inf;  bN = NaN;  bCo = NaN;
    for N = [5 10 20 40 80 160 320]
        for Co = [0.10 0.30 0.50 0.75 0.90 0.98]
            vc = bohr_chain_sim(N, VD, VT, Co, o);
            if abs(vp-vc)/vp <= 0.01
                [~, ic] = bohr_chain_sim(N, VD, VT, Co, o);
                cost = N*(ic.nI + ic.nE);
                if cost < best, best = cost; bN = N; bCo = Co; end
            end
        end
    end
    fprintf('    %6.2f %6d %8.2f %12.0f %22.0f\n', ...
            r, bN, bCo, best, ip.updates_per_breath);
end
fprintf('\n    The first four columns are Table 3.  The selected Courant number is\n');
fprintf('    the finding, not the cost ratio.  The parcel column is NOT in the\n');
fprintf('    table: it is that scheme at a fixed 400 steps per phase rather than\n');
fprintf('    at its own cheapest configuration, so the two are not a like-for-like\n');
fprintf('    ranking, which is why the paper does not print it.  See Supplement S2.\n');
end

fprintf('\nTotal %.0f s\n', toc(t0));
end

% -------------------------------------------------------------------------
function b = single_pass(N, VD, VT, Co)
%SINGLE_PASS  One expiration from a chain holding only inspirate.
Vc = VD/N;  PA = 40;  PI = 0;
n  = max(1, round(VT/(Co*Vc)));
dV = VT/n;
P  = zeros(1,N);  num = 0;  den = 0;
for k = 1:n
    num = num + P(1)*dV;  den = den + dV;
    P = fliplr(P);  F = [PA, P];  P = P - Co*(F(2:N+1)-F(1:N));  P = fliplr(P);
end
b = VT*(PA - num/den)/(PA - PI);
end
