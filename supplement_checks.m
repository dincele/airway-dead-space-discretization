function supplement_checks(which)
%SUPPLEMENT_CHECKS  Regenerate the supplement tables the isolated chain produces.
%
%   SUPPLEMENT_CHECKS        runs everything (about seventy seconds)
%   SUPPLEMENT_CHECKS(n)     runs one item:
%       1  S1.4  the effective Courant number of a half-sine
%       2  S4    the criterion read in reverse, r_max for a given N
%       3  S4.1  why two is a bound and 2.86 is not
%       4  S1.5  Proposition 3 against simulation
%
%   Companion to REPRODUCE_TABLES, which covers the tables in the main text
%   and the checks of Supplement S1.1 to S1.3.  Supplements S3 and S5 concern
%   the closed-loop reference model, which is not part of this package; see
%   README.md.
%
%   Baseline throughout: VD = 150 ml, PA = 40 mmHg, PI = 0, rate 15 min^-1,
%   I:E 1:2, constant alveolar boundary.

if nargin < 1, which = 0; end
sel = @(n) which == 0 || which == n;

VD = 0.150;  t0 = tic;

% ===================================================================== S1.4
if sel(1)
fprintf('\n===== S1.4: the effective Courant number of a half-sine =====\n');
fprintf('       N = 5, r = 0.50.  Which CONSTANT Courant number reproduces the\n');
fprintf('       error a half-sine delivers at a given peak?\n\n');
VT = 0.300;  m = min(VD,VT);
fprintf('    %14s %22s %14s\n','peak Co','equivalent constant Co','ratio to peak');
for p = [0.10 0.20 0.30 0.50 0.75 0.90 0.99]
    es = (m - bohr_chain_sim(5,VD,VT,p,struct('profile','halfsine')))/m;
    c  = match_square(5, VD, VT, es);
    fprintf('    %14.2f %22.3f %14.2f\n', p, c, c/p);
end
fprintf('\n    The time average of a half-sine is 2/pi = 0.637 of the peak and the\n');
fprintf('    flow-weighted average pi/4 = 0.785; neither matches.  The first two\n');
fprintf('    rows cover the operating point of the criterion, which is entered at\n');
fprintf('    Co = 0.10 and so corresponds to a half-sine peak of 0.20.  The ratio\n');
fprintf('    is flat at 0.58 there, which is why halving the peak is conservative.\n');
end

% ======================================================================= S4
if sel(2)
fprintf('\n===== S4: the criterion read in reverse =====\n');
fprintf('       Largest dead-space fraction a given count still meets, at\n');
fprintf('       Co = 0.10 under square flow against a constant boundary.\n\n');
fprintf('    %5s %22s %22s\n','N','r_max at one percent','r_max at five percent');
for N = [5 10 20 40]
    fprintf('    %5d', N);
    for tol = [0.01 0.05]
        fprintf('%22.3f', rmax(N, VD, tol));
    end
    fprintf('\n');
end
fprintf('\n    Each entry solves 2*e_disc(N,Co,r) = tol for r by bisection.  The\n');
fprintf('    factor two is the bound of Proposition 3, not a measured value, so\n');
fprintf('    the entries are conservative wherever the bound holds.\n');
end

% ===================================================================== S4.1
if sel(3)
fprintf('\n===== S4.1: why two is a bound and 2.86 is not =====\n');
fprintf('       K under a half-sine at a fixed step, peak Co = 0.10.\n');
fprintf('       Dashes are points whose single-pass deficit is below 1e-6 of the\n');
fprintf('       plug-flow value, where the ratio is formed from two negligibles.\n');
Copk = 0.10;
Ns   = [5 10 20 40 80];
rs   = [0.50 0.68 0.75 0.83 0.86 0.91 0.95];

es = nan(numel(Ns),numel(rs));
for i = 1:numel(Ns)
    for j = 1:numel(rs)
        VT = VD/rs(j);  m = min(VD,VT);
        es(i,j) = (m - bohr_chain_sim(Ns(i),VD,VT,Copk, ...
                                      struct('profile','halfsine')))/m;
    end
end

for frac = [1.00 0.50]
    if frac == 1
        fprintf('\n    (a) e_disc read at the PEAK.  This is the mismatched convention.\n\n');
    else
        fprintf('\n    (b) e_disc read at HALF the peak.  This is the rule the paper gives.\n\n');
    end
    fprintf('    %5s','N');  fprintf('%9.2f', rs);  fprintf('\n');
    best = 0;  where = '';
    for i = 1:numel(Ns)
        fprintf('    %5d', Ns(i));
        for j = 1:numel(rs)
            VT = VD/rs(j);  m = min(VD,VT);
            ed = (m - bohr_closed(Ns(i),VD,VT,Copk*frac))/m;
            if (m - bohr_closed(Ns(i),VD,VT,Copk))/m < 1e-6
                fprintf('%9s','--');  continue
            end
            K = es(i,j)/ed;
            fprintf('%9.3f', K);
            if K > best, best = K; where = sprintf('N=%d r=%.2f',Ns(i),rs(j)); end
        end
        fprintf('\n');
    end
    fprintf('      max = %.3f at %s\n', best, where);
end

fprintf('\n    (a) rises monotonically with N in every column with no sign of a\n');
fprintf('    limit, so no maximum measured on a finite grid can be quoted as a\n');
fprintf('    bound: truncated at N = 40 the same grid gives 2.71, carried to\n');
fprintf('    N = 80 it gives 2.856.  (b) tracks the square-flow values, which\n');
fprintf('    reach 2.000 and are the ones Proposition 3 bounds.  The remedy for a\n');
fprintf('    half-sine is to read the Courant number in the right place, not to\n');
fprintf('    inflate K.\n');
end

% ===================================================================== S1.5
if sel(4)
fprintf('\n===== S1.5: Proposition 3 against simulation =====\n');
fprintf('       Square flow, Co = 0.10.  Ce and R are the chain contents at the\n');
fprintf('       end of expiration and of inspiration, as alveolar-gas volumes.\n\n');
Ns = [5 10 20 40 80];  rs = [0.50 0.68 0.83 0.91 0.95];
d1 = 0;  d2 = 0;  dK = 0;  skipped = {};
for N = Ns
    for r = rs
        VT = VD/r;
        [VDb, in] = bohr_chain_sim(N, VD, VT, 0.10);
        [~, Cr] = bohr_closed(N, VD, VT, 0.10);
        tail = VD - bohr_closed(N, VD, VT, Cr);       % E[(X-VT)^+], single pass
        d1 = max(d1, abs((in.Ce - in.R) - VDb)/VDb);
        d2 = max(d2, abs((VD - in.Ce) - in.R)/VD);
        if tail < 1e-6*VD
            skipped{end+1} = sprintf('N=%d r=%.2f', N, r); %#ok<AGROW>
        else
            dK = max(dK, abs((VD - VDb)/tail - 2*in.R/tail));
        end
    end
end
fprintf('    step 1, VD_Bohr = Ce - R  : largest deviation %.1e of VD_Bohr\n', d1);
fprintf('    step 2, VD - Ce = R       : largest deviation %.1e of VD\n', d2);
fprintf('    K against 2R/E[(X-VT)^+]  : largest difference %.1e\n', dK);
fprintf('    K is not formed where the single-pass deficit is below 1e-6 of VD:\n');
fprintf('    %s\n', strjoin(skipped, ', '));

fprintf('\n    Departure |(VD - VD_Bohr) - 2R| / (VD - VD_Bohr), N = 5, I:E 1:2,\n');
fprintf('    Co = 0.10 at the inspiratory peak.\n\n');
Vc = VD/5;
fprintf('    %-42s', 'r');  fprintf('%10.2f', rs);  fprintf('\n');
names = {'square flow, constant Co, matched phases', ...
         'half-sine, fixed step, nE = 2 nI', ...
         'half-sine, phases matched at nI = nE', ...
         'square flow, phases at 0.05 and 0.20'};
for c = 1:4
    fprintf('    %-42s', names{c});
    for r = rs
        VT = VD/r;
        n  = ceil(pi*VT/(2*0.10*Vc));               % inspiratory peak at 0.10
        switch c
            case 1, o = struct();
            case 2, o = struct('profile','halfsine');
            case 3, o = struct('profile','halfsine','nI',n,'nE',n);
            case 4, o = struct('nI',max(1,round(VT/(0.05*Vc))), ...
                               'nE',max(1,round(VT/(0.20*Vc))));
        end
        [VDb, in] = bohr_chain_sim(5, VD, VT, 0.10, o);
        fprintf('%9.3g%%', 100*abs((VD - VDb) - 2*in.R)/(VD - VDb));
    end
    fprintf('\n');
end
[VDb, in] = bohr_chain_sim(40, VD, VD/0.68, 0.10, struct('profile','halfsine'));
fprintf('\n    Half-sine, fixed step, N = 40, r = 0.68: %.1f%%, so the second row\n', ...
        100*abs((VD - VDb) - 2*in.R)/(VD - VDb));
fprintf('    is a floor and not a ceiling.  The third row restores the identity\n');
fprintf('    on the same waveform that breaks it in the second.\n');
end

fprintf('\nTotal %.0f s\n', toc(t0));
end

% -------------------------------------------------------------------------
function c = match_square(N, VD, VT, target)
%MATCH_SQUARE  Constant Courant number whose square-flow error equals TARGET.
m = min(VD,VT);  lo = 1e-4;  hi = 0.9999;
for k = 1:60
    mid = (lo+hi)/2;
    e = (m - bohr_chain_sim(N,VD,VT,mid,struct('profile','square')))/m;
    if e > target, lo = mid; else, hi = mid; end
end
c = (lo+hi)/2;
end

% -------------------------------------------------------------------------
function r = rmax(N, VD, tol)
%RMAX  Largest dead-space fraction meeting TOL at N compartments, Co = 0.10.
lo = 0.05;  hi = 0.999;                  % error increases with r
for k = 1:60
    mid = (lo+hi)/2;  VT = VD/mid;  m = min(VD,VT);
    e = 2*(m - bohr_closed(N,VD,VT,0.10))/m;
    if e <= tol, lo = mid; else, hi = mid; end
end
r = (lo+hi)/2;
end
