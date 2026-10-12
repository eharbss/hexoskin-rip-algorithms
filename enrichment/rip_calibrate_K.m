function cal = rip_calibrate_K(method, varargin)
%RIP_CALIBRATE_K  Derive the RIP calibration constant K with a stated method.
%
% created 21-09-2026 by Eric Harbour (PLUS)
%
% MODEL
%   VT  is proportional to  dRC + K*dAB          with  K = gain_AB / gain_RC
%   %RC = dRC / (dRC + K*dAB) * 100
%
% WHY THIS FUNCTION EXISTS
%   %RC is steeply sensitive to K: about -7 percentage points per +1.0 in K
%   around K = 3, so K must be known to roughly +/-18% to keep %RC within
%   +/-5 points. On a cycling pilot two FVC-derived baselines from the SAME
%   recording gave K = 1.40 and K = 3.16 -- a 2.3x spread and a 14-point
%   swing in mean %RC. K therefore needs a method with a stated provenance,
%   not whatever window happened to be selected.
%
% METHODS (in descending order of trustworthiness)
%   'isovolume' - ISOCAL. Airway occluded (or glottis held open at a fixed
%                 volume) while the subject shifts volume between rib cage
%                 and abdomen. Lung volume is constant, so
%                     0 = dRC + K*dAB   ->   K = -dRC/dAB
%                 K is the negative reciprocal of the total-least-squares
%                 slope of AB on RC. Zimmerman et al. (1988) reported 6% VT
%                 error for ISOCAL vs 9-23% for posture-change calibration;
%                 Khemani et al. (2014) found ISOCAL clearly superior to QDC,
%                 and increasingly so as breathing became abnormal.
%   'lsq'       - Least-squares / multiple linear regression against a
%                 simultaneous volume reference (mobile spirometer):
%                     VT = a*dRC + b*dAB,  K = b/a
%                 Chadha et al. (1982): 93% of values within +/-10% of
%                 spirometry, the most accurate method tested. Needs breaths
%                 that SPAN a range of RC/AB partitioning -- instruct
%                 alternating chest-dominant, abdomen-dominant and natural
%                 breathing at varied depths, or the regression is
%                 ill-conditioned.
%   'qdc'       - Qualitative diagnostic calibration over ~5 min of quiet
%                 natural breathing (Sackner et al. 1989):
%                     K = SD(dRC) / SD(dAB)   over per-breath amplitudes
%                 NOTE this is the SD of BREATH AMPLITUDES. Using the SD of
%                 the raw samples instead makes K identically the amplitude
%                 ratio and forces mean %RC to 50%.
%   'fvc'       - Ratio of excursions during a forced vital capacity
%                 manoeuvre. Retained for backward compatibility only;
%                 RC and AB do not displace equal volumes during an FVC, so
%                 this is not a gain ratio. A warning is issued.
%
% CAVEATS THAT APPLY TO ALL METHODS
%   * K is posture specific. Barbosa et al. (2012) measured ISOCAL K of
%     1.6 supine, 2.5 sitting, 2.0 standing in the same subjects. Calibrate
%     in the posture you will measure in.
%   * K is not guaranteed to transfer from rest to exercise. Heyde et al.
%     (2014) found gains computed on different data subsets were not
%     equivalent. Repeat the manoeuvre AFTER the bout and compare; a change
%     over ~20% means belt slippage or a genuine gain change.
%   * Within-subject CHANGE in %RC is far more robust to K error than the
%     absolute level. Prefer reporting delta %RC against a reference condition.
%
% SYNTAX
%   cal = rip_calibrate_K('isovolume', tt)                 % tt = [thc abd] timetable
%   cal = rip_calibrate_K('lsq', dRC, dAB, VT)
%   cal = rip_calibrate_K('qdc', dRC, dAB)
%   cal = rip_calibrate_K('fvc', tt)
%   ... , 'fs', 128, 'plots', true
%
% OUTPUT cal
%   .K         calibration constant
%   .method    method used
%   .ci        95% bootstrap confidence interval for K
%   .quality   0-1 heuristic confidence in the estimate
%   .note      plain-language verdict
%   plus method-specific diagnostics (.r, .var_explained, .R2, .a, .b)

% --- argument parsing -------------------------------------------------
% Positional data arguments come first; name-value options follow.
optNames = {'fs','plots'};
kopt = numel(varargin) + 1;
for q = 1:numel(varargin)
    if (ischar(varargin{q}) || isstring(varargin{q})) && any(strcmpi(varargin{q}, optNames))
        kopt = q; break;
    end
end
data = varargin(1:kopt-1);
opts = varargin(kopt:end);
fs = []; plots = false;
for q = 1:2:numel(opts)-1
    switch lower(char(opts{q}))
        case 'fs',    fs = opts{q+1};
        case 'plots', plots = logical(opts{q+1});
    end
end
if ~any(strcmpi(method, {'isovolume','lsq','qdc','fvc'}))
    error('rip_calibrate_K: method must be isovolume, lsq, qdc or fvc.');
end
if isempty(data)
    error('rip_calibrate_K: no calibration data supplied.');
end

cal = struct('K',NaN,'method',lower(string(method)),'ci',[NaN NaN], ...
             'quality',0,'note','');

switch lower(method)

%% ---------------------------------------------------------------- ISOVOLUME
case 'isovolume'
    tt = data{1};
    [rc, ab, fs] = getpair(tt, fs);
    % Keep the slow manoeuvre, drop drift and high-frequency noise
    [sos,g] = rip_bandpass(fs, 0.03, 2.0, 4);
    rc = filtfilt(sos,g,detrend(rc));
    ab = filtfilt(sos,g,detrend(ab));

    [K, r, ve] = tls_iso(rc, ab);
    cal.K = K; cal.r = r; cal.var_explained = ve;
    cal.ci = boot_iso(rc, ab, 500);
    % A clean isovolume line is a tight, strongly negative relationship
    cal.quality = max(0, min(1, (abs(r)-0.5)/0.45)) * max(0, min(1,(ve-0.75)/0.20));
    if r > -0.5
        cal.note = sprintf(['POOR: RC and AB are not moving in opposition (r = %+.2f). ' ...
            'The subject probably did not hold lung volume constant. Repeat the manoeuvre.'], r);
    elseif ve < 0.85
        cal.note = sprintf('MARGINAL: isovolume loop is not a tight line (%.0f%% of variance on one axis).', ve*100);
    else
        cal.note = sprintf('GOOD: r = %+.2f, %.0f%% of variance on the isovolume line.', r, ve*100);
    end

    if plots
        figure('Name','rip_calibrate_K - isovolume','NumberTitle','off','Color','w');
        plot(rc, ab, '.', 'Color',[0.4 0.4 0.4], 'MarkerSize',4); hold on
        xl = xlim; plot(xl, -xl/K, 'r-', 'LineWidth',2);
        xlabel('RC (au)'); ylabel('AB (au)'); axis equal; grid on; box off;
        title(sprintf('Isovolume calibration: K = %.3f [%.3f, %.3f]', cal.K, cal.ci(1), cal.ci(2)));
        legend({'samples','isovolume line'},'Box','off','Location','best');
    end

%% --------------------------------------------------------------------- LSQ
case 'lsq'
    dRC = data{1}(:); dAB = data{2}(:); VT = data{3}(:);
    ok = all(isfinite([dRC dAB VT]),2);
    X = [dRC(ok) dAB(ok)]; y = VT(ok);
    if size(X,1) < 8
        warning('rip_calibrate_K:fewBreaths','Only %d usable breaths for LSQ.', size(X,1));
    end
    ab_ = X\y;
    cal.a = ab_(1); cal.b = ab_(2); cal.K = ab_(2)/ab_(1);
    yhat = X*ab_;
    cal.R2 = 1 - sum((y-yhat).^2)/sum((y-mean(y)).^2);
    % conditioning: does the data actually span a range of RC/AB partitioning?
    partition = X(:,1)./max(X(:,1)+X(:,2),eps);
    cal.partition_range = range(partition);
    cal.ci = boot_lsq(X, y, 500);
    cal.quality = max(0,min(1,cal.R2)) * max(0,min(1, cal.partition_range/0.20));
    if cal.a <= 0 || cal.b <= 0
        cal.note = 'FAIL: a negative compartment gain means the regression is ill-conditioned. Widen the range of RC/AB partitioning during calibration.';
        cal.quality = 0;
    elseif cal.partition_range < 0.10
        cal.note = sprintf(['WEAK: RC share varied by only %.0f percentage points during calibration. ' ...
            'Instruct alternating chest-dominant and abdomen-dominant breathing.'], cal.partition_range*100);
    else
        cal.note = sprintf('GOOD: R2 = %.2f over a %.0f point range of RC share.', cal.R2, cal.partition_range*100);
    end

    if plots
        figure('Name','rip_calibrate_K - LSQ','NumberTitle','off','Color','w');
        plot(yhat, y, 'o','MarkerFaceColor',[0.2 0.5 0.85],'MarkerEdgeColor','none'); hold on
        lims = [min([yhat;y]) max([yhat;y])]; plot(lims,lims,'k--');
        xlabel('RIP-estimated VT'); ylabel('Reference VT'); grid on; box off; axis equal
        title(sprintf('LSQ calibration: K = %.3f, R^2 = %.2f', cal.K, cal.R2));
    end

%% --------------------------------------------------------------------- QDC
case 'qdc'
    dRC = data{1}(:); dAB = data{2}(:);
    ok = isfinite(dRC) & isfinite(dAB);
    dRC = dRC(ok); dAB = dAB(ok);
    cal.K = std(dRC)/std(dAB);
    cal.n_breaths = numel(dRC);
    cal.ci = boot_qdc(dRC, dAB, 500);
    cal.cv_rc = std(dRC)/mean(dRC); cal.cv_ab = std(dAB)/mean(dAB);
    cal.quality = max(0, min(1, (numel(dRC)-30)/70)) * 0.6;   % capped: QDC is assumption-heavy
    cal.note = sprintf(['QDC assumes the two calibrated compartments contribute equally to the ' ...
        'variance of VT during quiet breathing. Measured on %d breaths. Treat it as a fallback: ' ...
        'on a cycling pilot with an OEP reference, QDC gave a mean %%RC error of 6.3 points against OEP, versus ' ...
        '2.0 points for a well-anchored K.'], numel(dRC));

%% --------------------------------------------------------------------- FVC
case 'fvc'
    tt = data{1};
    [rc, ab, fs] = getpair(tt, fs); %#ok<ASGLU>
    cal.K = range(rc)/range(ab);
    cal.quality = 0.2;
    cal.note = ['NOT RECOMMENDED: rib cage and abdomen do not displace equal volumes during a ' ...
        'forced VC, so this ratio is not a gain ratio. Two FVC windows from one recording ' ...
        'gave K = 1.40 and K = 3.16. Use ''isovolume'' or ''lsq''.'];
    warning('rip_calibrate_K:fvc','%s', cal.note);
end

fprintf('\n===== rip_calibrate_K (%s) =====\n', upper(char(cal.method)));
fprintf('K = %.4f   95%% CI [%.4f, %.4f]   quality %.2f\n', cal.K, cal.ci(1), cal.ci(2), cal.quality);
fprintf('%s\n', cal.note);
fprintf('Sensitivity: %%RC moves about %.1f points per 10%% error in K at %%RC = 50.\n', 0.1*25);
fprintf('================================\n\n');
end

%% ------------------------------------------------------------------ LOCALS
function [rc, ab, fs] = getpair(tt, fs)
if ~istimetable(tt) || width(tt) < 2
    error('rip_calibrate_K: expected a 2-column timetable [thc abd].');
end
if isempty(fs), fs = tt.Properties.SampleRate; end
if isempty(fs) || isnan(fs)
    fs = 1/seconds(tt.Time(2)-tt.Time(1));
end
rc = tt{:,1}; ab = tt{:,2};
end

function [K, r, ve] = tls_iso(rc, ab)
% Total least squares: the isovolume line is the FIRST principal axis of the
% (RC, AB) cloud. Ordinary least squares would be biased because both axes
% carry measurement noise.
X = [rc(:) ab(:)];
X = X - mean(X,1);
[V, D] = eig(cov(X));
[d, o] = sort(diag(D),'descend');
V = V(:,o);
slope = V(2,1)/V(1,1);        % dAB/dRC along the principal axis
K = -1/slope;                  % 0 = dRC + K*dAB
r = corr(X(:,1), X(:,2));
ve = d(1)/sum(d);
end

function ci = boot_iso(rc, ab, nb)
n = numel(rc); Ks = nan(nb,1);
% block bootstrap (samples are serially correlated)
blk = max(10, round(n/20));
nblk = floor(n/blk);
for b = 1:nb
    idx = [];
    for q = 1:nblk
        s = randi(n-blk+1);
        idx = [idx, s:s+blk-1]; %#ok<AGROW>
    end
    try, Ks(b) = tls_iso(rc(idx), ab(idx)); catch, end
end
ci = prctile(Ks(isfinite(Ks)), [2.5 97.5]);
if isempty(ci), ci = [NaN NaN]; end
end

function ci = boot_lsq(X, y, nb)
n = size(X,1); Ks = nan(nb,1);
for b = 1:nb
    i = randi(n, n, 1);
    c = X(i,:)\y(i);
    Ks(b) = c(2)/c(1);
end
ci = prctile(Ks(isfinite(Ks)), [2.5 97.5]);
end

function ci = boot_qdc(a, b_, nb)
n = numel(a); Ks = nan(nb,1);
for b = 1:nb
    i = randi(n, n, 1);
    Ks(b) = std(a(i))/std(b_(i));
end
ci = prctile(Ks(isfinite(Ks)), [2.5 97.5]);
end
