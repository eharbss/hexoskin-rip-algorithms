function [q, qsum] = rip_quality(rc, ab, fs, varargin)
%RIP_QUALITY  Two-step validity gate for dual-band RIP metrics (%RC and TAC).
%
% Band-power SNR alone is not sufficient to decide whether %RC and TAC are
% meaningful: it asks "is there energy in the breathing band?", and baseline
% drift satisfies that. On a reference running recording the abdominal channel
% reached +24.6 dB breathing-band SNR after filtering while its dominant
% in-band peak was 0.156 Hz drift, not the 0.461 Hz breathing of the thoracic
% channel. The metrics computed from it were meaningless.
%
% This function therefore asks two different questions:
%
%   STEP 1 (window level) -- Do BOTH bands contain the SAME respiratory
%       oscillation?  Magnitude-squared coherence between RC and AB at the
%       breathing frequency tracked from the RC channel. Bias-corrected for
%       the finite number of Welch segments.
%
%   STEP 2 (breath level) -- Is the phase relationship stable enough to quote
%       a number?  Circular concentration R of the instantaneous RC-AB phase
%       difference within each breath (R = 1 rigid coupling, R = 0 noise).
%       Returned by enrich_BP as bp.tac_R; computed here too when breath
%       times are supplied.
%
% Syntax:
%   q            = rip_quality(rc, ab, fs)
%   [q, qsum]    = rip_quality(rc, ab, fs, 'fr1', fr1, 'plots', true)
%
% Inputs:
%   rc, ab - thoracic and abdominal signals: numeric vectors, or 1-column
%            timetables (raw or filtered; filtering is applied internally)
%   fs     - sample rate (Hz); inferred if both inputs are timetables
%
% Name-value parameters:
%   'fc_range'   - [lo hi] Hz passband (default [0.1 1.5])
%   'order'      - Butterworth order (default 4)
%   'win'        - analysis window (default seconds(30))
%   'step'       - window hop (default seconds(10))
%   'seg'        - Welch segment inside each window (default seconds(10))
%   'coh_thresh' - bias-corrected coherence required to pass (default 0.60).
%                  Chosen from a sweep over five recordings:
%                    threshold                     0.45   0.60   0.70
%                    Cycling, OEP pilot (clean)    100%   100%   100%
%                    Seated, voluntary patterns    97%    97%    97%
%                    Fast running intervals        99%    98%    96%
%                    Mountain ultra, 14.8 h        95%    92%    89%
%                    Known-bad running recording   34%    17%    10%
%                  0.60 keeps at least 92% of every good recording while
%                  rejecting 83% of the known-bad one.
%   'R_thresh'   - per-breath concentration required to pass (default 0.70)
%   'fr1'        - timetable with .insp column; enables the breath-level step
%   'acc'        - trunk accelerometer (timetable, or vector with 'acc_fs').
%                  Enables the cadence cross-check below. Strongly
%                  recommended for running data.
%   'acc_fs'     - accelerometer sample rate when 'acc' is a bare vector
%   'plots'      - diagnostic figure (default false)
%
% CADENCE CROSS-CHECK (why 'acc' matters)
%   Coherence between RC and AB is high wherever the two belts share a
%   common mechanical driver -- and footfall is exactly such a driver.
%   Measured on a fast-running recording, RC~AB coherence at the STEP
%   frequency was 0.97, as high as at the breathing frequency. The gate is
%   therefore only safe while the tracked f_b really is breathing. If
%   cadence (or a subharmonic of it) falls inside the 0.1-1.5 Hz passband,
%   as it does at walking and slow-running stride rates, the PSD peak can be
%   footfall and the window would pass on pure artifact. Supplying 'acc'
%   lets the function detect that collision and fail the window instead.
%
% Outputs:
%   q    - window timetable:
%            .f_b        tracked breathing frequency (Hz)
%            .coh        bias-corrected MS coherence RC~AB at f_b (0-1)
%            .coh_raw    uncorrected coherence
%            .bfrac_rc   fraction of RC energy in the breathing band
%            .bfrac_ab   fraction of AB energy in the breathing band
%            .snr_rc_dB  RC breathing-band vs step-band (1.5-4.5 Hz) power
%            .snr_ab_dB  AB breathing-band vs step-band power
%            .pass       step-1 verdict
%   qsum - struct: median values, percentage of windows passing, per-breath R
%          (when 'fr1' given), and a plain-language verdict.
%
% Typical values: coherence 0.85-0.99 in clean seated/cycling recordings;
% 0.1-0.4 in a running recording where the abdominal belt is artifact-dominated.
%
% created by Eric Harbour (PLUS)

%% Parse
p = inputParser;
addRequired (p,'rc');
addRequired (p,'ab');
addRequired (p,'fs',          @(x) isempty(x) || (isnumeric(x) && x>0));
addParameter(p,'fc_range',  [0.1 1.5], @(x) numel(x)==2 && x(1)<x(2));
addParameter(p,'order',     4,           @(x) isnumeric(x) && x>=2);
addParameter(p,'win',       seconds(30), @isduration);
addParameter(p,'step',      seconds(10), @isduration);
addParameter(p,'seg',       seconds(10), @isduration);
addParameter(p,'coh_thresh',0.60,        @(x) x>=0 && x<=1);
addParameter(p,'R_thresh',  0.70,        @(x) x>=0 && x<=1);
addParameter(p,'fr1',       [],          @(x) isempty(x) || istimetable(x));
addParameter(p,'acc',       [],          @(x) isempty(x) || istimetable(x) || isnumeric(x));
addParameter(p,'acc_fs',    [],          @(x) isempty(x) || (isnumeric(x) && x>0));
addParameter(p,'plots',     false,       @islogical);
parse(p, rc, ab, fs, varargin{:});
o = p.Results;

%% Ingest
t0 = [];
if istimetable(rc), t0 = rc.Time; if isempty(fs), fs = rc.Properties.SampleRate; end; rc = rc{:,1}; end
if istimetable(ab), if isempty(t0), t0 = ab.Time; end; ab = ab{:,1}; end
rc = rc(:); ab = ab(:);
if numel(rc) ~= numel(ab)
    error('rip_quality: rc and ab must be the same length.');
end
if isempty(t0), t0 = seconds((0:numel(rc)-1)'/fs); end
N = numel(rc);

% accelerometer for the cadence cross-check
acc = []; acc_fs = o.acc_fs;
if ~isempty(o.acc)
    if istimetable(o.acc)
        if isempty(acc_fs), acc_fs = o.acc.Properties.SampleRate; end
        av = o.acc.Properties.VariableNames;
        pick = find(ismember(av, {'res','rms','resultant','mag'}), 1);
        if isempty(pick)
            acc = vecnorm(o.acc{:,1:min(3,width(o.acc))}, 2, 2);
        else
            acc = o.acc{:,pick};
        end
    else
        acc = o.acc(:);
    end
    if isempty(acc_fs) || isnan(acc_fs)
        error('rip_quality: supply acc_fs with a bare accelerometer vector.');
    end
    acc = fillmissing(detrend(fillmissing(acc,'linear')),'constant',0);
    [az,ap,ak] = butter(4, [0.6 min(5, 0.45*acc_fs)]/(acc_fs/2), 'bandpass');
    acc = filtfilt(zp2sos(az,ap,ak), 1, acc);
end

[sos,g] = rip_bandpass(fs, o.fc_range(1), o.fc_range(2), o.order);
rcf = filtfilt(sos, g, detrend(rc));
abf = filtfilt(sos, g, detrend(ab));

fprintf('\n===== rip_quality =====\n');

%% Step 1 -- window-level coherence
win_n  = round(seconds(o.win)  * fs);
step_n = round(seconds(o.step) * fs);
seg_n  = min(round(seconds(o.seg) * fs), win_n);
nfft   = max(4096, 2^nextpow2(seg_n));
n_win  = max(1, floor((N - win_n)/step_n) + 1);

% number of independent-ish Welch averages -> coherence bias floor 1/nd
nd      = max(1, floor((win_n - seg_n)/(seg_n/2)) + 1);
coh_bias = 1/nd;

f_b = nan(n_win,1); coh = nan(n_win,1); coh_raw = nan(n_win,1);
f_stride = nan(n_win,1); f_conflict = false(n_win,1);

% breath onsets, when supplied, give a better f_b than the PSD
ev_t = [];
if ~isempty(o.fr1) && ismember('insp', o.fr1.Properties.VariableNames)
    ev_t = o.fr1.Time(o.fr1.insp == 1);
end
bf_rc = nan(n_win,1); bf_ab = nan(n_win,1);
sn_rc = nan(n_win,1); sn_ab = nan(n_win,1);
wt  = t0(1) + seconds((0:n_win-1)' * seconds(o.step));

step_band = [1.5 4.5];
hw = hamming(seg_n);

for w = 1:n_win
    i1 = (w-1)*step_n + 1;
    i2 = min(i1 + win_n - 1, N);
    if (i2-i1+1) < seg_n, continue; end
    % f_b and coherence come from the FILTERED signals and exclude the
    % lowest part of the band. Taken from the raw window instead, baseline
    % drift wins the PSD: three hours into a mountain ultra the tracked f_b
    % read 0.19 Hz against an actual breathing rate of 0.34 Hz, which then
    % measured coherence at the wrong frequency and failed good windows.
    a = detrend(rcf(i1:i2)); b = detrend(abf(i1:i2));
    araw = detrend(rc(i1:i2)); braw = detrend(ab(i1:i2));
    if std(a) < eps || std(b) < eps, continue; end

    [pa,f] = pwelch(a, hw, round(seg_n/2), nfft, fs);
    [pb,~] = pwelch(b, hw, round(seg_n/2), nfft, fs);

    f_search = max(o.fc_range(1)*1.2, 0.12);
    inb = f >= f_search & f <= o.fc_range(2);
    if ~any(inb), continue; end
    [~,im] = max(pa(inb)); fband = f(inb);
    f_b(w) = fband(im);
    % If breath events are available, prefer them: the detector has already
    % resolved the breathing rate and is not confused by the slow torso
    % movement that dominates the PSD late in a long trail run.
    if ~isempty(ev_t)
        inw = ev_t >= t0(i1) & ev_t <= t0(i2);
        if nnz(inw) >= 3
            f_ev = 1/seconds(median(diff(ev_t(inw))));
            if f_ev >= o.fc_range(1) && f_ev <= o.fc_range(2)
                f_b(w) = f_ev;
            end
        end
    end
    inb = f >= o.fc_range(1) & f <= o.fc_range(2);

    ins = f >= step_band(1) & f <= step_band(2);
    [par,~] = pwelch(araw, hw, round(seg_n/2), nfft, fs);
    [pbr,~] = pwelch(braw, hw, round(seg_n/2), nfft, fs);
    bf_rc(w) = trapz(f(inb),par(inb))/max(trapz(f,par),eps);
    bf_ab(w) = trapz(f(inb),pbr(inb))/max(trapz(f,pbr),eps);
    if any(ins)
        sn_rc(w) = 10*log10(trapz(f(inb),par(inb))/max(trapz(f(ins),par(ins)),eps));
        sn_ab(w) = 10*log10(trapz(f(inb),pbr(inb))/max(trapz(f(ins),pbr(ins)),eps));
    end

    [cxy,fc] = mscohere(a, b, hw, round(seg_n/2), nfft, fs);
    sel = fc >= f_b(w)*0.85 & fc <= f_b(w)*1.15;
    if ~any(sel), sel = abs(fc - f_b(w)) == min(abs(fc - f_b(w))); end
    coh_raw(w) = max(cxy(sel));
    coh(w)     = max(0, (coh_raw(w) - coh_bias)/(1 - coh_bias));

    % cadence cross-check: is the "breathing" peak actually footfall?
    if ~isempty(acc)
        j1 = max(1, round((i1-1)*acc_fs/fs) + 1);
        j2 = min(numel(acc), round(i2*acc_fs/fs));
        if (j2-j1+1) > round(acc_fs*4)
            [pq,fq] = pwelch(acc(j1:j2), hamming(min(round(acc_fs*8), j2-j1+1)), [], 4096, acc_fs);
            ina = fq >= 0.6 & fq <= 5;
            if any(ina)
                fqa = fq(ina); [~,ia] = max(pq(ina));
                f_stride(w) = fqa(ia);
                % footfall lands in band at the stride rate and its halves
                cand = f_stride(w) ./ [1 2 4];
                cand = cand(cand >= o.fc_range(1) & cand <= o.fc_range(2));
                if any(abs(cand - f_b(w)) ./ f_b(w) < 0.10)
                    f_conflict(w) = true;
                end
            end
        end
    end
end

pass = (coh >= o.coh_thresh) & ~f_conflict;
q = timetable(wt, f_b, coh, coh_raw, bf_rc, bf_ab, sn_rc, sn_ab, f_stride, double(f_conflict), double(pass), ...
    'VariableNames', {'f_b','coh','coh_raw','bfrac_rc','bfrac_ab','snr_rc_dB','snr_ab_dB', ...
                      'f_stride','f_conflict','pass'});
q.Properties.DimensionNames{1} = 'Time';

fprintf('Step 1  coherence RC~AB at f_b: median %.2f (bias floor %.2f, %d Welch averages)\n', ...
    median(coh,'omitnan'), coh_bias, nd);
fprintf('        windows passing (coh >= %.2f): %.0f%% of %d\n', ...
    o.coh_thresh, 100*mean(pass,'omitnan'), n_win);
fprintf('        median f_b %.2f Hz (%.0f br/min) | median SNR  RC %.1f dB  AB %.1f dB\n', ...
    median(f_b,'omitnan'), median(f_b,'omitnan')*60, ...
    median(sn_rc,'omitnan'), median(sn_ab,'omitnan'));
if ~isempty(acc)
    fprintf('        cadence cross-check: median stride %.2f Hz | %.0f%% of windows rejected for f_b/cadence collision\n', ...
        median(f_stride,'omitnan'), 100*mean(f_conflict));
elseif median(f_b,'omitnan') > 0.9
    fprintf('        NOTE: tracked f_b is high. Supply ''acc'' so cadence collisions can be ruled out.\n');
end

%% Step 2 -- per-breath phase concentration
R = []; 
if ~isempty(o.fr1) && ismember('insp', o.fr1.Properties.VariableNames)
    ti  = o.fr1.Time(o.fr1.insp == 1);
    dphi = angle(hilbert(rcf)) ; %#ok<NASGU>
    d = wrapToPi(unwrap(angle(hilbert(rcf))) - unwrap(angle(hilbert(abf))));
    nb = numel(ti) - 1;
    R = nan(nb,1);
    for k = 1:nb
        m = t0 >= ti(k) & t0 < ti(k+1);
        if nnz(m) < 5, continue; end
        R(k) = abs(mean(exp(1i*d(m)),'omitnan'));
    end
    fprintf('Step 2  per-breath phase concentration R: median %.2f | breaths passing (R >= %.2f): %.0f%%\n', ...
        median(R,'omitnan'), o.R_thresh, 100*mean(R >= o.R_thresh,'omitnan'));
end

%% Verdict
pct = 100*mean(pass,'omitnan');
if pct >= 80
    verdict = 'GOOD - %RC and TAC can be reported for this recording.';
elseif pct >= 40
    verdict = 'PARTIAL - gate per window; report %RC and TAC only where pass == 1.';
else
    verdict = 'FAIL - the two bands do not share a common respiratory oscillation. Do not report %RC or TAC.';
end
fprintf('VERDICT: %s\n', verdict);
fprintf('=======================\n\n');

qsum = struct('coh_median',median(coh,'omitnan'), 'coh_bias_floor',coh_bias, ...
    'pct_windows_pass',pct, 'f_b_median',median(f_b,'omitnan'), ...
    'snr_rc_dB_median',median(sn_rc,'omitnan'), 'snr_ab_dB_median',median(sn_ab,'omitnan'), ...
    'R_breath',R, 'R_median',median(R,'omitnan'), ...
    'coh_thresh',o.coh_thresh, 'R_thresh',o.R_thresh, ...
    'pct_f_conflict',100*mean(f_conflict), 'verdict',verdict);

%% Plot
if o.plots
    figure('Name','rip_quality','NumberTitle','off','Color','w');
    tl = tiledlayout(3,1,'TileSpacing','tight','Padding','tight');
    ax1 = nexttile; plot(t0, rcf, 'Color',[0.20 0.50 0.85]); hold on
    plot(t0, abf, 'Color',[0.85 0.30 0.10]); legend({'RC','AB'},'Box','off');
    ylabel('filtered (au)'); title('Band-limited signals'); box off; grid on; ax1.XColor='none';
    ax2 = nexttile; plot(q.Time, q.coh, 'k-', 'LineWidth',1.4); hold on
    yline(o.coh_thresh,'r--','threshold'); ylim([0 1]);
    ylabel('coherence at f_b'); title('Step 1: do both bands share one respiratory oscillation?');
    box off; grid on; ax2.XColor='none';
    ax3 = nexttile; plot(q.Time, q.snr_rc_dB,'-','Color',[0.20 0.50 0.85]); hold on
    plot(q.Time, q.snr_ab_dB,'-','Color',[0.85 0.30 0.10]); yline(0,'k:');
    ylabel('SNR (dB)'); xlabel('Time'); legend({'RC','AB'},'Box','off');
    title('Breathing-band vs step-band power'); box off; grid on;
    linkaxes([ax1 ax2 ax3],'x');
    title(tl, sprintf('rip\\_quality - %s', verdict), 'FontWeight','bold');
end
end
