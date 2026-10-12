function [fr, sig, snr] = segment_RIP(rip, varargin)
%SEGMENT_RIP  Segment RIP data (e.g. Hexoskin) into flow-reversal moments.
%
% updated 21-09-2026  (v2)
% created by Eric Harbour (PLUS)
% cite as: Harbour, E., Lasshofer, M., Genitrini, M., & Schwameder, H. (2021).
%   Enhanced breathing pattern detection during running using wearable
%   sensors. Sensors, 21(16), 5606.
%
% CHANGES IN v2
%   * One filter for everything. The passband is obtained from RIP_BANDPASS
%     and applied IDENTICALLY to the thoracic and abdominal channels. In v1
%     this function used an adaptive upper cutoff while segment_ABD used a
%     different one, which biased the thoraco-abdominal phase angle.
%   * 4th-order Butterworth in second-order-section form (v1: 2nd-order,
%     transfer-function form). On a reference running recording this raised
%     abdominal breathing-band SNR from +1.6 dB to +24.6 dB.
%   * sig now returns the two channels SEPARATELY as well as the fusion.
%     v1 returned only the fusion, and the pipeline then passed it to
%     enrich_BP as if it were the thoracic channel. Detection still runs on
%     the fusion (breath timing should reference the whole chest wall), but
%     %RC and TAC must be computed from the separate bands.
%   * The fusion weights come from breathing-band SNR, not total in-band
%     variance. v1 weighted by variance over 0.1-2 Hz, which includes step
%     energy, so a noisier abdominal channel earned a HIGHER weight.
%
% Syntax:
%   [fr,sig,snr] = segment_RIP(rip)
%   [fr,sig,snr] = segment_RIP(rip, 'fc_range',[0.1 1.5], 'order',4, ...)
%
% Inputs:
%   rip  - Regular timetable with 1 column (RC) or 2 columns (RC, AB).
%
% Name-value parameters:
%   'fc_range'    - [low high] Hz passband (default [0.1 1.5]; 1.5 Hz = 90 br/min)
%   'order'       - Butterworth order (default 4)
%   'adaptive_hi' - shrink the upper cutoff to the tracked breathing
%                   frequency (default false). Improves SNR during running
%                   but will clip genuinely fast breathing, so it is off.
%   'tb_thresh'   - fraction of est. Tb for the short-breath gate (default 0.6)
%   'amp_thresh'  - min abs(signal) fraction for the ZC amplitude gate (default 0.10)
%   'debug_mode'  - print debugging information (default false)
%   'plots'       - show the detection plot (default false)
%
% Outputs:
%   fr  - timetable: .insp (1 = start of inspiration), .exp (1 = start of expiration)
%         NOTE on polarity: the Hexoskin RIP signal is INVERTED with respect
%         to lung volume -- a signal maximum is end-expiration. The detector
%         places .insp at maxima, which is therefore correct. This is
%         verified at runtime and reported (see sig.Properties.CustomProperties).
%   sig - timetable:
%           .thc_raw .abd_raw  unfiltered inputs
%           .thc     .abd      band-limited, IDENTICAL filter on both
%           .fused             SNR-weighted combination (detection only)
%           .det               fused + Savitzky-Golay, the detector input
%   snr - struct: .f_breath_hz .f_step_hz .SNR_dB .breath_frac .step_frac
%                 .w_thc .w_abd .win_tt (30 s windows, 10 s hop)

%% Parse inputs
p = inputParser;
addRequired (p, 'rip',         @istimetable);
addParameter(p, 'fc_range',    [0.1 1.5], @(x) numel(x)==2 && x(1)<x(2));
addParameter(p, 'order',       4,         @(x) isnumeric(x) && x>=2);
addParameter(p, 'adaptive_hi', false,     @islogical);
addParameter(p, 'tb_thresh',   0.6,       @(x) x>0.1 && x<0.9);
addParameter(p, 'amp_thresh',  0.10,      @(x) x>0 && x<1);
addParameter(p, 'debug_mode',  false,     @islogical);
addParameter(p, 'plots',       false,     @islogical);
parse(p, rip, varargin{:});

fc_range   = p.Results.fc_range;
ford       = p.Results.order;
adaptive   = p.Results.adaptive_hi;
tb_thresh  = p.Results.tb_thresh;
amp_thr    = p.Results.amp_thresh;
debug_mode = p.Results.debug_mode;
do_plot    = p.Results.plots;

%% Validate
if ~isregular(rip) || ~issorted(rip) || any(any(ismissing(rip)))
    error('segment_RIP: rip must be a regular, sorted timetable with no missing values.');
end
n_chan = width(rip);
if n_chan > 2
    error('segment_RIP: rip may have at most 2 columns (RC, AB).');
end
rip.Properties.DimensionNames{1} = 'Time';
fs = rip.Properties.SampleRate;

fprintf('\n===== segment_RIP (v2) =====\n');

%% Step 1: filter both channels with ONE design
x_thc = rip{:,1};
x_abd = [];
if n_chan == 2, x_abd = rip{:,2}; end

% Optional adaptive upper cutoff, estimated from the thoracic channel
fc_hi = fc_range(2);
if adaptive
    f_pre = estimate_br_psd(x_thc, fs, max(fc_range(1)*1.2,0.12), fc_range(2));
    if isfinite(f_pre), fc_hi = min(fc_range(2), max(0.5, 2.5*f_pre)); end
end
[sos, gg, finfo] = rip_bandpass(fs, fc_range(1), fc_hi, ford);

y_thc = filtfilt(sos, gg, x_thc);
if n_chan == 2, y_abd = filtfilt(sos, gg, x_abd); else, y_abd = []; end

% Estimate the dominant breathing frequency from the FILTERED signal and
% exclude the bottom of the band. Run on the raw signal, baseline drift wins
% the PSD on long recordings: three hours into a mountain ultra this read
% 0.12 Hz (8 br/min) against an actual 20 br/min, which then mis-set the
% state machine's expected breath period.
[f_dom, ~, psd_info] = estimate_br_psd(y_thc, fs, max(fc_range(1)*1.2,0.12), fc_hi); %#ok<ASGLU>

est_tb = 1 / max(f_dom, eps);
fprintf('Filter: %.2f-%.2f Hz, order %d, SOS (max pole radius %.4f)\n', ...
    finfo.fc_lo, finfo.fc_hi, finfo.order, finfo.max_pole_radius);
fprintf('Dominant breathing frequency: %.2f Hz (%.0f br/min)\n', f_dom, f_dom*60);

%% Step 2: SNR-weighted fusion (detection signal only)
% Each channel is first normalised to unit in-band SD, so the weighting
% reflects signal QUALITY rather than belt gain, then weighted by the ratio
% of breathing-band to out-of-band power.
if n_chan == 2
    r1 = band_quality(x_thc, fs, fc_range(1), fc_hi);
    r2 = band_quality(x_abd, fs, fc_range(1), fc_hi);
    w1 = r1 / max(r1 + r2, eps);
    w2 = 1 - w1;
    s1 = std(y_thc); s2 = std(y_abd);
    fused = w1 * (y_thc / max(s1,eps)) + w2 * (y_abd / max(s2,eps));
    fused = fused * s1;                       % back to roughly thoracic units
    fprintf('Fusion by breathing-band SNR: %s w=%.2f | %s w=%.2f\n', ...
        rip.Properties.VariableNames{1}, w1, rip.Properties.VariableNames{2}, w2);
else
    w1 = 1; w2 = 0; fused = y_thc;
    fprintf('Single sensor: %s\n', rip.Properties.VariableNames{1});
end

% Savitzky-Golay smoother for detection only (preserves peak timing).
sg_frame = 2*round(0.10 * est_tb * fs) + 1;
sg_frame = max(sg_frame, 5);
det = sgolayfilt(fused, 3, sg_frame);

%% Step 3: SNR estimation (breathing vs locomotor artifact)
step_range = [1.5 4.5];
bw_step    = 0.4;
snr = struct();
try
    [pxx_r, f_r] = pwelch(x_thc, hamming(round(fs*4)), [], [], fs);
    in_br    = f_r >= fc_range(1) & f_r <= fc_hi;
    in_step  = f_r >= step_range(1) & f_r <= step_range(2);
    in_total = f_r >= fc_range(1);
    [~, sp_idx] = max(pxx_r(in_step));
    f_step_all  = f_r(in_step);
    f_step      = f_step_all(sp_idx);
    in_step_band = f_r >= (f_step - bw_step) & f_r <= (f_step + bw_step);

    P_breath = trapz(f_r(in_br),        pxx_r(in_br));
    P_step   = trapz(f_r(in_step_band), pxx_r(in_step_band));
    P_total  = trapz(f_r(in_total),     pxx_r(in_total));

    snr.f_breath_hz = f_dom;
    snr.f_step_hz   = f_step;
    snr.SNR_dB      = 10*log10(P_breath / max(P_step, eps));
    snr.breath_frac = P_breath / max(P_total, eps);
    snr.step_frac   = P_step   / max(P_total, eps);

    fprintf('SNR: %.1f dB | breath %.1f%% vs step %.2f Hz (%.1f%%) of total energy\n', ...
        snr.SNR_dB, snr.breath_frac*100, f_step, snr.step_frac*100);
    if snr.SNR_dB < -6
        warning('segment_RIP: low breathing SNR (%.1f dB). Run rip_quality before trusting %%RC or TAC.', snr.SNR_dB);
    end
catch ME
    warning('segment_RIP: could not calculate SNR features (%s).', ME.message);
    f_step = NaN;
end
snr.w_thc = w1; snr.w_abd = w2;

% Rolling-window SNR (30 s window, 10 s hop)
win_s = 30; step_s = 10;
win_n = round(win_s*fs); step_n = round(step_s*fs);
n_wins = max(1, floor((height(rip) - win_n)/step_n) + 1);
win_SNR = NaN(n_wins,1); win_bfrac = NaN(n_wins,1);
win_times = rip.Time(1) + seconds((0:n_wins-1)' * step_s);
wlen = min(round(fs*8), win_n);
nfft = max(4096, 2^nextpow2(wlen));
hw   = hamming(wlen);
try
    for ww = 1:n_wins
        i1 = (ww-1)*step_n + 1;
        i2 = min(i1 + win_n - 1, height(rip));
        seg = detrend(x_thc(i1:i2), 'constant');
        if numel(seg) < round(fs*4), continue; end
        if numel(seg) >= wlen, wu = hw; else, wu = hamming(numel(seg)); end
        [pw, fw] = pwelch(seg, wu, [], nfft, fs);
        fw = fw(:); pw = pw(:);
        pb = safe_bandpower(fw, pw, fc_range(1), fc_hi);
        ps = safe_bandpower(fw, pw, f_step - bw_step, f_step + bw_step);
        pt = safe_bandpower(fw, pw, fc_range(1), max(fw));
        if isfinite(pb) && isfinite(ps) && isfinite(pt) && pb > 0 && pt > 0
            win_SNR(ww)   = 10*log10(pb / max(ps, eps));
            win_bfrac(ww) = pb / pt;
        end
    end
catch ME
    warning('segment_RIP: could not calculate all SNR features (%s).', ME.message);
end
snr.win_tt = timetable(win_times, win_SNR, win_bfrac, 'VariableNames', {'SNR_dB','breath_frac'});

%% Step 4: flow-reversal detection (state machine on the smoothed fusion)
min_tb_possible = 1 / fc_range(2);
buf_s     = min(max(0.25, 0.35*est_tb), 0.40*min_tb_possible);
min_gap_s = min(est_tb*tb_thresh, 0.5*min_tb_possible);
min_gap   = round(min_gap_s * fs);

fr = breath_statemachine(det, rip.Time, fs, est_tb, amp_thr, tb_thresh, buf_s, min_gap, debug_mode);

%% Step 5-6: cleanup
det_tt = timetable(rip.Time, det, 'VariableNames', {'det'});
fr = local_clean(fr);
fr = recover_missed(fr, det_tt, rip.Time, est_tb);
fr = local_clean(fr);

fprintf('Detected: %d inspirations | %d expirations\n', sum(fr.insp==1), sum(fr.exp==1));

%% Step 7: assemble the signal timetable
if n_chan == 2
    sig = timetable(rip.Time, x_thc, x_abd, y_thc, y_abd, fused, det, ...
        'VariableNames', {'thc_raw','abd_raw','thc','abd','fused','det'});
else
    sig = timetable(rip.Time, x_thc, y_thc, fused, det, ...
        'VariableNames', {'thc_raw','thc','fused','det'});
end
sig.Properties.DimensionNames{1} = 'Time';
sig = addprop(sig, {'filter_info','fusion_w','polarity_note'}, {'table','table','table'});
sig.Properties.CustomProperties.filter_info = finfo;
sig.Properties.CustomProperties.fusion_w    = [w1 w2];
sig.Properties.CustomProperties.polarity_note = ...
    'Hexoskin RIP is inverted w.r.t. lung volume: signal maximum = end-expiration, so .insp markers sit at maxima.';

fprintf('============================\n\n');

%% Step 8: plot
if do_plot
    plot_fr_detection(sig, fr);
end
end

%% ------------------------------------------------------------------ LOCALS
function r = band_quality(x, fs, lo, hi)
% Ratio of breathing-band power to everything else in the signal.
[pxx,f] = pwelch(detrend(x), hamming(round(fs*8)), [], 4096, fs);
inb = f >= lo & f <= hi;
pb  = trapz(f(inb), pxx(inb));
pt  = trapz(f, pxx);
r   = pb / max(pt - pb, eps);
end

function [max_fc, fc_all, info] = estimate_br_psd(rip_raw, fs, fc_lo, fc_hi)
% Dominant respiratory frequency from the PSD. Whole signal if < 5 min,
% otherwise 2-min non-overlapping windows; returns the 90th percentile.
x = rip_raw(:); n = numel(x); dur_s = n/fs;
seg_len_s = 120; seg_len_n = round(seg_len_s*fs);
welch_len = max(round(8*fs), 256);
nfft = max(4096, 2^nextpow2(welch_len));
fc_all = []; seg_idx = [];
if dur_s < 300
    segments = [1 n];
else
    n_seg = floor(n/seg_len_n);
    segments = zeros(n_seg,2);
    for k = 1:n_seg
        segments(k,:) = [(k-1)*seg_len_n+1, k*seg_len_n];
    end
end
for k = 1:size(segments,1)
    xw = detrend(x(segments(k,1):segments(k,2)), 'constant');
    if numel(xw) < welch_len, continue; end
    [pxx,f] = pwelch(xw, hamming(welch_len), [], nfft, fs);
    in_band = (f >= fc_lo) & (f <= fc_hi);
    if ~any(in_band), continue; end
    [~,idx] = max(pxx(in_band));
    fb = f(in_band);
    fc_all(end+1,1) = fb(idx);   %#ok<AGROW>
    seg_idx(end+1,1) = k;        %#ok<AGROW>
end
if isempty(fc_all), max_fc = NaN; else, max_fc = prctile(fc_all, 90); end
info = struct('fc_lo',fc_lo,'fc_hi',fc_hi,'seg_len_s',seg_len_s, ...
    'welch_len',welch_len,'nfft',nfft,'segments_used',seg_idx,'duration_s',dur_s);
end

function val = safe_bandpower(f, p, f1, f2)
f = f(:); p = p(:);
if ~isfinite(f1) || ~isfinite(f2), val = NaN; return; end
lo = max(min(f), min(f1,f2)); hi = min(max(f), max(f1,f2));
if hi <= lo, val = NaN; return; end
mask = (f >= lo) & (f <= hi);
if nnz(mask) < 2, val = NaN; return; end
val = trapz(f(mask), p(mask));
end

function fr = breath_statemachine(sig, t, fs, est_tb, amp_thr, tb_thresh, buf_s, min_gap, debug_mode) %#ok<INUSD>
% History-buffer state machine with adaptive amplitude thresholds.
% At each sample the oldest point in a trailing buffer is checked: if it is
% the buffer maximum (insp) or minimum (exp) it is a confirmed extremum.
% State alternates strictly INSP -> EXP -> INSP.
sig = sig(:); n = numel(sig);
buf_n = max(3, round(buf_s*fs));
max_events = ceil(n/max(min_gap,1)) + 2;
insp_idx = zeros(max_events,1); exp_idx = zeros(max_events,1);
n_insp = 0; n_exp = 0; last_ev_oi = 1;

max_fwd = movmax(sig, [0 buf_n]);
min_fwd = movmin(sig, [0 buf_n]);
max_fwd(max(1,n-buf_n+1):n) =  Inf;
min_fwd(max(1,n-buf_n+1):n) = -Inf;
is_peak   = (sig >= max_fwd);
is_trough = (sig <= min_fwd);

adapt_n = round(5*est_tb*fs);
hi = movmax(sig, [adapt_n adapt_n]);   % cheap robust local range
lo = movmin(sig, [adapt_n adapt_n]);
local_range = hi - lo;

boot_end = min(n, round(3*est_tb*fs) + buf_n);
state = 1; last_val = sig(min(buf_n+1,n));
for bi = buf_n+1 : boot_end
    oi = bi - buf_n;
    thr = amp_thr * max(local_range(oi), eps);
    if is_peak(oi) && abs(sig(oi)-last_val) > thr
        last_val = sig(oi); state = 2; break;
    elseif is_trough(oi) && abs(sig(oi)-last_val) > thr
        last_val = sig(oi); state = 1; break;
    end
end

i = buf_n + 1;
while i <= n
    oi = i - buf_n;
    if (oi - last_ev_oi) < min_gap, i = i + 1; continue; end
    thr = amp_thr * max(local_range(oi), eps);
    if state == 1
        if is_peak(oi) && (sig(oi) - last_val) > thr
            n_insp = n_insp + 1; insp_idx(n_insp) = oi;
            last_val = sig(oi); last_ev_oi = oi; state = 2;
            i = oi + buf_n + max(1, round(min_gap*0.5)); continue;
        end
    else
        if is_trough(oi) && (last_val - sig(oi)) > thr
            n_exp = n_exp + 1; exp_idx(n_exp) = oi;
            last_val = sig(oi); last_ev_oi = oi; state = 1;
            i = oi + buf_n + max(1, round(min_gap*0.5)); continue;
        end
    end
    i = i + 1;
end

insp_idx = insp_idx(1:n_insp); exp_idx = exp_idx(1:n_exp);
all_times = [t(insp_idx); t(exp_idx)];
all_insp  = [ones(n_insp,1);  zeros(n_exp,1)];
all_exp   = [zeros(n_insp,1); ones(n_exp,1)];
[all_times, si] = sort(all_times);
all_insp = all_insp(si); all_exp = all_exp(si);
if isempty(all_times)
    fr = timetable(datetime.empty(0,1), zeros(0,1), zeros(0,1), 'VariableNames', {'insp','exp'});
else
    fr = timetable(all_times, all_insp, all_exp, 'VariableNames', {'insp','exp'});
end
fr.Properties.DimensionNames{1} = 'Time';
end

function fr = local_clean(fr)
fr = fillmissing(fr, 'constant', 0);
fr = fr(any(fr.Variables, 2), :);
keep = true(height(fr),1);
for k = 1:height(fr)-1
    if     fr.insp(k)==1 && fr.insp(k+1)==1, keep(k) = false;
    elseif fr.exp(k)==1  && fr.exp(k+1)==1,  keep(k) = false;
    end
end
fr = fr(keep,:);
end

function fr = recover_missed(fr, filt, t_vec, est_tb)
% Insert missed flow reversals in gaps exceeding 1.5 x est_tb.
MAX_PASSES = 5;
zci_fn = @(v) find(diff(sign(v(:))));
min_gap_s = est_tb * 1.5;
for pass = 1:MAX_PASSES
    new_rows = {}; n_before = height(fr); k = 1;
    while k < height(fr)
        gap_s = seconds(fr.Time(k+1) - fr.Time(k));
        if gap_s <= min_gap_s, k = k + 1; continue; end
        n_missing = max(1, round(gap_s/est_tb) - 1);
        if fr.insp(k) == 1, first_type = 'exp'; else, first_type = 'insp'; end
        mask = t_vec > fr.Time(k) & t_vec < fr.Time(k+1);
        gap_sig = filt{mask,1}; gap_t = t_vec(mask);
        if numel(gap_sig) < 5, k = k + 1; continue; end
        d = gradient(gap_sig); d2 = gradient(d);
        zc = zci_fn(d); zc = zc(zc+1 <= numel(d));
        zc_peak   = zc(d2(zc) < 0);
        zc_trough = zc(d2(zc) > 0);
        n_gap = numel(gap_sig);
        fp = (zc_peak-1)/max(n_gap-1,1); ft = (zc_trough-1)/max(n_gap-1,1);
        zc_peak   = zc_peak(fp > 0.15 & fp < 0.85);
        zc_trough = zc_trough(ft > 0.15 & ft < 0.85);
        cur_type = first_type; last_frac = 0;
        for m = 1:n_missing
            if strcmp(cur_type,'insp'), cands = zc_peak; else, cands = zc_trough; end
            frac_c = (cands-1)/max(n_gap-1,1);
            cands = cands(frac_c > last_frac + 0.05);
            if isempty(cands), break; end
            [~,bi] = max(abs(gap_sig(cands)));
            best_idx = cands(bi);
            last_frac = (best_idx-1)/max(n_gap-1,1);
            new_rows{end+1} = timetable(gap_t(best_idx), ...
                double(strcmp(cur_type,'insp')), double(strcmp(cur_type,'exp')), ...
                'VariableNames', {'insp','exp'}); %#ok<AGROW>
            if strcmp(cur_type,'insp'), cur_type = 'exp'; else, cur_type = 'insp'; end
        end
        k = k + 1;
    end
    if ~isempty(new_rows)
        fr = sortrows([fr; vertcat(new_rows{:})]);
        fr = local_clean(fr);
    end
    if height(fr) == n_before, break; end
end
end

function plot_fr_detection(sig, fr)
figure('Name','segment_RIP','NumberTitle','off','Color','w');
plot(sig.Time, sig.thc, '-', 'Color',[0.20 0.50 0.85], 'DisplayName','RC (filtered)');
hold on
if ismember('abd', sig.Properties.VariableNames)
    plot(sig.Time, sig.abd, '-', 'Color',[0.85 0.30 0.10], 'DisplayName','AB (filtered)');
end
plot(sig.Time, sig.det, '-', 'Color',[0.35 0.35 0.35], 'LineWidth',1.1, 'DisplayName','detection signal');
ti = fr.Time(fr.insp==1); te = fr.Time(fr.exp==1);
plot(ti, sig.det(ti), 'b^', 'MarkerSize',7, 'DisplayName','start inspiration');
plot(te, sig.det(te), 'rv', 'MarkerSize',7, 'DisplayName','start expiration');
hold off
legend('Location','best','Box','off'); grid on; box off;
ylabel('Signal (au)'); xlabel('Time');
title('segment\_RIP v2 - detected flow reversals');
end
