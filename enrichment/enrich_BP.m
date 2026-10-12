function bp = enrich_BP(sig, fr1, varargin)
%ENRICH_BP  Breath-by-breath breathing pattern from segmented RIP data.
%
% updated 21-09-2026  (v2)
% created by Eric Harbour (PLUS)
%
% CHANGES IN v2
%   * TAC phase angle is now a CIRCULAR mean. v1 used mean(dphi_wrapped),
%     an arithmetic mean of angles, which collapses near antiphase: a
%     synthetic perfectly paradoxical signal read +48.9 deg instead of
%     +175.7 deg. One-line consequence, but it fires exactly in the
%     paradoxical case that matters clinically.
%   * bp.tac_R reports the per-breath circular concentration (0-1). This is
%     the confidence attached to every TAC value and drives the gate.
%   * Takes the sig timetable from segment_RIP so the thoracic and abdominal
%     channels are guaranteed to be identically filtered, and so the
%     THORACIC channel is used for %RC rather than the RC+AB fusion that
%     v1 received from the pipeline.
%   * The calibration constant K must now be supplied or derived from an
%     explicit calibration. v1 fell back to K = std(RC)/std(AB) of the raw
%     samples, which forces mean %RC to ~50% by construction (measured 48.8%
%     on a cycling pilot with an OEP reference). That fallback is removed; %RC returns NaN with
%     a warning instead. Use RIP_CALIBRATE_K.
%   * Optional validity gating from RIP_QUALITY.
%
% Syntax:
%   bp = enrich_BP(sig, fr1, 'fr2', fr2, 'K', 3.16)
%   bp = enrich_BP(sig, fr1, 'fr2', fr2, 'K', K, 'quality', q, 'plots', true)
%   bp = enrich_BP(sig_filt_1, fr1, 'sig_filt_2', sig_filt_2, ...)   % legacy
%
% Required inputs:
%   sig  - timetable from segment_RIP (.thc, .abd), or a 1-column filtered
%          thoracic timetable in legacy mode
%   fr1  - thoracic flow reversals (.insp, .exp)
%
% Name-value parameters:
%   'fr2'          - abdominal flow reversals from segment_ABD
%   'K'            - calibration constant, K = gain_AB / gain_RC, so that
%                    VT is proportional to dRC + K*dAB. Get it from
%                    RIP_CALIBRATE_K. Required for %RC.
%   'sig_filt_2'   - legacy: separate abdominal timetable
%   'baseline'     - timetable [thc abd] of a maximal manoeuvre, used only
%                    for depth normalisation (NOT for K -- see RIP_CALIBRATE_K)
%   'quality'      - window timetable from RIP_QUALITY; gates %RC and TAC
%   'R_thresh'     - min per-breath phase concentration for a valid TAC (default 0.70)
%   'rcp_mode'     - 'range'  : %RC from per-channel peak-to-trough (default)
%                    'signed' : %RC from signed inspiratory volume change.
%                               Physiologically correct when a compartment
%                               moves paradoxically (can exceed 0-100), but
%                               noisier because it depends on exact event
%                               timing. Validate before making it default.
%   'rcp_ref'      - reference %RC (scalar). Adds bp.rcp_d = rcp - rcp_ref.
%                    Within-subject CHANGE in %RC is far more robust to K
%                    error than the absolute level.
%   'rb_win'       - lookback window for reset breath detection (default minutes(1))
%   'rb_amp_thresh'- reset breath amplitude multiplier (default 2)
%   'rb_tb_thresh' - reset breath duration multiplier (default 2)
%   'plots'        - summary plots (default false)
%
% Output bp (row times = thoracic inspiration onsets):
%   Timing : .ti .te .tb .fr .fr_smo .dc
%   Depth  : .thc_bd_filt .thc_bd_fiqr .thc_bd_frms (+ abd_ and total bd_)
%   Derived: .rcp (%), .rcp_d, .lov, .drive
%   TAC    : .tac_ph (signed deg, + = RC leads), .tac_R (concentration 0-1)
%   Gate   : .gate  0 = ok, 1 = failed coherence, 2 = failed concentration
%            .rcp_g, .tac_ph_g  gated copies; the raw columns are never
%            overwritten, so the gate can be revisited without re-running
%
% NOTE ON tac_R AND PARADOXICAL BREATHING
%   tac_R is the circular CONCENTRATION of the within-breath phase
%   difference, so it measures how STABLE the phase relationship is, not how
%   large it is. Antiphase does not depress it: on a synthetic sweep R stayed
%   at 1.00 for every imposed phase difference from 0 to 180 deg, and on the
%   reverse-breathing recordings median R was 0.99 -- higher than in the
%   30-90 deg bins (0.96). Across 180 real breaths corr(|TAC|, R) was +0.19,
%   i.e. the more paradoxical breaths had slightly MORE stable phase.
%
%   What does depress R is one compartment becoming nearly immobile, because
%   the analytic phase of a near-zero signal is noise. R held above 0.99 down
%   to an abdominal excursion of 5% of thoracic and fell only below about 2%.
%   Every breath the R gate removed from a voluntary-pattern dataset had an abdominal /
%   thoracic amplitude ratio of 0.14 (median) against 0.26 overall. So the
%   gate preferentially removes breaths in the most extreme patterns -- which
%   is defensible, since the phase angle there genuinely is not well defined,
%   but it is worth checking bp.abd_bd_filt ./ bp.thc_bd_filt when a lot of
%   breaths are being gated.
%   Sighs  : .rb    0 normal, 1 sigh (big+long), 2 big only

%% Parse
p = inputParser;
addRequired (p, 'sig',            @istimetable);
addRequired (p, 'fr1',            @istimetable);
addParameter(p, 'fr2',            [],   @(x) isempty(x) || istimetable(x));
addParameter(p, 'K',              [],   @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
addParameter(p, 'sig_filt_2',     [],   @(x) isempty(x) || istimetable(x));
addParameter(p, 'baseline',       [],   @(x) isempty(x) || istimetable(x));
addParameter(p, 'quality',        [],   @(x) isempty(x) || istimetable(x));
addParameter(p, 'R_thresh',       0.70, @(x) x>=0 && x<=1);
addParameter(p, 'rcp_mode',       'range', @(x) any(strcmpi(x,{'range','signed'})));
addParameter(p, 'rcp_ref',        [],   @(x) isempty(x) || isscalar(x));
addParameter(p, 'rb_win',         minutes(1), @isduration);
addParameter(p, 'rb_amp_thresh',  2,    @(x) x>1);
addParameter(p, 'rb_tb_thresh',   2,    @(x) x>1);
addParameter(p, 'plots',          false,@islogical);
parse(p, sig, fr1, varargin{:});
o = p.Results;

%% Resolve channels
vn = sig.Properties.VariableNames;
if all(ismember({'thc','abd'}, vn))
    x1 = sig.thc; x2 = sig.abd; t_vec = sig.Time; has2 = true;
elseif ismember('thc', vn)
    x1 = sig.thc; t_vec = sig.Time; x2 = []; has2 = false;
else
    % legacy: sig is a single filtered thoracic column
    if width(sig) ~= 1
        error('enrich_BP: sig must come from segment_RIP, or be a 1-column thoracic timetable.');
    end
    x1 = sig{:,1}; t_vec = sig.Time; x2 = []; has2 = false;
end
if ~has2 && ~isempty(o.sig_filt_2)
    if width(o.sig_filt_2) ~= 1
        error('enrich_BP: sig_filt_2 must have exactly 1 column.');
    end
    x2 = o.sig_filt_2{:,1}; has2 = true;
    warning('enrich_BP:legacyInput', ...
        ['Channels supplied separately. Pass the sig timetable from segment_RIP ' ...
         'instead, so both channels are guaranteed to share one filter.']);
end
if ~ismember('insp', fr1.Properties.VariableNames) || ~ismember('exp', fr1.Properties.VariableNames)
    error('enrich_BP: fr1 must contain .insp and .exp columns.');
end
has_fr2  = ~isempty(o.fr2) && istimetable(o.fr2);
has_base = ~isempty(o.baseline) && istimetable(o.baseline);
if has_base
    o.baseline.Properties.VariableNames(1:2) = {'thc','abd'};
end

fprintf('\n===== enrich_BP (v2) =====\n');

%% Calibration constant
K = o.K;
if has2
    if isempty(K)
        warning('enrich_BP:noK', ...
            ['No K supplied -- %%RC will be NaN. Derive K with RIP_CALIBRATE_K ' ...
             '(isovolume manoeuvre preferred). The v1 fallback K = std(RC)/std(AB) ' ...
             'is NOT used because it forces mean %%RC to ~50%% by construction.']);
    else
        fprintf('Calibration constant K = %.4f  (VT ~ dRC + K*dAB)\n', K);
    end
end

%% Breath table
in_list  = fr1.Time(fr1.insp == 1);
exp_list = fr1.Time(fr1.exp  == 1);
n = numel(in_list) - 1;
if n < 1
    error('enrich_BP: fewer than two inspiration events.');
end
bp = timetable(in_list(1:n));
bp.Properties.DimensionNames{1} = 'Time';

%% Per-breath loop
for ii = 1:n
    next_exp  = min(exp_list(exp_list > in_list(ii)));
    next_insp = in_list(ii+1);
    if isempty(next_exp) || next_exp >= next_insp, continue; end

    bp.ti(ii) = next_exp  - in_list(ii);
    bp.te(ii) = next_insp - next_exp;
    bp.tb(ii) = next_insp - in_list(ii);
    bp.fr(ii) = 60 / seconds(bp.tb(ii));
    bp.dc(ii) = bp.ti(ii) / bp.tb(ii);

    m = t_vec >= in_list(ii) & t_vec < next_insp;

    bp.thc_bd_filt(ii) = range(x1(m));
    bp.thc_bd_fiqr(ii) = iqr(x1(m));
    bp.thc_bd_frms(ii) = rms(detrend(x1(m)));

    if has2
        bp.abd_bd_filt(ii) = range(x2(m));
        bp.abd_bd_fiqr(ii) = iqr(x2(m));
        bp.abd_bd_frms(ii) = rms(detrend(x2(m)));

        % signed inspiratory excursion, both channels on the SAME boundaries.
        % Hexoskin polarity is inverted: signal maximum = end-expiration,
        % so inspiratory excursion = value(start) - value(end of inspiration).
        i0 = find(m,1,'first');
        [~, rel] = min(abs(t_vec(m) - next_exp));
        i1 = i0 + rel - 1;
        bp.thc_dv(ii) = x1(i0) - x1(i1);
        bp.abd_dv(ii) = x2(i0) - x2(i1);

        bp.bd_filt(ii) = bp.thc_bd_filt(ii) + bp.abd_bd_filt(ii);
        bp.bd_fiqr(ii) = bp.thc_bd_fiqr(ii) + bp.abd_bd_fiqr(ii);
        bp.bd_frms(ii) = bp.thc_bd_frms(ii) + bp.abd_bd_frms(ii);

        bp.lov(ii)   = mean(x1(m) + x2(m));
        bp.drive(ii) = bp.bd_filt(ii) / max(seconds(bp.ti(ii)), eps);
    else
        bp.bd_filt(ii) = bp.thc_bd_filt(ii);
    end
end

%% Rib cage percentage
if has2
    if isempty(K)
        bp.rcp = nan(height(bp),1);
    elseif strcmpi(o.rcp_mode,'signed')
        bp.rcp = bp.thc_dv ./ (bp.thc_dv + K*bp.abd_dv) * 100;
    else
        bp.rcp = bp.thc_bd_filt ./ (bp.thc_bd_filt + K*bp.abd_bd_filt) * 100;
    end
    if ~isempty(o.rcp_ref)
        bp.rcp_d = bp.rcp - o.rcp_ref;
    end
end

%% Depth normalisation to a maximal manoeuvre
if has_base && has2
    base_bd_filt = range(o.baseline.thc) + range(o.baseline.abd);
    base_bd_fiqr = iqr(o.baseline.thc)   + iqr(o.baseline.abd);
    base_bd_frms = rms(detrend(o.baseline.thc)) + rms(detrend(o.baseline.abd));
    bp.thc_bd_filt_n = bp.thc_bd_filt / range(o.baseline.thc) * 100;
    bp.abd_bd_filt_n = bp.abd_bd_filt / range(o.baseline.abd) * 100;
    bp.bd_filt_n     = bp.bd_filt / base_bd_filt * 100;
    bp.bd_fiqr_n     = bp.bd_fiqr / base_bd_fiqr * 100;
    bp.bd_frms_n     = bp.bd_frms / base_bd_frms * 100;
    fprintf('Baseline bd (filt/fiqr/frms): %.4f / %.4f / %.4f\n', base_bd_filt, base_bd_fiqr, base_bd_frms);
end

%% Thoraco-abdominal coordination -- CIRCULAR statistics
if has2 && has_fr2
    [tac_ph, tac_R] = calc_tac_hilbert(x1, x2, t_vec, in_list(1:n+1));
    bp.tac_ph = tac_ph;
    bp.tac_R  = tac_R;
    fprintf('TAC: %d of %d breaths with a valid phase angle | median concentration R = %.2f\n', ...
        sum(~isnan(bp.tac_ph)), height(bp), median(bp.tac_R,'omitnan'));
end

%% Validity gate
bp.gate = zeros(height(bp),1);
if ~isempty(o.quality) && ismember('pass', o.quality.Properties.VariableNames)
    qt = t2num(o.quality.Time);
    bt = t2num(bp.Time);
    if max(qt) < min(bt) || min(qt) > max(bt)
        warning('enrich_BP:gateNoOverlap', ...
            ['The quality timetable does not overlap the breath times in time, so no ' ...
             'gating was applied. Pass rip_quality the sig TIMETABLE (e.g. sig(:,''thc'')), ' ...
             'not a bare vector, so its window times share the recording clock.']);
    elseif numel(qt) == 1
        bp.gate(o.quality.pass(1) < 0.5) = 1;
    else
        passv = interp1(qt, o.quality.pass, bt, 'previous', 'extrap');
        bp.gate(passv < 0.5) = 1;
    end
end
if ismember('tac_R', bp.Properties.VariableNames)
    bp.gate(bp.gate == 0 & bp.tac_R < o.R_thresh) = 2;
end
% The gate is ADVISORY. bp.rcp and bp.tac_ph are never overwritten, because
% throwing the values away means you cannot audit what was removed or revisit
% the threshold without re-running. Gated copies are provided alongside for
% drop-in averaging:
%     bp.rcp_g     = rcp    with failed-coherence breaths set to NaN
%     bp.tac_ph_g  = tac_ph with failed-coherence OR low-concentration breaths NaN
if ismember('rcp', bp.Properties.VariableNames)
    bp.rcp_g = bp.rcp;
    bp.rcp_g(bp.gate == 1) = NaN;
end
if ismember('tac_ph', bp.Properties.VariableNames)
    bp.tac_ph_g = bp.tac_ph;
    bp.tac_ph_g(bp.gate > 0) = NaN;
end
if any(bp.gate > 0)
    fprintf('Gate: %d breaths failed coherence, %d failed phase concentration (of %d)\n', ...
        sum(bp.gate==1), sum(bp.gate==2), height(bp));
    fprintf('      raw .rcp / .tac_ph are preserved; use .rcp_g / .tac_ph_g for gated means.\n');
end

%% Sigh detection
if has2
    bp.rb = find_rb(bp, o.rb_win, o.rb_amp_thresh, o.rb_tb_thresh);
    fprintf('Reset breaths: %d sigh (big+long) | %d big\n', sum(bp.rb==1), sum(bp.rb==2));
end

%% Smooth breathing rate (10 s rear-facing window)
fr_smo = smoothdata(bp(:,'fr'), 'movmean', [seconds(10) seconds(0)], 'omitnan');
bp.fr_smo = fr_smo.fr;

%% Clean
bp = standardizeMissing(bp, seconds(0));
bp(isnan(bp.fr) | bp.fr == 0, :) = [];
fprintf('Enriched: %d breaths\n', height(bp));
fprintf('==========================\n\n');

%% Plots
if o.plots
    plot_summary(bp, has2, has_fr2, x1, x2, t_vec);
end
end

%% ------------------------------------------------------------------ LOCALS
function [tac_ph, tac_R] = calc_tac_hilbert(x1, x2, t_vec, in_list)
%CALC_TAC_HILBERT  Signed thoraco-abdominal phase angle per breath.
%
% Instantaneous phase difference between the two band-limited channels,
% aggregated over each breath with a CIRCULAR mean:
%     phi = angle( mean( exp(1i*dphi) ) )
% and the concentration
%     R   = abs(  mean( exp(1i*dphi) ) )
% R = 1 means the two compartments hold a rigid phase relationship over the
% breath; R near 0 means the "phase angle" is not a meaningful quantity.
%
% Convention: positive = rib cage LEADS abdomen; 0 = synchronous;
% +/-180 = paradoxical. The full +/-180 range is available, unlike the
% asin(m/s) approximation which saturates at +/-90 (Hammer & Newth 2009).
%
% References:
%   Agostoni E, Mognoni P. J Appl Physiol. 1966;21(6):1827-1832.
%   Konno K, Mead J. J Appl Physiol. 1967;22(3):407-422.
%   Hammer J, Newth CJL. Paediatr Respir Rev. 2009;10(2):75-80.

phi1 = unwrap(angle(hilbert(x1)));
phi2 = unwrap(angle(hilbert(x2)));
d    = wrapToPi(phi1 - phi2);

nb = numel(in_list) - 1;
tac_ph = nan(nb,1); tac_R = nan(nb,1);
for k = 1:nb
    m = t_vec >= in_list(k) & t_vec < in_list(k+1);
    if nnz(m) < 5, continue; end
    zbar = mean(exp(1i*d(m)), 'omitnan');
    tac_ph(k) = angle(zbar) * 180/pi;
    tac_R(k)  = abs(zbar);
end
end

function rb = find_rb(bp, rb_win, amp_thresh, tb_thresh)
% Reset breaths (sighs) via a rear-facing window.
% Inspired by Vlemincx et al. (2013). Biol Psychol, 93(1), 24-32.
rb = zeros(height(bp),1);
first_idx = find(bp.Time > bp.Time(1) + rb_win, 1);
if isempty(first_idx)
    warning('find_rb: recording shorter than rb_win -- no reset breaths detected.');
    return
end
for ii = first_idx : height(bp)-1
    w = bp.Time >= bp.Time(ii) - rb_win & bp.Time < bp.Time(ii) - seconds(0.1);
    upper_bd = mean(bp.bd_filt(w),'omitnan') * amp_thresh;
    upper_tb = mean(seconds(bp.tb(w)),'omitnan') * tb_thresh;
    if isnan(upper_bd) || isnan(upper_tb), continue; end
    is_deep = bp.bd_filt(ii) >= upper_bd;
    is_long = seconds(bp.tb(ii)) >= upper_tb;
    if is_deep && is_long, rb(ii) = 1; elseif is_deep, rb(ii) = 2; end
end
end

function plot_summary(bp, has2, has_fr2, x1, x2, t_vec)
x_tm = minutes(bp.Time - bp.Time(1));
f = figure('Name','enrich_BP - breathing pattern summary','NumberTitle','off','Color','w');
tiledlayout('flow','TileSpacing','tight','Padding','tight');
ax = gobjects(6,1); ii = 0;

ii = ii+1; ax(ii) = nexttile;
plot(x_tm, bp.fr_smo, 'b-', 'LineWidth',1.5);
ylabel('f_R (bpm)'); grid on; box off; title('Breathing rate (smoothed)'); ax(ii).XColor='none';

ii = ii+1; ax(ii) = nexttile;
if has2
    sm = smoothdata(bp(:,'bd_filt'),'movmean',[seconds(10) seconds(0)],'omitnan');
    plot(x_tm, sm{:,1}, 'Color',[0 0.45 0.74],'LineWidth',1.2,'DisplayName','bd (total)'); hold on
    plot(x_tm, bp.thc_bd_filt, 'b:', 'DisplayName','thc');
    plot(x_tm, bp.abd_bd_filt, 'r:', 'DisplayName','abd'); hold off
    legend('Box','off','Location','best'); ylabel('Breath depth (au)');
else
    sm = smoothdata(bp(:,'thc_bd_filt'),'movmean',[seconds(10) seconds(0)],'omitnan');
    plot(x_tm, sm{:,1}, 'b-'); ylabel('Thoracic bd (au)');
end
grid on; box off; title('Breath depth'); ax(ii).XColor='none';

if ismember('rcp', bp.Properties.VariableNames) && any(isfinite(bp.rcp))
    ii = ii+1; ax(ii) = nexttile;
    sm = smoothdata(bp(:,'rcp'),'movmean',[seconds(10) seconds(0)],'omitnan');
    plot(x_tm, sm{:,1}, 'Color',[0.85 0.33 0.10],'LineWidth',1.2); hold on
    yline(50,'b--','HandleVisibility','off');
    text(min(x_tm), 45, 'AB dominant','FontSize',9,'Color','b','BackgroundColor','w');
    text(min(x_tm), 55, 'RC dominant','FontSize',9,'Color','b','BackgroundColor','w'); hold off
    ylabel('RC_P (%)'); ylim([0 100]); grid on; box off;
    title('Rib cage percentage'); ax(ii).XColor='none';
end

if has_fr2 && ismember('tac_ph', bp.Properties.VariableNames)
    ii = ii+1; ax(ii) = nexttile;
    sm = smoothdata(bp(:,'tac_ph'),'movmean',[seconds(10) seconds(0)],'omitnan');
    plot(x_tm, sm{:,1}, 'Color',[0.49 0.18 0.56],'LineWidth',1.2);
    yline(0,'b--','HandleVisibility','off');
    ylabel('\phi_{TAC} (\circ)'); grid on; box off;
    title('TAC phase angle (+ = rib cage leads)'); ax(ii).XColor='none';

    ii = ii+1; ax(ii) = nexttile;
    plot(x_tm, bp.tac_R, '.', 'Color',[0.35 0.35 0.35]); hold on
    yline(0.7,'r--','HandleVisibility','off'); hold off
    ylabel('R'); ylim([0 1]); grid on; box off;
    title('Per-breath phase concentration (confidence in TAC)');
end

xlabel('Time (min)');
linkaxes(ax(1:ii),'x');

% Lissajous
if has2
    figure('Name','enrich_BP - Lissajous','NumberTitle','off','Color','w');
    rc_n = rescale(x1,-1,1); ab_n = rescale(x2,-1,1);
    scatter(rc_n, ab_n, 1, (1:numel(rc_n))'/numel(rc_n), 'filled', ...
        'MarkerFaceAlpha',0.3,'HandleVisibility','off');
    colormap(gca,'cool'); cb = colorbar; cb.Label.String = 'Elapsed time (normalised)';
    xlabel('RC (normalised)'); ylabel('AB (normalised)');
    title('Lissajous: RC vs AB'); axis equal; grid on; box off; hold on
    plot([-1 1],[-1  1],'b--','LineWidth',1.2,'DisplayName','\theta = 0\circ (sync)');
    plot([-1 1],[ 1 -1],'r--','LineWidth',1.2,'DisplayName','\theta = 180\circ (paradox)');
    if ismember('tac_ph', bp.Properties.VariableNames)
        phm = angle(mean(exp(1i*deg2rad(bp.tac_ph(~isnan(bp.tac_ph))))))*180/pi;
        text(0.02,0.95, sprintf('\\phi_{circ mean} = %.1f\\circ', phm), ...
            'Units','normalized','FontSize',11,'Color',[0.49 0.18 0.56]);
    end
    legend('Location','best','Box','off'); hold off
end
sgtitle(f,'Breathing pattern summary');
end

function s = t2num(t)
% Common numeric time axis for duration or datetime row times.
if isdatetime(t)
    s = seconds(t - datetime(2000,1,1,'TimeZone',t.TimeZone));
else
    s = seconds(t);
end
end
