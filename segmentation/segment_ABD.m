function [fr2, info] = segment_ABD(sig, fr1, varargin)
%SEGMENT_ABD  Guided abdominal breath segmentation from thoracic events.
%
% updated 21-09-2026  (v2)
% created by Eric Harbour (PLUS)
% cite as: Harbour, E., Lasshofer, M., Genitrini, M., & Schwameder, H. (2021).
%   Sensors, 21(16), 5606.
%
% CHANGES IN v2
%   * Takes the sig timetable from segment_RIP instead of a raw abdominal
%     channel, so the abdominal signal is guaranteed to have passed through
%     exactly the same filter as the thoracic one. v1 designed its own
%     passband ([0.1, min(1.5, est_BR*2.5)]) which differed from the one
%     segment_RIP used, biasing the thoraco-abdominal phase angle.
%   * SNR bookkeeping moved to RIP_QUALITY, which answers the question that
%     actually matters (do both bands share one respiratory oscillation?).
%
% Uses the detected thoracic flow reversals (fr1) to define half-cycle
% search windows for the abdominal signal, then finds exactly one abdominal
% event (peak or trough) per window.
%
% Syntax:
%   fr2          = segment_ABD(sig, fr1)
%   [fr2, info]  = segment_ABD(sig, fr1, 'plots', true)
%
% Inputs:
%   sig  - timetable from segment_RIP; must contain .thc and .abd
%          (a 2-column timetable [thc abd] is also accepted)
%   fr1  - timetable from segment_RIP with .insp and .exp
%
% Name-value parameters:
%   'amp_thresh'   - min window amplitude as a fraction of the global IQR
%                    below which the window is flagged low quality (default 0.15)
%   'phase_margin' - extends the trailing edge of each search window (default 0.30)
%   'plots'        - show the comparison plot (default false)
%
% Outputs:
%   fr2  - timetable: .insp (1 = abdominal end of the inspiratory window),
%                     .exp  (1 = abdominal end of the expiratory window),
%                     .quality (1 = high-amplitude window, found inside the
%                     nominal window rather than in the margin)
%   info - struct with window counts and the high-quality percentage

%% Parse
p = inputParser;
addRequired (p, 'sig',          @istimetable);
addRequired (p, 'fr1',          @istimetable);
addParameter(p, 'amp_thresh',   0.15, @(x) x>0 && x<1);
addParameter(p, 'phase_margin', 0.30, @(x) x>=0 && x<1);
addParameter(p, 'plots',        false, @islogical);
parse(p, sig, fr1, varargin{:});
amp_thr      = p.Results.amp_thresh;
phase_margin = p.Results.phase_margin;
do_plot      = p.Results.plots;

%% Validate / resolve channels
vn = sig.Properties.VariableNames;
if all(ismember({'thc','abd'}, vn))
    thc = sig.thc; abd = sig.abd;
elseif width(sig) == 2
    thc = sig{:,1}; abd = sig{:,2};
    warning('segment_ABD:legacyInput', ...
        ['sig has no .thc/.abd columns; using columns 1 and 2. Pass the timetable ' ...
         'returned by segment_RIP to guarantee identical filtering.']);
else
    error('segment_ABD: sig must contain .thc and .abd (output of segment_RIP).');
end
if ~ismember('insp', fr1.Properties.VariableNames) || ~ismember('exp', fr1.Properties.VariableNames)
    error('segment_ABD: fr1 must contain .insp and .exp columns.');
end
t_vec = sig.Time;

fprintf('\n===== segment_ABD (v2) =====\n');
if isprop(sig.Properties.CustomProperties,'filter_info')
    fi = sig.Properties.CustomProperties.filter_info;
    fprintf('Using the segment_RIP passband: %.2f-%.2f Hz, order %d (identical on both channels)\n', ...
        fi.fc_lo, fi.fc_hi, fi.order);
end

%% Build ordered, alternating thoracic event list
t_insp = fr1.Time(fr1.insp == 1);
t_exp  = fr1.Time(fr1.exp  == 1);
ev_times = [t_insp; t_exp];
ev_types = [ones(numel(t_insp),1); 2*ones(numel(t_exp),1)];
[ev_times, si] = sort(ev_times);
ev_types = ev_types(si);
n_ev = numel(ev_times);
if n_ev < 4
    error('segment_ABD: fr1 contains too few events for guided detection.');
end

%% Guided half-cycle detection
global_iqr = iqr(abd);
zci = @(v) find(diff(sign(v(:))));

abd_insp_t = t_vec([]); abd_exp_t = t_vec([]);
abd_insp_q = logical.empty(0,1); abd_exp_q = logical.empty(0,1);
n_win_insp = 0; n_miss_insp = 0; n_multi_insp = 0;
n_win_exp  = 0; n_miss_exp  = 0; n_multi_exp  = 0;

for k = 1:n_ev-1
    t_start = ev_times(k); t_end = ev_times(k+1);
    ev_a = ev_types(k); ev_b = ev_types(k+1);
    half_tb_k = seconds(t_end - t_start);
    t_end_ext = min(t_end + seconds(half_tb_k*phase_margin), t_vec(end));

    m = t_vec >= t_start & t_vec <= t_end_ext;
    sig_win = abd(m); t_win = t_vec(m);
    if numel(sig_win) < 3, continue; end

    is_quality = (max(sig_win) - min(sig_win)) >= amp_thr * global_iqr;
    d_win = gradient(sig_win);
    zc_win = zci(d_win); zc_win = zc_win(zc_win+1 <= numel(d_win));

    if ev_a == 2 && ev_b == 1
        n_win_insp = n_win_insp + 1;
        peaks = zc_win(d_win(zc_win+1) - d_win(zc_win) < 0);
        if isempty(peaks)
            [~,best] = max(sig_win); n_miss_insp = n_miss_insp + 1;
        elseif isscalar(peaks)
            best = peaks;
        else
            [~,bi] = max(sig_win(peaks)); best = peaks(bi); n_multi_insp = n_multi_insp + 1;
        end
        abd_insp_t(end+1,1) = t_win(best);                        %#ok<AGROW>
        abd_insp_q(end+1,1) = is_quality && t_win(best) <= t_end; %#ok<AGROW>
    elseif ev_a == 1 && ev_b == 2
        n_win_exp = n_win_exp + 1;
        troughs = zc_win(d_win(zc_win+1) - d_win(zc_win) > 0);
        if isempty(troughs)
            [~,best] = min(sig_win); n_miss_exp = n_miss_exp + 1;
        elseif isscalar(troughs)
            best = troughs;
        else
            [~,bi] = min(sig_win(troughs)); best = troughs(bi); n_multi_exp = n_multi_exp + 1;
        end
        abd_exp_t(end+1,1) = t_win(best);                         %#ok<AGROW>
        abd_exp_q(end+1,1) = is_quality && t_win(best) <= t_end;  %#ok<AGROW>
    end
end

%% Assemble
n_insp = numel(abd_insp_t); n_exp = numel(abd_exp_t);
insp_tbl = timetable(abd_insp_t, ones(n_insp,1), zeros(n_insp,1), double(abd_insp_q), ...
    'VariableNames', {'insp','exp','quality'});
exp_tbl  = timetable(abd_exp_t, zeros(n_exp,1), ones(n_exp,1), double(abd_exp_q), ...
    'VariableNames', {'insp','exp','quality'});
fr2 = sortrows([insp_tbl; exp_tbl]);
fr2.Properties.DimensionNames{1} = 'Time';

pct_quality = 100 * sum(fr2.quality) / max(height(fr2),1);
fprintf('Insp windows: %d  |  missed (inserted max): %d  |  disambiguated: %d\n', n_win_insp, n_miss_insp, n_multi_insp);
fprintf('Exp  windows: %d  |  missed (inserted min): %d  |  disambiguated: %d\n', n_win_exp,  n_miss_exp,  n_multi_exp);
fprintf('High-quality windows: %.0f%%\n', pct_quality);
fprintf('============================\n\n');

info = struct('n_win_insp',n_win_insp,'n_miss_insp',n_miss_insp,'n_multi_insp',n_multi_insp, ...
              'n_win_exp',n_win_exp,'n_miss_exp',n_miss_exp,'n_multi_exp',n_multi_exp, ...
              'pct_quality',pct_quality);

%% Plot
if do_plot
    figure('Name','segment_ABD','NumberTitle','off','Color','w');
    ax1 = subplot(2,1,1);
    plot(t_vec, thc, '-', 'Color',[0.20 0.50 0.85], 'DisplayName','thoracic filtered'); hold on
    plot(fr1.Time(fr1.insp==1), thc(ismember(t_vec, fr1.Time(fr1.insp==1))), 'b^', 'MarkerSize',8, 'DisplayName','thc insp');
    plot(fr1.Time(fr1.exp==1),  thc(ismember(t_vec, fr1.Time(fr1.exp==1))),  'rv', 'MarkerSize',8, 'DisplayName','thc exp');
    legend('Location','best','Box','off'); title('Thoracic reference events');
    ylabel('Thoracic (au)'); grid on; box off;

    ax2 = subplot(2,1,2);
    plot(t_vec, abd, '-', 'Color',[0.85 0.30 0.10], 'DisplayName','abdominal filtered'); hold on
    hi_q = fr2.quality == 1;
    sel = @(tt) abd(ismember(t_vec, tt));
    if any(hi_q)
        plot(fr2.Time(fr2.insp==1 &  hi_q), sel(fr2.Time(fr2.insp==1 &  hi_q)), 'b^','MarkerSize',8,'MarkerFaceColor','b','DisplayName','insp (HQ)');
        plot(fr2.Time(fr2.exp==1  &  hi_q), sel(fr2.Time(fr2.exp==1  &  hi_q)), 'rv','MarkerSize',8,'MarkerFaceColor','r','DisplayName','exp (HQ)');
    end
    if any(~hi_q)
        plot(fr2.Time(fr2.insp==1 & ~hi_q), sel(fr2.Time(fr2.insp==1 & ~hi_q)), 'b^','MarkerSize',8,'LineStyle','none','DisplayName','insp (LQ)');
        plot(fr2.Time(fr2.exp==1  & ~hi_q), sel(fr2.Time(fr2.exp==1  & ~hi_q)), 'rv','MarkerSize',8,'LineStyle','none','DisplayName','exp (LQ)');
    end
    hold off; legend('Location','best','Box','off');
    ylabel('Abdominal (au)'); xlabel('Time'); grid on; box off;
    title('Abdominal detected events (filled = high quality)');
    linkaxes([ax1 ax2],'x');
end
end
