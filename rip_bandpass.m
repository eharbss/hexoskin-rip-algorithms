function [sos,g,info] = rip_bandpass(fs, fc_lo, fc_hi, ord)
%RIP_BANDPASS  The single definition of the respiratory bandpass.
%
% Every function in the dual-band RIP chain (segment_RIP, segment_ABD,
% enrich_BP, rip_quality) obtains its filter here so
% that the thoracic and abdominal channels are ALWAYS filtered identically.
% Unequal passbands between the two channels bias the thoraco-abdominal phase
% angle even though filtfilt is zero-phase, because the passband changes the
% waveform shape and therefore the Hilbert instantaneous phase.
%
% Syntax:
%   [sos,g]      = rip_bandpass(fs)
%   [sos,g,info] = rip_bandpass(fs, fc_lo, fc_hi, ord)
%
% Inputs:
%   fs     - sample rate (Hz)
%   fc_lo  - lower cutoff (Hz), default 0.1
%   fc_hi  - upper cutoff (Hz), default 1.5   (1.5 Hz = 90 br/min)
%   ord    - Butterworth order, default 4
%
% Outputs:
%   sos, g - second-order-section representation and gain. Use as
%              y = filtfilt(sos, g, x);    % offline, zero phase
%              y = g * sosfilt(sos, x);    % causal
%            SOS form is used rather than [b,a] because a 4th-order
%            0.1-1.5 Hz design at 128 Hz is poorly conditioned in
%            transfer-function form.
%   info   - struct with the realised design parameters
%
% Note on order: a 2nd-order design at [0.1 1.5] leaves substantial
% locomotor energy in the passband. On a reference running recording the
% abdominal breathing-band SNR was +1.6 dB with butter(2,...) and +24.6 dB
% with a 4th-order design.
%
% created by Eric Harbour (PLUS)

if nargin < 2 || isempty(fc_lo), fc_lo = 0.1; end
if nargin < 3 || isempty(fc_hi), fc_hi = 1.5; end
if nargin < 4 || isempty(ord),   ord   = 4;   end

fn = fs/2;
if fc_hi >= 0.95*fn
    fc_hi = 0.95*fn;
    warning('rip_bandpass:hiClamped','fc_hi clamped to %.3f Hz for fs = %.1f Hz.',fc_hi,fs);
end
if fc_lo <= 0 || fc_lo >= fc_hi
    error('rip_bandpass: require 0 < fc_lo < fc_hi.');
end

[z,p,k] = butter(ord, [fc_lo fc_hi]/fn, 'bandpass');
[sos,g] = zp2sos(z,p,k);

info = struct('fs',fs,'fc_lo',fc_lo,'fc_hi',fc_hi,'order',ord, ...
              'max_pole_radius',max(abs(p)));
end
