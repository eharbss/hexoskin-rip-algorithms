# hexoskin-rip-algorithms

MATLAB functions for breath segmentation and breath-by-breath breathing-pattern metrics from dual-band respiratory inductance plethysmography (RIP), developed for the Hexoskin smart shirt but applicable to any two-belt (rib cage + abdomen) RIP system.

The two headline outputs are:

- **Thoraco-abdominal coordination (TAC)** — the signed rib cage–abdomen phase angle per breath, from the Hilbert analytic signal, aggregated with circular statistics so that paradoxical (antiphase) breathing is represented correctly.
- **Rib cage contribution (%RC)** — the rib cage share of breath depth, given an explicit calibration constant K.

## Functions

| file | role |
|---|---|
| `rip_bandpass.m` | The single respiratory bandpass definition (4th-order Butterworth, SOS form, default 0.1–1.5 Hz). Every other function calls it, so the two belts are always filtered identically. |
| `segmentation/segment_RIP.m` | Flow-reversal detection on the SNR-weighted fusion of both belts. Returns the breath events and the separately filtered thoracic and abdominal signals. |
| `segmentation/segment_ABD.m` | Abdominal events found inside half-cycle windows defined by the thoracic events. |
| `enrichment/enrich_BP.m` | Breath-by-breath timing, depth, %RC, TAC and per-breath phase concentration, with optional validity gating. |
| `enrichment/rip_calibrate_K.m` | Derives K by isovolume manoeuvre, least squares against a spirometer, or QDC, with bootstrap confidence intervals and a quality score. |
| `rip_quality.m` | Two-step validity gate: rib cage–abdomen coherence at the breathing frequency (with an optional accelerometer cadence cross-check), then per-breath phase concentration. |

## Requirements

MATLAB with the Signal Processing Toolbox and the Statistics and Machine Learning Toolbox. Developed and tested on R2026a. Add the repository root and both subfolders to the path:

```matlab
addpath(genpath('hexoskin-rip-algorithms'))
```

## Typical use

```matlab
% rip: regular timetable with columns [thc abd] at the device sample rate
[fr1, sig]  = segment_RIP(rip);
fr2         = segment_ABD(sig, fr1);

% K from a 20 s occluded-airway isovolume manoeuvre (preferred)
cal         = rip_calibrate_K('isovolume', rip(isoTR, :));

% validity gate; pass the accelerometer for running data
[q, qsum]   = rip_quality(sig(:,'thc'), sig(:,'abd'), fs, 'fr1', fr1, 'acc', acc);

bp          = enrich_BP(sig, fr1, 'fr2', fr2, 'K', cal.K, 'quality', q);
```

`bp.tac_ph` and `bp.rcp` hold the values as computed. The gate never overwrites them; gated copies are in `bp.tac_ph_g` and `bp.rcp_g`, and the reason for each exclusion is in `bp.gate` (0 ok, 1 failed coherence, 2 failed phase concentration).

## Notes on interpretation

**K drives absolute %RC, and it is hard to get right.** %RC moves about 7 percentage points per unit change in K near K ≈ 3, so K needs to be known to roughly ±18% for ±5-point accuracy. The *within-subject change* in %RC is much less sensitive to K than its absolute level, so report Δ%RC against a reference condition where possible. There is deliberately no automatic fallback for K: the common shortcut K = SD(RC)/SD(AB) over raw samples forces mean %RC to about 50% by construction.

**Range vs signed %RC.** The default uses each belt's peak-to-trough excursion, which is robust but cannot represent paradox. `'rcp_mode','signed'` uses the signed inspiratory excursion instead; when the abdomen moves inward on inspiration it returns values above 100%, which is the physiologically correct reading.

**Hexoskin polarity.** The Hexoskin RIP signal falls as lung volume rises, so a signal maximum marks end-expiration. `.insp` events are placed accordingly.

## Citation

If you use the segmentation, please cite:

> Harbour, E., Lasshofer, M., Genitrini, M., & Schwameder, H. (2021). Enhanced breathing pattern detection during running using wearable sensors. *Sensors*, 21(16), 5606.

## License

MIT — see [LICENSE](LICENSE).
