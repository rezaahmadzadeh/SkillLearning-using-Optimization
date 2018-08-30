function performanceMeasures = evaluate_reproductions(demos,repros)

[performanceMeasures.SEA.mean, performanceMeasures.SEA.std] = compute_SEA_stats(demos,repros);
[performanceMeasures.SSE.mean, performanceMeasures.SSE.std] = compute_SSE_stats(demos,repros);
[performanceMeasures.DTWD.mean, performanceMeasures.DTWD.std] = compute_DWTD_stats(demos,repros);