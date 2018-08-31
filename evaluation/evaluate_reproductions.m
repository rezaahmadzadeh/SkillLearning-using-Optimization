function performanceMeasures = evaluate_reproductions(demos,repros)

[performanceMeasures.SEA.mean, performanceMeasures.SEA.std, performanceMeasures.SEA.list] = compute_SEA_stats(demos,repros);
[performanceMeasures.SSE.mean, performanceMeasures.SSE.std, performanceMeasures.SSE.list] = compute_SSE_stats(demos,repros);
[performanceMeasures.DTWD.mean, performanceMeasures.DTWD.std, performanceMeasures.DTWD.list] = compute_DTWD_stats(demos,repros);