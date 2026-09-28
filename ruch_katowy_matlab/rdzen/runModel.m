function R = runModel(A, cfg, onProgress)
%RUNMODEL Pełne obliczenie: 6DOF + teoria liniowa + podsumowanie.
%   R = runModel(A, cfg, onProgress)
%   R.n, R.ch (przebiegi 6DOF), R.post (teoria liniowa, postprocess), R.meta, R.summary,
%   R.cfg (konfiguracja użyta w obliczeniu), R.derived (m, Ix, IT, d, S, kp, dXcg, p0, QE), R.P.
if nargin < 2, cfg = struct(); end
if nargin < 3, onProgress = []; end
R = simulate(A, cfg, onProgress);
t0 = tic;
R.post = postprocess(R);
R.summary = summarize(R);
R.summary.post_ms = 1000 * toc(t0);
R.cfg = R.P.cfg;
R.derived = struct('m', R.P.m, 'Ix', R.P.Ix, 'IT', R.P.IT, 'd', R.P.d, 'S', R.P.S, 'kp', R.P.kp, ...
  'dXcg', R.P.dXcg, 'p0', R.P.p0, 'QE', R.P.QE);
end
