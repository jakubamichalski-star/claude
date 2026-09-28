function u = uniformXi(R, t0, t1)
%UNIFORMXI Próbki ξ = α_tot·e^{iφ_A} z przedziału [t0, t1] bez próbek zdarzeń (równoodległe).
%   u.ts, u.z (zespolone), u.dt, u.uniform (czy odstępy są równe).
ch = R.ch;
sel = ch.t >= t0 & ch.t <= t1 & ~ch.ev;
u.idx = find(sel);
u.ts = ch.t(sel); u.z = ch.alphaTot(sel) .* exp(1i * ch.phiA(sel));
u.dt = (u.ts(end) - u.ts(1)) / (numel(u.ts) - 1);
dev = max(abs(u.ts - u.ts(1) - (0:numel(u.ts) - 1).' * u.dt));
u.uniform = dev < 1e-6 * u.dt;
end
