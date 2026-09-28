function cfg = defaultConfig()
%DEFAULTCONFIG Konfiguracja domyślna obliczeń.
%   norm   – normowanie CYpa, Cnpa, Clp, Cmq, CNq: 'pd/2V' (domyślne) albo 'pd/V'
%   mass   – m [kg], Ix, IT [kg·m²], XCG_mm – położenie środka masy od wierzchołka [mm]
%   mult   – mnożniki Ix, IT, p0 do badania wrażliwości
%   f, fd  – f_C* i f_ΔC* (B.9); ΔC = σ z tabel błędów pliku
%   fire   – u0 [m/s], QE (jednostka QEunit: 'deg' | 'mil6400' | 'mil6000'), tc – skok gwintu [kal.],
%            bore_mm – kaliber lufy, hM – wysokość stanowiska [m], q0, r0 [rad/s],
%            alphaTot0_deg, phiA0_deg – początkowe odchylenie osi od lufy
%   earth  – g, coriolis (true/false), lat_deg, az_deg – azymut strzału
%   wind   – base.W [W1 W2 W3] albo base.table [h W1 W2 W3] (N×4, N ≥ 2), gust – podmuch:
%            enabled, mode ('time' | 'height'), h1, branch ('asc' | 'desc' | 'both' | 'apex'),
%            A, azimuth_deg (dokąd wieje, od osi 1 ku osi 3), W2, shape ('step' | 'trapez' | 'cos'),
%            tr [s], td [s], dh [m], HL [m]
%   met    – warstwy METCM: enabled, hMDP ([] → h_M), dirUnit ('mil' | 'deg'),
%            zones (1×32: dir – skąd wieje, spd [m/s], dT [K], dRho [%])
%   num    – stepsPerPeriod, samplesPerPeriod, tMax [s], alphaWarn_deg, alphaStop_deg,
%            thetaStop_deg, alphaEWarn_deg, sigmaEps [1/s], wkb
%   test   – frozen (współczynniki zamrożone: g = 0, C_D = 0, C_spin = 0, ρ = const), gravity
Z = metcmZones();
cfg.norm = 'pd/2V';
cfg.mass = struct('m', 40.3, 'Ix', 0.161, 'IT', 1.89439, 'XCG_mm', 593.02);
cfg.mult = struct('Ix', 1, 'IT', 1, 'p0', 1);
cfg.f = struct('CD0', 1, 'CD2', 1, 'CLa', 1, 'CMa', 1, 'CMq', 1, 'CNq', 1, 'Cmagf', 1, 'Cmagm', 1, 'Cspin', 1);
cfg.fd = struct('CD0', 0, 'CLa', 0, 'CMa', 0, 'CMq', 0, 'Cmagm', 0, 'Cspin', 0);
cfg.fullCmqTransfer = false;
cfg.fire = struct('u0', 800, 'QE', 30, 'QEunit', 'deg', 'tc', 20, 'bore_mm', 155, 'hM', 0, 'q0', 1.0, 'r0', 0, ...
  'alphaTot0_deg', 0, 'phiA0_deg', 0);
cfg.earth = struct('g', 9.80665, 'coriolis', false, 'lat_deg', 52.2, 'az_deg', 0);
cfg.wind.base = struct('W', [0 0 0], 'table', []);
cfg.wind.gust = struct('enabled', false, 'mode', 'time', 'h1', 1000, 'branch', 'asc', 'A', 10, 'azimuth_deg', 90, ...
  'W2', 0, 'shape', 'step', 'tr', 0.05, 'td', 0, 'dh', 50, 'HL', 0);
cfg.met = struct('enabled', false, 'hMDP', [], 'dirUnit', 'mil');
cfg.met.zones = repmat(struct('dir', 0, 'spd', 0, 'dT', 0, 'dRho', 0), 1, numel(Z.n));
cfg.num = struct('stepsPerPeriod', 100, 'samplesPerPeriod', 20, 'tMax', 400, 'alphaWarn_deg', 10, 'alphaStop_deg', 15, ...
  'thetaStop_deg', 85, 'alphaEWarn_deg', 5, 'sigmaEps', 0.02, 'wkb', true);
cfg.test = struct('frozen', false, 'gravity', false);
end
