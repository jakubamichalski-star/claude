%% Ruch kątowy pocisku stabilizowanego obrotowo
% Model 6DOF i teoria liniowa w notacji AEP-96, dane aerodynamiczne Spin2004.
% Wersja MATLAB aplikacji „ruch_katowy_etap4_metcm.html”.
%
% Rdzeń obliczeniowy całkuje równania ruchu bryły sztywnej osiowo-symetrycznej (6DOF, układ zero-roll)
% z siłami i momentami wg AEP-96 Ed. A v1, dodatek B. W każdej zapisanej chwili liczy wielkości
% liniowej teorii ruchu kątowego (sekcja 4 specyfikacji), kąt równowagi ᾱ_e wg STANAG 4355,
% rozwiązanie szczególne ξ_e oraz ramiona epicykliczne K_F, K_S do nałożenia na wynik 6DOF.
%
% Układ ogniowy (1, 2, 3): 1 w kierunku celu, 2 w górę, 3 w prawo. Wysokości nad poziomem wylotu;
% h_M służy tylko do atmosfery. ξ = α_tot·e^{iφ_A}: Re ξ – nos w górę względem V̄, Im ξ – nos w prawo.
%
% Skrypt jest podzielony na sekcje (Ctrl+Enter uruchamia bieżącą sekcję).

clear; close all; clc;
katalog = fileparts(mfilename('fullpath'));
addpath(fullfile(katalog, 'rdzen'), fullfile(katalog, 'wykresy'), fullfile(katalog, 'narzedzia'));

%% 1. Dane aerodynamiczne Spin2004
% Plik wyników Spin2004 (ASCII). Tabele są rozpoznawane po nagłówkach kolumn.
plikSpin = fullfile(katalog, 'dane', '155OFdMKM.txt');
[~, nazwa, roz] = fileparts(plikSpin);
Spin = parseSpin2004(fileread(plikSpin), [nazwa roz]);
Wer = verifySpin(Spin);                     % testy wczytania (tolerancje z zaokrągleń zapisu)
A = buildAero(Spin);                        % przypisanie kolumn do oznaczeń AEP-96
pokazDane(Spin, A, Wer);

%% 2. Parametry obliczeń
% Wartości domyślne – zob. help defaultConfig. Poniżej najczęściej zmieniane pola.
cfg = defaultConfig();
cfg.norm = 'pd/2V';                         % normowanie CYpa, Cnpa, Clp, Cmq, CNq: 'pd/2V' albo 'pd/V'
% masa i bezwładność
cfg.mass.m = 40.3;                          % [kg]
cfg.mass.Ix = 0.161;                        % [kg·m²]
cfg.mass.IT = 1.89439;                      % [kg·m²]
cfg.mass.XCG_mm = 593.02;                   % środek masy od wierzchołka [mm]
% strzał
cfg.fire.u0 = 800;                          % prędkość wylotowa [m/s]
cfg.fire.QE = 30;  cfg.fire.QEunit = 'deg'; % kąt podniesienia ('deg', 'mil6400', 'mil6000')
cfg.fire.tc = 20;                           % skok gwintu [kalibry]
cfg.fire.bore_mm = 155;                     % kaliber lufy [mm] (tylko do p₀)
cfg.fire.hM = 0;                            % wysokość stanowiska n.p.m. [m]
cfg.fire.q0 = 1.0;  cfg.fire.r0 = 0;        % zaburzenie wylotowe [rad/s]
% Ziemia
cfg.earth.coriolis = false;  cfg.earth.lat_deg = 52.2;  cfg.earth.az_deg = 0;
% wiatr bazowy: stały [W1 W2 W3] albo tabela [h W1 W2 W3] (N×4)
cfg.wind.base.W = [0 0 0];
cfg.wind.base.table = [];
% podmuch (opcjonalnie)
cfg.wind.gust.enabled = false;
cfg.wind.gust.mode = 'time';                % 'time' – od przecięcia h₁, 'height' – warstwa wysokości
cfg.wind.gust.h1 = 1000;                    % [m]
cfg.wind.gust.branch = 'asc';               % 'asc', 'desc', 'both', 'apex'
cfg.wind.gust.A = 10;  cfg.wind.gust.azimuth_deg = 90;  cfg.wind.gust.W2 = 0;
cfg.wind.gust.shape = 'step';               % 'step', 'trapez', 'cos'
cfg.wind.gust.tr = 0.05;  cfg.wind.gust.td = 0;  cfg.wind.gust.dh = 50;  cfg.wind.gust.HL = 0;
% numeryka
cfg.num.stepsPerPeriod = 100;               % kroków na okres najszybszej składowej u wylotu
cfg.num.samplesPerPeriod = 20;              % próbek zapisu na ten okres
cfg.num.tMax = 400;                         % [s]

%% 3. Meteo: warstwy komunikatu METCM (opcjonalnie)
% Strefy jak w komunikacie METCM (STANAG 4082): linia 00 – powierzchnia MDP, strefy 01–31 do 30 km nad MDP.
% Dla każdej strefy: wiatr (kierunek, z którego wieje, i prędkość) oraz odchyłki temperatury (wirtualnej)
% i gęstości od ISO 2533 na wysokości środka strefy. Między środkami stref – interpolacja PCHIP.
% Warstwy zastępują wiatr bazowy W̄₀; podmuch działa dalej.
uzyjMeteo = false;
if uzyjMeteo
  cfg.met.dirUnit = 'mil';                  % kierunek w tysięcznych (6400) albo 'deg'
  plikMet = fullfile(katalog, 'dane', 'przyklad_METCM.txt');
  M = parseMETCM(fileread(plikMet), cfg.fire.hM);
  for k = 1:numel(M.warnings), fprintf('METCM: %s\n', M.warnings{k}); end
  fprintf('METCM: %d stref, MDP = %g m n.p.m. (%s)\n', numel(M.zones), M.hMDP, M.hMDPsource);
  cfg.met = metcmToConfig(M, cfg.met);      % włącza warstwy
  % ręcznie: cfg.met.enabled = true; cfg.met.zones(k) = struct('dir', …, 'spd', …, 'dT', …, 'dRho', …);
  plotMeteo(cfg);
end

%% 4. Wielkości u wylotu (teoria liniowa, α_tot → 0) – bez całkowania toru
pokazWylot(A, cfg);

%% 5. Obliczenie toru: 6DOF + teoria liniowa + podsumowanie
postep = @(t, h) fprintf('  t = %6.2f s, h = %6.0f m\n', t, h);
R = runModel(A, cfg, postep);
pokazWyniki(R);
% tabela przebiegu w wybranych chwilach jako obiekt table:  T = pokazWyniki(R);

%% 6. Wykresy
% opcje: tRange – zakres czasu, cursor – chwila dla rozety/diagramu/toru, rosWin – okno rozety [s]
plotResults(R, struct('tRange', [0 5], 'cursor', 0.3, 'rosWin', 1));
% cały lot:  plotResults(R);

%% 7. Eksport przebiegów do CSV
% exportCSV(R, fullfile(katalog, 'przebieg.csv'));             % „,” i kropka dziesiętna (MATLAB)
% exportCSV(R, fullfile(katalog, 'przebieg_excel.csv'), 'excel'); % „;” i przecinek (Excel PL)
% zapis całego wyniku:  save(fullfile(katalog, 'wynik.mat'), 'R');

%% 8. Scenariusze (a)–(f)
% Każdy scenariusz ustawia parametry (wartości domyślne plus zmiany scenariusza; masa, bezwładność,
% środek masy i normowanie zostają z cfg), liczy tor, wypisuje wynik i rysuje wykresy.
%   'a'  – lot wzorcowy: zanikanie zaburzenia wylotowego
%   'b'  – zbyt wolny obrót: S_g < 1 (wariant 1: t_c = 25 kal.)
%   'c1' – niestateczność dynamiczna, składowa wolna (cykl graniczny)
%   'c2' – niestateczność dynamiczna, składowa szybka
%   'd'  – przejście przez zakres okołodźwiękowy
%   'e'  – duże QE: przestabilizowanie w pobliżu wierzchołka (wariant 1: podmuch w wierzchołku)
%   'f1' – podmuch boczny na zadanej wysokości, skokowy
%   'f2' – podmuch boczny na zadanej wysokości, łagodny (uruchom po 'f1', by dostać stosunek)
% Rs = uruchomScenariusz(A, 'a', [], cfg);
% Rs = uruchomScenariusz(A, 'e', 1, cfg);

%% 9. Testy modelu (sekcja 10 specyfikacji + testy dodatkowe; kilka minut)
% wynikiTestow = runTests(A, cfg);
