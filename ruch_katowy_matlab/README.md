# Ruch kątowy pocisku stabilizowanego obrotowo – wersja MATLAB

Port aplikacji `ruch_katowy_etap4_metcm.html` (model 6DOF i teoria liniowa w notacji AEP-96,
dane aerodynamiczne Spin2004, warstwy meteo METCM) do MATLAB-a. Rdzeń obliczeniowy przepisano
funkcja po funkcji. Nazwy funkcji i oznaczenia są takie jak w oryginale, więc łatwo porównać oba kody.

## Uruchomienie

Otwórz `ruch_katowy.m` i uruchom go w całości (F5) albo sekcjami (Ctrl+Enter). Skrypt:

1. wczytuje plik Spin2004 (`dane/155OFdMKM.txt`), weryfikuje go (V1–V6) i buduje model aerodynamiczny,
2. ustawia parametry (`cfg`, pola opisane w `help defaultConfig`),
3. opcjonalnie wczytuje komunikat METCM (`uzyjMeteo = true`),
4. wypisuje wielkości u wylotu (S_g, S_d, ν, σ, krok całkowania),
5. liczy tor (6DOF + teoria liniowa) i wypisuje wyniki,
6. rysuje wykresy,
7. eksportuje przebiegi do CSV (sekcja zakomentowana),
8. uruchamia scenariusze (a)–(f) i testy modelu (sekcje zakomentowane).

Minimalny przykład:

```matlab
addpath('rdzen', 'wykresy', 'narzedzia');
A = buildAero(parseSpin2004(fileread('dane/155OFdMKM.txt')));
cfg = defaultConfig();  cfg.fire.QE = 45;
R = runModel(A, cfg);
pokazWyniki(R);  plotResults(R);
```

## Odpowiedniki zakładek aplikacji HTML

| Zakładka HTML | MATLAB |
|---|---|
| Dane aerodynamiczne | `parseSpin2004`, `verifySpin`, `buildAero`, `pokazDane` |
| Parametry, „Wielkości u wylotu” | `defaultConfig`, `deepMerge`, `pokazWylot` |
| Meteo (METCM) | `parseMETCM`, `metcmToConfig`, `metProfile`, `plotMeteo` |
| Scenariusze | `scenariusze`, `uruchomScenariusz(A, id, wariant, cfg)` |
| Wyniki | `runModel` → `R.summary`, `pokazWyniki` (także tabela przebiegu w wybranych chwilach) |
| Wykresy | `plotResults`: `plotTimeSeries`, `plotRosette`, `plotStabDiagram`, `plotTrajectory` |
| Eksport CSV | `exportCSV(R, plik, 'matlab' \| 'excel')` – te same kolumny co w aplikacji |
| Testy | `runTests(A, cfg)` – T1–T14 |
| Teoria, O modelu | komentarze w nagłówkach funkcji rdzenia |

## Struktura

```
ruch_katowy.m        skrypt główny
dane/                155OFdMKM.txt (plik wbudowany w aplikację), przyklad_METCM.txt (komunikat syntetyczny)
rdzen/               rdzeń obliczeniowy (bez grafiki)
  atmosphere, atmAt                  atmosfera ISO 2533 (+ odchyłki METCM)
  parseSpin2004, verifySpin, buildAero, sigmaAt, sigmaCnpaAt
  aeroCoef, aeroCoefVec              współczynniki (ℳ, α_tot), κ_p, f_C, f_ΔC
  metcmZones, metProfile, metW, metdW, parseMETCM, metcmToConfig, pchipMake/Val/Der
  Wiatr                              wiatr bazowy, podmuch (klasa handle – stan podmuchu)
  defaultConfig, deepMerge, makeParams
  quatFromEuler, quatToT, TfromEuler, eulerFromT, zrAxesF
  initialState, rhs, projectZeroRoll, diagnose, simulate     6DOF, RK4, zdarzenia
  quadRoots, orderRoots, linPoint, stabBands, alphaEStanag   teoria liniowa, ᾱ_e (STANAG 4355)
  cmagmAtMach, describingMagnus, postprocess                 ξ_e, ramiona K_F, K_S, funkcja opisująca
  statusBits, statusClass, statusText, summarize, runModel, gustExcitation
  prony2, uniformXi, runTests
wykresy/             wykresy
narzedzia/           wydruki, eksport CSV, scenariusze
```

## Struktura wyniku `R = runModel(A, cfg)`

- `R.ch` – przebiegi 6DOF (wektory kolumnowe): `t, X1, X2, X3, U1, U2, U3, V, M, rho, p, q, r, psi, theta,
  alphaTot, phiA, W1, W2, W3, …` (pełna lista: `simChannels`). Kąty w radianach.
- `R.post` – teoria liniowa w każdej próbce: `Sg, Sd0, SdFA, SdSA, nuF0, nuS0, sigF0, sigS0, sFreA, …,
  xiRe, xiIm, xieRe, xieIm, aeAbs, KFm, KSm, KFabs, KSabs, status, …` (lista: `postChannels`).
- `R.meta` – krok, zdarzenia (`R.meta.events`), wierzchołek, punkt upadku.
- `R.summary` – podsumowanie i komunikaty (`R.summary.messages`), przedziały stanów.
- `R.cfg` – użyta konfiguracja, `R.derived` – m, Ix, IT, d, S, κ_p, ΔX_CG, p₀, QE.

## Zgodność z oryginałem

Kod sprawdzono w GNU Octave 8.4 na tych samych danych co aplikację HTML (Node.js):

- lot domyślny (u₀ = 800 m/s, QE = 30°, 57 s lotu, 125 396 kroków RK4, 25 082 próbki): donośność,
  derywacja, wierzchołek, ekstrema, wielkości u wylotu, amplitudy ramion, S_d, ξ_e, ᾱ_e i kody stanu
  są identyczne z oryginałem co do 10–12 cyfr znaczących;
- lot z podmuchem warstwowym, warstwami METCM i przyspieszeniem Coriolisa oraz lot z podmuchem
  trapezowym dają te same zdarzenia, przebiegi i komunikaty;
- testy wczytania V1–V6 i testy modelu T1–T14 dają te same wartości co w aplikacji.

## Czego nie przeniesiono

- Interfejsu HTML (formularze, kursor czasu, odtwarzanie) – zastępują go struktura `cfg` i funkcje wykresów.
- Widoku 3D z animacją pocisku.
- Nakładania przebiegu referencyjnego z CSV na wykresy. Można go dorysować zwykłym `plot` na osiach z `plotTimeSeries`.
- Zapisu i odczytu konfiguracji w JSON. W MATLAB-ie służą do tego `jsonencode(cfg)` i `jsondecode`
  (po `jsondecode` pola puste i tablice stref trzeba sprawdzić).

## Wymagania

MATLAB R2016b lub nowszy, bez dodatkowych toolboksów. Kod działa też w GNU Octave ≥ 8. W Octave
`pokazWyniki` nie zwraca obiektu `table`, a obliczenia są kilkadziesiąt razy wolniejsze.
