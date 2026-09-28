function T = TfromEuler(psi, theta, phi)
%TFROMEULER Macierz (A.8) wprost z kątów Eulera – do testu zgodności z kwaternionem (wektor 1×9, wierszami).
cps = cos(psi); sps = sin(psi); cth = cos(theta); sth = sin(theta); cph = cos(phi); sph = sin(phi);
T = [cth * cps, cth * sps, -sth, ...
     sph * sth * cps - cph * sps, cph * cps + sph * sth * sps, sph * cth, ...
     cph * sth * cps + sph * sps, cph * sth * sps - sph * cps, cph * cth];
end
