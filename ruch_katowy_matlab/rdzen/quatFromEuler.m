function l = quatFromEuler(psi, theta, phi)
%QUATFROMEULER Kwaternion λ̄ = [λ0 λ1 λ2 λ3] z kątów Eulera ψ, θ, φ (sekwencja 3-2-1).
cps = cos(psi / 2); sps = sin(psi / 2); cth = cos(theta / 2); sth = sin(theta / 2);
cph = cos(phi / 2); sph = sin(phi / 2);
l = [cph * cth * cps + sph * sth * sps, sph * cth * cps - cph * sth * sps, ...
     cph * sth * cps + sph * cth * sps, cph * cth * sps - sph * sth * cps];
end
