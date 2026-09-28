function [psi, theta, phi] = eulerFromT(T)
%EULERFROMT Kąty Eulera ψ, θ, φ z macierzy T (wektor 1×9, wierszami).
psi = atan2(T(2), T(1));
theta = asin(max(-1, min(1, -T(3))));
phi = atan2(T(6), T(9));
end
