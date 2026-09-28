function o = diagnose(P, wind, t, y)
%DIAGNOSE Wielkości pochodne stanu (zapis przebiegów).
%   o – wiersz 1×41 w kolejności simChannels(); pola 'ev' i 'reinit' uzupełnia simulate.
if P.frozen, hA = P.hM; else, hA = P.hM + y(2); end
[Tair, ~, rho, aS] = atmAt(P, hA);
T = quatToT(y(10:13));
[psi, theta] = eulerFromT(T);
X2 = y(2); Ux = y(4); Uy = y(5); Uz = y(6);
[W, f] = wind.W(t, X2);
WI = [W(1), W(3), -W(2)];
WZ = [T(1) * WI(1) + T(2) * WI(2) + T(3) * WI(3), T(4) * WI(1) + T(5) * WI(2) + T(6) * WI(3), T(7) * WI(1) + T(8) * WI(2) + T(9) * WI(3)];
Vx = Ux - WZ(1); Vy = Uy - WZ(2); Vz = Uz - WZ(3);
V = sqrt(Vx * Vx + Vy * Vy + Vz * Vz);
UI = [T(1) * Ux + T(4) * Uy + T(7) * Uz, T(2) * Ux + T(5) * Uy + T(8) * Uz, T(3) * Ux + T(6) * Uy + T(9) * Uz];
U1 = UI(1); U2 = -UI(3); U3 = UI(2);
gr = wind.baseGrad(X2);
% ξ̇ z prawej strony równań: dV̄/dt (składowe zero-roll) = dŪ/dt − (T·dW̄_I/dt − Ω̄ × W̄)
[dy, clamp] = rhs(P, wind, t, y);
fr = wind.dfdt(t, X2, U2); dWg = wind.dW;
dWF = gr * U2 + fr * dWg;
dWI = [dWF(1), dWF(3), -dWF(2)];
tth = tan(theta); q = y(8); r = y(9); Ox = -r * tth;
dWZ = [T(1) * dWI(1) + T(2) * dWI(2) + T(3) * dWI(3) - (q * WZ(3) - r * WZ(2)), ...
       T(4) * dWI(1) + T(5) * dWI(2) + T(6) * dWI(3) - (r * WZ(1) - Ox * WZ(3)), ...
       T(7) * dWI(1) + T(8) * dWI(2) + T(9) * dWI(3) - (Ox * WZ(2) - q * WZ(1))];
dVx = dy(4) - dWZ(1); dVy = dy(5) - dWZ(2); dVz = dy(6) - dWZ(3);
Vd = (Vx * dVx + Vy * dVy + Vz * dVz) / V;
alpha = atan2(Vz, Vx);                              % (B.2)
alphaY = atan2(Vy, Vx);                             % (B.3)
beta = asin(max(-1, min(1, Vy / V)));               % (B.4)
alphaTot = atan2(sqrt(Vy * Vy + Vz * Vz), Vx);      % (B.5) = acos(V_x/V)
phiA = atan2(-Vy, Vz);                              % (B.6), w układzie zero-roll
dxiRe = dVz / V - Vz * Vd / (V * V);                % ξ ≈ (V_z − i·V_y)/V
dxiIm = -dVy / V + Vy * Vd / (V * V);
o = [t, y(1), X2, y(3), U1, U2, U3, sqrt(Ux * Ux + Uy * Uy + Uz * Uz), V, V / aS, rho, aS, y(7), q, r, ...
     psi, theta, y(14), y(15), Vx, Vy, Vz, alpha, alphaY, beta, alphaTot, phiA, W(1), W(2), W(3), ...
     gr(1) * U2, gr(2) * U2, gr(3) * U2, f, wind.st.piece, dxiRe, dxiIm, 0, 0, clamp, Tair];
end
