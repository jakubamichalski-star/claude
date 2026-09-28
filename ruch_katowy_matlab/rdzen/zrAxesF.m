function ax = zrAxesF(psi, theta)
%ZRAXESF Wersory x̄, ȳ, z̄ układu zero-roll w składowych układu ogniowego (φ = 0).
%   ax.x, ax.y, ax.z – wektory 1×3.
cps = cos(psi); sps = sin(psi); cth = cos(theta); sth = sin(theta);
ax.x = [cth * cps, sth, cth * sps];
ax.y = [-sps, 0, cps];
ax.z = [sth * cps, -cth, sth * sps];
end
