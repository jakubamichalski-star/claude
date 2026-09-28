function r = quadRoots(B, C)
%QUADROOTS Pierwiastki s² + B·s + C = 0 (B, C zespolone), r – wektor 2×1.
%   Pierwiastek z wyróżnika: Re ≥ 0, znak Im jak znak Im wyróżnika (dla Im = 0: dodatni).
disc = B * B - 4 * C;
ad = abs(disc); dr = real(disc);
sq = complex(sqrt(max(0, (ad + dr) / 2)), sqrt(max(0, (ad - dr) / 2)));
if imag(disc) < 0, sq = conj(sq); end
r = [(-B + sq) / 2; (-B - sq) / 2];
end
