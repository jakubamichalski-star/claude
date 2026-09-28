function ST = statusBits()
%STATUSBITS Bity kodu stanu (R.post.status).
%   GYRO – S_g < 1; DYN_F, DYN_S – σ > σ_eps przy bieżącej amplitudzie; WEAK – |ᾱ_e| > próg;
%   AWARN – α_tot > próg ostrzeżenia; DYN_F0, DYN_S0 – σ > 0 przy α_tot → 0;
%   MACH – ℳ poza zakresem tabel; LIMIT – cykl graniczny.
ST = struct('GYRO', 1, 'DYN_F', 2, 'DYN_S', 4, 'WEAK', 8, 'AWARN', 16, 'DYN_F0', 32, 'DYN_S0', 64, 'MACH', 128, 'LIMIT', 256);
end
