function k = statusClass(s, linear)
%STATUSCLASS Klasa stanu z kodu stanu (R.post.status).
%   k = statusClass(s)         – przy bieżącej amplitudzie (funkcja opisująca):
%                                'gyro' | 'dynF' | 'dynS' | 'dynFS' | 'limit' | 'weak' | 'ok'
%   k = statusClass(s, true)   – teoria liniowa (α_tot → 0): 'gyro' | 'dynF' | 'dynS' | 'dynFS' | 'ok'
ST = statusBits();
if nargin < 2, linear = false; end
if bitand(s, ST.GYRO), k = 'gyro'; return; end
if linear, bF = ST.DYN_F0; bS = ST.DYN_S0; else, bF = ST.DYN_F; bS = ST.DYN_S; end
f = bitand(s, bF) > 0; g = bitand(s, bS) > 0;
if f && g, k = 'dynFS'; return; end
if f, k = 'dynF'; return; end
if g, k = 'dynS'; return; end
if ~linear && bitand(s, ST.LIMIT), k = 'limit'; return; end
if ~linear && bitand(s, ST.WEAK), k = 'weak'; return; end
k = 'ok';
end
