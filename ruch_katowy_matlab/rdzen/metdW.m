function dW = metdW(met, hr)
%METDW Gradient wiatru z warstw METCM dW/dh [1/s] na wysokości hr nad MDP.
dW = [pchipDer(met.f1, hr); 0; pchipDer(met.f3, hr)];
end
