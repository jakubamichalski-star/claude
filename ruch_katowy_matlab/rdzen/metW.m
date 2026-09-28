function W = metW(met, hr)
%METW Wiatr z warstw METCM [W1; W2; W3] (układ ogniowy) na wysokości hr nad MDP.
W = [pchipVal(met.f1, hr); 0; pchipVal(met.f3, hr)];
end
