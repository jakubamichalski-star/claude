function key = sigmaKey(name)
%SIGMAKEY Nazwa współczynnika z tabeli błędów jako nazwa pola struktury (np. 'Cnpa@1' → 'Cnpa_1').
key = regexprep(name, '[^A-Za-z0-9_]', '_');
end
