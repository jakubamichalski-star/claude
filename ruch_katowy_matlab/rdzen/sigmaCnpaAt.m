function s = sigmaCnpaAt(A, M, aDeg)
%SIGMACNPAAT Odchylenie standardowe Cnpa dla ℳ i α [°]: interpolacja między 1°, 3° i 5° (wektorowo w α).
s1 = sigmaAt(A, 'Cnpa@1', M); s3 = sigmaAt(A, 'Cnpa@3', M); s5 = sigmaAt(A, 'Cnpa@5', M);
s = s3 + (s5 - s3) * (aDeg - 3) / 2;
lo = aDeg <= 3; s(lo) = s1 + (s3 - s1) * (aDeg(lo) - 1) / 2;
s(aDeg <= 1) = s1; s(aDeg >= 5) = s5;
end
