function r = orderRoots(r)
%ORDERROOTS Kolejność pierwiastków: najpierw F (większe |Im s|), przy równych – większa część rzeczywista.
a = abs(imag(r(1))); b = abs(imag(r(2))); eps_ = 1e-12 * (a + b + 1e-300);
if a > b + eps_ || (abs(a - b) <= eps_ && real(r(1)) >= real(r(2))), return; end
r = r([2 1]);
end
