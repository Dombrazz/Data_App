function out = co2GasLineFitFunction(w,center,w_g,w_l,a1,a2,a3,c0,c1)

persistent gasLines gasLinesFreq
if isempty(gasLines)
    m = load('CO2_gas_lines.mat');
    gasLinesFreq = m.freq;
    gasLines = m.gasLinesCorrected;
end

anh = 12;

gl = interp1(gasLinesFreq,gasLines,w(:));

dw = w(2)-w(1);
if dw < 0 
    error('only works for increasing frequencies. Please flip the input data!')
end
%out = @(w,center,w_g,w_l,a1,a2,a3) ...
out = a1.*voigt(w,center,w_g,w_l) + ...
    a1.*a2.*voigt(w,center - anh,w_g,w_l) + ...
    a3.*gl + ...
    c0 + ...
    c1.*(w-center);
