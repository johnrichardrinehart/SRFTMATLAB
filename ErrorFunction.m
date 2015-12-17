function ErrorValue = ErrorFunction(break_freqs, function_freqs, loadfunc, Ds, Rdc, gaingoal)
R = RsFromDs(break_freqs, function_freqs, Ds, Rdc);
X = XsFromDs(break_freqs, function_freqs, Ds);
Z = R + 1i*X;
z = loadfunc(function_freqs);
TG = TransducerGain(Z,z);
ErrorValue = gaingoal-TG;
ErrorValue = sum(ErrorValue.^2);