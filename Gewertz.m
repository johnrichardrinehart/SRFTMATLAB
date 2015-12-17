function [num,den] = Gewertz(R0, resistance_denominator_coefficients)
%%% Gewertz(res_num_coeff,res_den_coeff) takes two arguments:
%%% 1) Coefficients of the numerator of the resistance polynomial R. These
%%% must be specified in array where the coefficents are ordered greatest
%%% to least.
%%% 2) Coefficients of the denominator of the resistance polynomial R.
%%% These must be specified in the same way as the numerator.
%%% e.g. R(w^2) = R(-s^2) = (2.2*w^2+1)/(1+2*w^2+3*w^4-5*w^6)
%%% = (-2.2*s^2+1)/(1-2*s^2+3*s^4+5*s^6).
%%% res_num_coeff = [-2.2 0 1]; res_denom_coeff = [5, 0, 3, 0, -2, 0, 1]
% The resistance_numerator_cofficients are the coefficients of the
% numerator polynomial of r(-s^2), the resistance_denominator coefficients
% are the coefficients of the denominator polynomial of r(-s^2). These
% coefficients must adhere to the MATLAB standard of ordering coefficients
% from greatest power to least power
% If the number of numerator coeffecients is n, the number of denominator
% coefficients should be 2n-1 (i.e., numerator is 3rd order
%(4 coefficients), denominator is 6th order (7 coefficients)).

resistance_denominator_coefficients = wCoeffsTosCoeffs(resistance_denominator_coefficients);
ResistanceSpectralRoots = roots(resistance_denominator_coefficients); % Find all the roots of the denominator
ResistanceSpectralRoots = RemoveSmallRealPart(ResistanceSpectralRoots);
NegSpectralRoots = BuildNegPolyArray(ResistanceSpectralRoots); % Remove the zero entries corresponding to RHP roots

if isrow(poly(NegSpectralRoots))
    NegSpectralPoly = poly(NegSpectralRoots)'; % PosSpectralDen is the vector of polynomial coefficients
end
NegSpectralPoly = sqrt(resistance_denominator_coefficients(1))*NegSpectralPoly;
ds = NegSpectralPoly;

if length(NegSpectralRoots) ~= length(ResistanceSpectralRoots)/2
    disp('Bad spectral factorization');
end

% Form the matrix of d coefficients
D = zeros(length(NegSpectralPoly),length(NegSpectralPoly));
for i = 1:length(NegSpectralPoly) % row
    for j = 1:length(NegSpectralPoly) % column
        Dval = (2*i-1)-(j-1);
        if Dval > 0 && Dval <= length(ds)
            D(i,j) = ds(Dval);
        end
    end
end

SignFlipper = eye(size(D));
for i = 1 : length(D)
    SignFlipper(i,i) = -1*(1i)^(i+i);
end

D = D*SignFlipper;

ResistanceNumeratorPolynomial = zeros(1,length(D));
if isrow(ResistanceNumeratorPolynomial)
    ResistanceNumeratorPolynomial = ResistanceNumeratorPolynomial';
end
ResistanceNumeratorPolynomial(end) = R0;
num = D\ResistanceNumeratorPolynomial; %Conforms to MATLAB standard with highest coeff on top
den = NegSpectralPoly; % Conform to MATLAB standard, put highest coeff on top

%%%% Test code
num_star = zeros(length(num),1);
den_star = zeros(length(den),1);
for i = 1:length(num)
    num_star(end-(i-1)) = (-1)^(i-1)*num(end-(i-1));
end
for i = 1:length(den)
    den_star(end-(i-1)) = (-1)^(i-1)*den(end-(i-1));
end
reconstr_num = .5*(conv(num,den_star)+conv(num_star,den));
reconstr_den = conv(den,den_star);

end

function NewArray = RemoveSmallRealPart(OrigArray)
NewArray = size(OrigArray);
for i = 1:length(OrigArray)
    value = abs(real(OrigArray(i))/abs(OrigArray(i)));
    if value < 1e-7
        NewArray(i) = 1i*imag(OrigArray(i));
    else
        NewArray(i) = OrigArray(i);
    end
end
end

function NegPolyArray = BuildNegPolyArray(ResistanceRoots)
NegPolyArray = zeros(length(ResistanceRoots),1);
for i = 1:length(NegPolyArray)
    if real(ResistanceRoots(i)) < 0
        NegPolyArray(i) = ResistanceRoots(i);
    elseif real(ResistanceRoots(i)) == 0 && imag(ResistanceRoots(i)) < 0
        NegPolyArray(i) = ResistanceRoots(i);
    end
end
NegPolyArray(NegPolyArray == 0) = [];
end

function sCoeffs = wCoeffsTosCoeffs(Coeffs)
% This function converts coefficients of a function of w^2 to coefficients
% of a function of (-s)^2 so that spectral decomposition makes sense. It
% assumes that Coeffs(end) is the 0th order term.

sCoeffs = zeros(size(Coeffs));
for i = 1:length(Coeffs)
    power = i-1;
    powerIdx = length(Coeffs)-(i-1);
    if mod(power-2,4) == 0 % is it the 3rd term (second order) or 7th term (sixth order)
        sCoeffs(powerIdx) = -1*Coeffs(powerIdx);
    else
        sCoeffs(powerIdx) = Coeffs(powerIdx);
    end
end
end