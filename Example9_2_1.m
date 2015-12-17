clear;
% number of breakpoints and fit order below
numbrkpts = 6;
%fitorder = 6; % make sure it's even or things will break
% load description below
zload = @(w)((1+1i*w*1.2).^-1+1i*w*2.3); % g = 1, Yc = 1i*w*1.2, Zl = i*w*2.3
% Break frequencies below
omega_max = 1.25;
break_omegas = linspace(0,omega_max,numbrkpts)';
% Frequencies over which to fit the polynomial
%optim_freqs = (0:1/20:omega_max+1)';
optim_freqs = break_omegas;
omega_one_idx = find(optim_freqs>1,1);
% Transducer gain below
gaingoal = ones(length(optim_freqs),1);
%gaingoal(omega_one_idx:end) = linspace(1,0,length(gaingoal)-omega_one_idx+1)';
gaingoal(omega_one_idx:end) = 0*ones(length(optim_freqs)-omega_one_idx+1,1);

R0 = 2.2;
% TODO: Make an initial guess function based on reactive component
% cancellation
Dinit = [-.3; -.4; -.3; -.6; -.6];

objective = @(D)ErrorFunction(break_omegas,optim_freqs,zload,D,R0,gaingoal);
sum_constraint = -1*ones(1,length(Dinit)); % Make sure the sum is less than R0
% x = fmincon(@myfun,x0,A (inequality),b,Aeq (equality),beq,lb,ub,@mycon)
[Ds, TG] = fmincon(objective,Dinit,sum_constraint,R0-.3,[],[],[],[],@CheckDecreasing); %TODO: Investigate other optimization scheme
%Ds = Dinit;
R = RsFromDs(break_omegas,optim_freqs,Ds,R0);

if (any(R < 0))
    warning(['Your fit didn''t determine no resistance at the highest '...
        'frequency. Check your fit.']);
    for i = 1:length(R)
        if R(i) < 0
            R(i) = 0.1;
        end
    end
end

X = XsFromDs(break_omegas, optim_freqs, Ds); % Obtain X to check TG

z = zload(optim_freqs);
Z = R + 1i*X;
fig=figure(1);
plot(optim_freqs,gaingoal,'g',optim_freqs,R,'k',optim_freqs,TransducerGain(Z,z),'b'); % Check TG
title(['Goal TG (green), Calculated TG (blue), '...
    'Resistance (black)']);

total_freqs = [-1*flipud(optim_freqs(2:end)); optim_freqs];
%normalized_R = R / R(1); % Turn the numerator into a "one".
total_R = [flipud(R(2:end)); R]; % Make the double-sided R, assuming even
%total_R = [fliplr(normalized_R(2:end)),normalized_R]; % Make the double-sided R, assuming even
InverseR = 1./total_R; % We assume the denominator of R is a polynomial

% Perform Polynomial Fitting
%[xData, yData] = prepareCurveData( total_omegas, InverseR );
%ft = fittype(sprintf('poly%d',fitorder));\
%[fitresult, gof] = fit( xData, yData, ft );

% Perform Custom Fitting
[xData, yData] = prepareCurveData( total_freqs, InverseR);
ft = fittype('a8*x.^8+a7*x.^7+a6*x.^6+a5*x.^5+a4*x.^4+a3*x.^3+a2*x.^2+a1*x.^1+a0','independent','x');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
[fitresult, gof] = fit(xData, yData, ft, opts);
% Below: First term is made highest order, consistent with roots and poly
resistance_den_coeffs = flipud(coeffvalues(fitresult)');

% Below: Check to make sure resistance_den_coeffs are good.
figure(2);plot(optim_freqs,R,'b',optim_freqs,1./polyval(resistance_den_coeffs,optim_freqs),'r');
title('Polynomial fit resistance (red) and Optimized Resistance (blue)');

% Make sure that the polynomial is non-negative
%figure(3);plot(-50:.01:50,polyval(resistance_den_coeffs,-50:.01:50));
% December 5: Good up to here.

% coeffs need to be ordered from greatest to least power
[numz, denz] = Gewertz(1, resistance_den_coeffs);
[components,~] = PoleExtractionAtInfinity(numz,denz);

% Print the component values below
for i = 1:length(components)
    fprintf('----------------------------------------------------\n')
    fprintf('|                                                  |\n')
    if i == length(components)
        if mod(i,2) == 0
            fprintf('| Component %d, Inductor Value (Henries):  %f |\n',i,components{i}(1))
            fprintf('| Resistor Value (Ohms):                  %f |\n',components{i}(2))
            fprintf('|                                                  |\n')
            fprintf('-----------------------------------------------------\n')
        else
            fprintf('| Component %d, Capacitor Value (Farads):  %f |\n',i,components{i}(1))
            fprintf('| Resistor Value (Ohms):                  %f |\n', 1/components{i}(2));
            fprintf('|                                                  |\n')
            fprintf('----------------------------------------------------\n')
        end
    elseif mod(i,2) == 0
        fprintf('| Component %d, Inductor Value (Henries):  %f |\n',i,components{i}(1))
        fprintf('|                                                  |\n')
    elseif mod(i,2) ~= 0
        fprintf('| Component %d,  Capacitor Value (Farads): %f |\n',i,components{i}(1));
        fprintf('|                                                  |\n')
    end
end