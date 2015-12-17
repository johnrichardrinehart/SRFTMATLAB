function [ineq, eq] = CheckDecreasing(item)
% This function is used to ensure the D values are not increasing
if isrow(item)
    item = item';
end
ineq = item; % fmincon will try to keey all of ineq <= 0... so, we're good
%ineq = 0;
ineq(end) = 0; % The end won't make sense
eq = 0;