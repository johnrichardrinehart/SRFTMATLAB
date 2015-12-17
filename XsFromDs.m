% function X = XsFromDs(break_freqs, function_freqs, Ds)
%
% F = @(w,wk)( (w+wk)*log(abs(w+wk)) + (w-wk)*log(abs(w-wk)) );
% % The following bk expressions cover the cases when w = w_k or w = w_{k-1}
% % 1) w ~= w_k or w ~= w_{k-1}
% bk = @(w,ws,widx)( (pi*(ws(widx)-ws(widx-1)))^-1*( F(w,ws(widx)) - F(w,ws(widx-1))) );
%
% X = zeros(1,length(Ds)+1); % Make sure X is zero before we construct it
% X(1) = 0; % bk(0,wk) = 0;
% for j = 2:length(omegas) % j is going to store the omega index
%     w = omegas(j);
%     for k = 1:length(Ds)
%         wk = omegas(k+1);
%         wkm1 = omegas(k);
%         if w == wk
%             bk1 = (pi*(wk-wkm1))^-1*( ...
%                 2*wk*log(2*wk) - ...
%                 (wk + wkm1)*log(wk + wkm1) - ...
%                 (wk-wkm1)*log(wk-wkm1) ) ;
%             X(j) = X(j) + Ds(k)*bk1;
%         elseif w == wkm1 % 4th term is zero
%             wkm2 = omegas(k-1);
%             bk2 = (pi*(wkm1-wkm2))^-1*( ...
%                 2*wkm1*log(2*wkm1) - ...
%                 (wkm1+wkm2)*log(wkm1+wkm2) - ...
%                 (wkm1-wkm2)*log(wkm1-wkm2) ) ;
%             X(j) = X(j) + Ds(k)*bk2;
%         else
%             bk4 = bk(w,omegas,k+1);
%             X(j) = X(j)+Ds(k)*bk4;
%         end
%     end
% end

function X = XsFromDs(break_freqs, function_freqs, Ds)
X = zeros(length(function_freqs),1);
bs = BsFromBreakFreqs(break_freqs, function_freqs);
for i = 1:length(Ds)
    X = X + Ds(i)*bs{i};
end
X(1)=0;
end

function b = BsFromBreakFreqs(break_freqs, function_freqs)
b = cell(length(break_freqs)-1,1);
for i = 1:length(break_freqs)-1
    b{i,1} = zeros(length(function_freqs),1);
       wk = break_freqs(i+1);
       wkm1 = break_freqs(i);
    for j = 1:length(function_freqs)
        w = function_freqs(j);
        if w ~=wk && w~=wkm1
            b{i,1}(j,1) = (pi*(wk-wkm1))^-1*( ...
                (w+wk)*log(abs(w+wk)) + ...
                (w-wk)*log(abs(w-wk)) - ...
                (w+wkm1)*log(abs(w+wkm1)) - ...
                (w-wkm1)*log(abs(w-wkm1)) );
        elseif w == wkm1
            b{i,1}(j,1) = (pi*(wk-wkm1))^-1*( ...
                (w+wk)*log(abs(w+wk)) + ...
                (w-wk)*log(abs(w-wk)) - ...
                (w+wkm1)*log(abs(w+wkm1)) );
        elseif w == wk
            b{i,1}(j,1) = (pi*(wk-wkm1))^-1*( ...
                (w+wk)*log(abs(w+wk)) - ...
                (w+wkm1)*log(abs(w+wkm1)) - ...
                (w-wkm1)*log(abs(w-wkm1)) );
        end
    end
end
end