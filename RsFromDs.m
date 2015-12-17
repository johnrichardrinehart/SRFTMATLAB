% function R = RsFromDs(optim_freqs,Ds,R0)
% R = zeros(1,length(Ds)+1);
% R(1) = R0;
% for i = 1:length(Ds)
%     R(i+1) = R(i) + Ds(i);
% end
% end

function R = RsFromDs(break_freqs, function_freqs, Ds, R0)
R = R0*ones(length(function_freqs),1);
a = AsFromBreakFreqs(break_freqs,function_freqs);
for i = 1:length(Ds)
    R = R + Ds(i)*a{i,1}(:);
end
end

function a = AsFromBreakFreqs(break_freqs,function_freqs)
a = cell(length(break_freqs)-1,1);
for i = 1:length(break_freqs)-1
    a{i,1} = zeros(length(function_freqs),1);
    for j = 1:length(function_freqs)
        w = function_freqs(j);
        wk = break_freqs(i+1);
        wkm1 = break_freqs(i);
        if w >= wk
            a{i,1}(j) = 1;
        elseif w >=wkm1 && w < wk
            a{i,1}(j) = (w - wkm1)/(wk-wkm1);
        else
            a{i,1}(j) = 0;
        end
    end
end
end