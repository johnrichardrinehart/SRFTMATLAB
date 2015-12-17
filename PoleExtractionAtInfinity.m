function [q, r] = PoleExtractionAtInfinity(NumPolynomial,DenomPolynomial)
NumPolynomial = RemoveSmallValues(NumPolynomial);
DenomPolynomial = RemoveSmallValues(DenomPolynomial);
NumPolynomial = RemoveLeadingZeros(NumPolynomial);
DenomPolynomial = RemoveLeadingZeros(DenomPolynomial);
el = 1;
[q{el}, r{el}] = deconv(DenomPolynomial,NumPolynomial); %Assuming den is of greater  or equal order than num
q{el} = RemoveLeadingZeros(q{el}); r{el} = RemoveLeadingZeros(r{el});
el = el + 1;
[q{el}, r{el}] = deconv(NumPolynomial,r{el-1});
q{el} = RemoveLeadingZeros(q{el}); r{el} = RemoveLeadingZeros(r{el});
while any(r{end})
    el = el + 1;
    r{el-2} = RemoveLeadingZeros(r{el-2});
    r{el-1} = RemoveLeadingZeros(r{el-1})   ;
    [q{el}, r{el}] = deconv(r{el-2},r{el-1});
end
end

function item = RemoveLeadingZeros(item)
zerosidx = find(item,1,'first');
item(1:zerosidx-1) = [];
end

function item = RemoveSmallValues(item)
item = item(abs(item) > 10^-8);
end