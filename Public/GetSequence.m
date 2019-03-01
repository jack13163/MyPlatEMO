function [ squence ] = GetSequence(x0,len)
% 获取混沌序列
    squence = [];
    %判断是否已经生成
    if isempty(squence)
        squence = zeros(1,len);
        xn = x0;
        for i=1:len
            xn = 4*xn*(1-xn);
            squence(i) = xn;
        end
    end
end
