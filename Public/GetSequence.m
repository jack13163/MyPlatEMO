function [ squence ] = GetSequence(x0,len)
% ��ȡ��������
    squence = [];
    %�ж��Ƿ��Ѿ�����
    if isempty(squence)
        squence = zeros(1,len);
        xn = x0;
        for i=1:len
            xn = 4*xn*(1-xn);
            squence(i) = xn;
        end
    end
end
