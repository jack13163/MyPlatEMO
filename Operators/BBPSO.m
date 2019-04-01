function [ Offspring ] = BBPSO( Global,Parent )
% <operator> <real>
% 混沌骨干粒子群算法

    %% 产生混沌序列
    persistent squence;
    x0 = 0.5546;
    len = 1000;
    %判断是否已经生成过
    if isempty(squence)
        squence = zeros(1,len);
        xn = x0;
        for i=1:len
            xn = 4*xn*(1-xn);
            squence(i) = xn;
        end
    end
    
    %% 取出参数：Population,Pbest,Gbest
    n = length(Parent);
    Particles      = Parent([1:end,1:ceil(end/3)*3-end]);
    PopulationDec = Particles(1:n/3).decs;
    Pbest   = Particles(n/3+1:n/3*2).decs;
    Gbest   = Particles(n/3*2+1:end).decs;
    
	%% 按照常规进行搜索，产生种群
	for i = 1:n/3
        %粒子位置移动
        for j = 1:Global.D
            c1 = rand;
            c2 = 1 - c1;
            if rand < 0.5
                a = c1 * Pbest(i,j) + c2 * Gbest(i,j);
                b = abs(Pbest(i,j) - Gbest(i,j));
                PopulationDec(i,j) = a + b * randn;     %正态分布
            else
              	PopulationDec(i,j) = Pbest(i,j);
            end
         	%越界处理
            if PopulationDec(i,j) < Global.lower(j)
                PopulationDec(i,j) = Global.lower(j);
          	elseif PopulationDec(i,j)  > Global.upper(j)
                PopulationDec(i,j)  = Global.upper(j);
            end
        end
	end
    
    Offspring = INDIVIDUAL(PopulationDec);
end