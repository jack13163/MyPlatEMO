%% 计算解的适应度
function [ eff,scheduleplans ] = fat_2( pop )
    %变量取值范围[0,1]
    %决策变量维数w
    %目标函数个数M=4
    [popsize, w] = size(pop);   %计算种群的规模和决策变量的维数
    eff = pop(:,1:w);
    scheduleplans = cell(1,popsize);    %调度计划
    %遍历粒子种群
    for p = 1:popsize    
        %参数设置
        [RT,DSFR,PIPEFR,FPORDER,c1,c2,TKS,FP] = parameter_setting();
        
        ET = [];              %空闲供油罐集合
        DSFET = [0,0,0];   %蒸馏塔最后一次的炼油结束时间
        PET = 0;                %管道最后一次转运的结束时间(停运也记录在内)
        %调度信息记录表
        scheduleplan = [];
        %获取粒子
        x = pop(p,1:w);    
        %统计初始时空油罐的集合
        for i=1:size(TKS,1)             
            if TKS(i,3) == 0
                ET = [ET,i];
            end
        end
        %未实现炼油计划的蒸馏塔集合
        UD = [1,2,3];
        %首先判断当前状态是否可以调度
        if schedulable(FPORDER, DSFR, PIPEFR, RT, TKS, PET, FP, UD) 
            %根据炼油顺序随机选择一个非空的供油罐为各个蒸馏塔供油
            for i = 1:3     %i代表蒸馏塔
                for j = 1:2     %j代表进料包
                    for k = 1:size(TKS,1)       %k代表油罐
                        if FPORDER(i,j) == TKS(k,2)
                            Temp = DSFET(i);
                            %计算供油罐的供油结束时间，并添加到供油罐对应的供油记录表中
                            DSFET(i) = DSFET(i) + TKS(k,3) / DSFR(i);%DSFET(i)代表蒸馏塔i的最后一次供油结束时间，TKS(k,3)代表供油罐k的当前油量，DSFR(i)代表蒸馏塔i的炼油速率
                            scheduleplan = [scheduleplan; i, k, Temp, DSFET(i), TKS(k,2)];
                            TKS(k, 4) = i;   %记录蒸馏塔
                            TKS(k, 5) = Temp;   %记录供油开始时间
                            TKS(k, 6) = DSFET(i);   %记录供油结束时间
                        end
                    end
                end
            end
            %回溯调度，并修正不可行解
            [ Flag,x,~,~,~,~,~,~,~,scheduleplan ] = schedule( x,1,ET,UD,PET,DSFET,FP,TKS,scheduleplan,DSFR,PIPEFR,FPORDER,RT );
        else
            fprintf('初始库存不满足指定条件，系统不可调度！\n');
        end
        %输出方案
        scheduleplans{p} = scheduleplan;
        scheduleplan=sortrows(scheduleplan,1);

        if Flag
            %计算适应度函数
            f1 = gnum(scheduleplan);            %供油罐个数
            f2 = gchange(scheduleplan);        %蒸馏塔的供油罐切换次数
            f3 = gdmix(scheduleplan, c1);      %管道混合成本
            f4 = gdimix(scheduleplan, c2);     %罐底混合成本
        else
            f1 = inf;
            f2 = inf;
            f3 = inf;
            f4 = inf;
        end

        eff(p,1:w) = x(1:w);   %更新
        eff(p,w+1) = f1;
        eff(p,w+2) = f2;
        eff(p,w+3) = f3;
        eff(p,w+4) = f4;
    end
end

%% 参数设置
function [RT,DSFR,PIPEFR,FPORDER,c1,c2,TKS,FP] = parameter_setting()
    %% 参数设置(不变参数)    
   	RT = 6;                                                %驻留时间
	DSFR = [333.3, 291.7, 625];                  %蒸馏塔炼油速率(可计算蒸馏塔个数)
    PIPEFR = 1250;                                   %几组不同的管道输油速率
    FPORDER = [5,1;6,2;4,3];                     %蒸馏塔的炼油顺序
	c1 = [0 11 12 13 7 15;
            10 0 9 12 13 7;
            13 8 0 7 12 13;
            13 12 7 0 11 12;
            7 13 12 11 0 11;
            15 7 13 12 11 0];           %管道混合成本
  	c2 = [0 11 12 13 10 15;
            11 0 11 12 13 10;
            12 11 0 10 12 13;
            13 12 10 0 11 12;
            10 13 12 11 0 11;
            15 10 13 12 11 0];          %罐底混合成本
    
	%% 参数设置(可变参数)
	TKS = [
            16000,5,8000,0,0,0;
            34000,5,30000,0,0,0;
            34000,4,30000,0,0,0;
            34000,inf,0,0,0,0;
            34000,3,30000,0,0,0;
            16000,1,16000,0,0,0;
            20000,6,16000,0,0,0;
            16000,6,5000,0,0,0;
            16000,inf,0,0,0,0;
            30000,inf,0,0,0,0
	];                                        %油罐初始库存，后面会自动添加蒸馏塔和供油开始时间和结束时间

	FP = [25992,49008,90000,0,0,0];        %进料包(可计算原油种类)
end

%% 回溯调度(只需要找到一个可行解即可)
function [ Flag,X,step,ET,UD,PET,DSFET,FP,TKS,SchduleTable ] = schedule( X,step,ET,UD,PET,DSFET,FP,TKS,SchduleTable,DSFR,PIPEFR,FPORDER,RT )
% 只需要依次遍历下去，若出现不可调度的地方就回溯
% 然后接着遍历，直到FP为空
    Flag = false;
    %保存旧值
    Old_step = step;
    Old_ET = ET;
    Old_UD = UD;
    Old_PET = PET;
    Old_DSFET = DSFET;
    Old_FP = FP;
    Old_TKS = TKS;
    Old_SchduleTable = SchduleTable;
    
    %记录遍历的足迹
    footprint = zeros(1,length(UD)+1);
        
    %结束条件
    if sum(FP) == 0 || isempty(UD)
        Flag = true;
    elseif step>25   %步长不得超过25
        Flag = false;
    else
        %深度优先遍历
        while ~all(footprint)&&~Flag&&step<=25%找到解，就返回
            %供油罐
            TK_NO=getInt(X(2*step-1), size(ET, 2));
            %蒸馏塔
            DS = UD(getInt(X(2*step), size(UD, 2)-1)+1);
            %可行性判断
            f1 = false;
            if TK_NO==0
                PipeStoptime  = roundn(getPipeStoptime(FPORDER, DSFR, RT, TKS, PET, FP, UD),-6);%四舍五入，取6位小数
                if PipeStoptime > 0 
                    f1 = true;
                    %停运
                    [SchduleTable,TKS,ET,PET] = stop(ET,PET,TKS,SchduleTable,PipeStoptime,DSFR);
                end
            else
                %试调度
                [Flag,PET,TKS,FP,ET,UD,DSFET,SchduleTable] = tryschedule(ET(TK_NO),DS,DSFET,PET,PIPEFR,RT,ET,UD,DSFR,TKS,FP,FPORDER,SchduleTable);
                %判断试调度是否成功，下一状态是否可以调度
               	if Flag && schedulable(FPORDER, DSFR, PIPEFR, RT, TKS, PET, FP, UD)
                    f1 = true;
                end
            end
            %判断试调度是否成功，下一状态是否可以调度
            if f1
                step = step + 1;
                [ Flag,X,step,ET,UD,PET,DSFET,FP,TKS,SchduleTable ] = schedule( X,step,ET,UD,PET,DSFET,FP,TKS,SchduleTable,DSFR,PIPEFR,FPORDER,RT );
            else
                Flag = false;
                
                %数据回滚
                [step,ET,UD,PET,DSFET,FP,TKS,SchduleTable] = rollback(Old_step,Old_ET,Old_UD,Old_PET,Old_DSFET,Old_FP,Old_TKS,Old_SchduleTable);
                
                %更改策略（先试试是否可以停运）
                for i = 1:length(footprint)
                    if footprint(i) == 0
                        %footprint中，第1位代表管道停运
                        if i == 1
                            while 0 ~= getInt(X(2*step-1),size(ET, 2))
                                X(2*step-1) = rand;
                            end
                            footprint(1) = 1;%停运标记
                        elseif isempty(ET)    %已经尝试过停运不行，并且无空罐
                            footprint = ones(size(footprint));
                        else    %尝试过停运不行，但有空罐
                            while 0 == getInt(X(2*step-1),size(ET, 2))
                                X(2*step-1) = rand;
                            end
                            while (i-1) ~= getInt(X(2*step),size(UD, 2)-1)+1
                                X(2*step) = rand;
                            end
                            footprint(i) = 1;%蒸馏塔标记
                        end
                        break;
                    end
                end
            end
        end
    end
end

%% 数据回滚
function [step,ET,UD,PET,DSFET,FP,TKS,SchduleTable] = rollback(Old_step,Old_ET,Old_UD,Old_PET,Old_DSFET,Old_FP,Old_TKS,Old_SchduleTable)
        step = Old_step;
        ET = Old_ET;
        UD = Old_UD;
        PET = Old_PET;
        DSFET = Old_DSFET;
        FP = Old_FP;
        TKS = Old_TKS;
        SchduleTable = Old_SchduleTable;
end

%% 计算可用原油总量
function [total] = gettotal(DSFR,TKS,PET)
	EO = TKS(TKS(:,2)~=inf,:);             %所有的非空油罐
   	total = zeros(1,length(DSFR));   
	for i=1:size(EO,1)
        %蒸馏塔
        ds = EO(i, 4);
        %统计各个蒸馏塔的库存油量
        if EO(i, 5) < PET && EO(i, 6) > PET
          	total(ds) = total(ds) + EO(i,3) - (PET - EO(i, 5)) * DSFR(ds);
        else
        	total(ds) = total(ds) + EO(i,3);
        end
	end
end

%% 停运
function [SchduleTable,TKS,ET,PET] = stop(ET,PET,TKS,SchduleTable,PipeStoptime,DSFR)
	PETOLD = PET;
    tk = 0;
    [feedendtimes,ind] = sort(TKS(:,6));
	%停运期间，各个蒸馏塔是否炼油结束
	for i = 1:size(TKS,1)
        if feedendtimes(i) > PET && feedendtimes(i) <= PET + PipeStoptime
            tk = i;
            break;
        end
	end
	%有油罐释放
	if tk ~= 0
        %更新转运结束时间为供油罐可用的时间
        PET = feedendtimes(tk);
        %将停运期间释放的供油罐添加到ET中
        ET = [ET, ind(tk)];      %将当前油罐添加到ET中
        %供油罐信息
        TKS(ind(tk), 2) = inf;
        TKS(ind(tk), 3) = 0;
        TKS(ind(tk), 4) = 0;
        TKS(ind(tk), 5) = 0;
        TKS(ind(tk), 6) = 0;
    else
        %正常停运
        PET = PET + PipeStoptime;
	end
	%更新转运记录表，不需要挑选供油罐供油
	SchduleTable = [SchduleTable; length(DSFR) + 1, 0, PETOLD, PET, 0];
end

%% 安全油量
function [V] = getVolume(DS,DSFET,PET,RT,UD,total,DSFR,PIPEFR)
	V = (DSFET(DS) - PET - RT) * PIPEFR;    %满足驻留时间约束
	VSec = inf;     %对于其他原油转运的安全体积
	for i=1:size(UD,2)
        if UD(i) ~= DS
            if DSFR(UD(i)) == max(DSFR)
                mincapacity = 2 * size(UD, 2) * RT * DSFR(UD(i));
            else
                mincapacity = size(UD, 2) * RT * DSFR(UD(i));
            end
            if VSec > PIPEFR / DSFR(UD(i)) * (total(UD(i)) - mincapacity)
                VSec = PIPEFR / DSFR(UD(i)) * (total(UD(i)) - mincapacity);
            end
        end
	end
	if VSec < V
        V = VSec;
	end
end

%% 试调度
function [Flag,PET,TKS,FP,ET,UD,DSFET,SchduleTable] = tryschedule(TK,DS,DSFET,PET,PIPEFR,RT,ET,UD,DSFR,TKS,FP,FPORDER,SchduleTable)
	%原油类型
	if FP(FPORDER(DS, 1)) == 0
        COT = FPORDER(DS, 2);      %原油类型2
	else
        COT = FPORDER(DS, 1);      %原油类型1
	end
	%库存总油量
	total = gettotal(DSFR,TKS,PET);    
    %安全的原油体积
    V = roundn(getVolume(DS,DSFET,PET,RT,UD,total,DSFR,PIPEFR),-6);
	%体积不可以低于特定的阈值
	if V >= 150
        if TKS(TK, 1) < V
            V = TKS(TK, 1);
        end
        if FP(COT) < V
            V = FP(COT);
        end

        %进料包
        FP(COT) = FP(COT) - V;

        %供油罐集合
        ET(ET==TK) = [];       %删除TK

        %转运期间释放的供油罐添加到供油罐集合
        for i = 1:size(TKS,1)
            if (TKS(i,6) > PET && TKS(i,6) <= PET + V / PIPEFR)
                %将当前油罐添加到ET中
                ET = [ET, i];
                %供油罐信息清空
                TKS(i, 2) = inf;
                TKS(i, 3) = 0;
                TKS(i, 4) = 0;
                TKS(i, 5) = 0;
                TKS(i, 6) = 0;
            end
        end

        %管道转运结束时间
        PETOLD = PET;
        PET = PET + V / PIPEFR;

        %结束供油时间
        DSFETOLD = DSFET(DS);
        DSFET(DS) = DSFETOLD + V / DSFR(DS);

        %供油罐的状态信息
        TKS(TK, 2) = COT;
        TKS(TK, 3) = V;
        TKS(TK, 4) = DS;
        TKS(TK, 5) = DSFETOLD;
        TKS(TK, 6) = DSFET(DS);
        Flag = true;
        
    	%更新转运记录和炼油记录
        SchduleTable = [SchduleTable;  length(DSFR)+1, TK, PETOLD, PET, COT];
    	SchduleTable = [SchduleTable;  DS, TK, DSFETOLD, DSFET(DS), COT];
        %判断DS是否炼油成功
     	if roundn(DSFET(DS),-6) == 240
            UD(UD==DS) = [];        %删除DS
     	end
    else
        Flag = false;
	end
end

%% 供油罐切换次数
function [ x ] = gchange( a )
    K = 0;
    
    %获取转运记录
    for i=1:size(a,1)
        if a(i,1) >= max(a(:,1))
            K = i;
            break;
        end
    end
    x = K - 1;
end

%% 罐底混合次数
function [ x ] = gdimix( a, c2 )
%c2：罐底混合成本
    m2 = zeros(6,6);
    K = 0;
    
    %获取转运记录
    for i=1:size(a,1)
        if a(i,1) >= max(a(:,1))
            K = i;
            break;
        end
    end
    K = K - 1;
    record = a(1:K,:);
    record = sortrows(record,2);
    
    %统计各个油罐的使用次数
    A = record(:,2)';
    B = unique(A);
    tks = [];
    for i=1:size(B,2)
        tks(i) = sum(A(:)==B(i));
    end
    
	%根据炼油记录表分析罐底混合次数
    for i=1:size(tks,2)     %i代表油罐
        tmp = 0;
        for j=1:(size(record,1) - 1)      %j代表炼油记录
            if record(j,2) == i
                tmp =  tmp + 1;
                if tmp >= tks(i)     %tks(i)：油罐i中的使用次数
                    break;
                end
                if record(j,5) ~= record(j+1,5)
                    m2(record(j,5), record(j+1,5)) = m2(record(j,5), record(j+1,5)) + 1;
                end
            end
        end
    end
    x = sum(sum(m2.*c2));          %罐底混合次数
end

%% 获取管道混合次数
function [ x ] = gdmix( a, c1 )
%c1：管道混合成本
    m1 = zeros(6,6);
    K = 0;
    
    %获取转运记录
    for i=1:size(a,1)
        if a(i,1) >= max(a(:,1))
            K = i;
            break;
        end
    end
    PIPE = a(K:size(a,1),:);
    PIPE = PIPE(PIPE(:,5)~=0,:);   %删除其中的停运记录
    
    for i=1:size(PIPE,1)-1
        if PIPE(i,5) ~= PIPE(i+1,5)
            m1(PIPE(i,5),PIPE(i+1,5)) = m1(PIPE(i,5),PIPE(i+1,5)) +1;
        end
    end
    x = sum(sum(m1.*c1));       %管道混合次数
end

%% 供油罐个数
function [ x ] = gnum( a )
    TKSet = [];
    for i=1:size(a,1)
        if a(i,2) ~= 0 && ~ismember(a(i,2), TKSet)
            TKSet = [TKSet, a(i,2)];
        end
    end
    x = size(TKSet,2);      %供油罐的个数
end

%% 根据[0~1]之间的随机数，返回一个[0~length]整数
function [ r ] = getInt( a, length )
%% 产生混沌序列
    persistent squence;
    x0 = 0.29583;
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
    if a ~= 0
        a = squence(ceil(a * len));
    end
    %a：[0~1]
    if a == 1
        r = length;
    else
        r = fix((length + 1) * a);
    end
end

%% 获取管道停运安全时间
function [ t ] = getPipeStoptime(FPORDER, DSFR, RT, TKS, PET, FP, UD)
    t=inf;
    %1.0000  333.3000
    %3.0000  625.0000
    UDFR = sortrows([UD;DSFR(UD)]',2);   %蒸馏塔的炼油速率(升序排列)
    available = zeros(size(DSFR, 2), 2);       %库存中各个蒸馏塔目前可用的油量
    for i = 1:size(UD, 2)
        DSN = UDFR(i,1);    %蒸馏塔
        COTN1 = FPORDER(DSN, 1);      %原油类型1
        COTN2 = FPORDER(DSN, 2);      %原油类型1
        for j = 1:size(TKS, 1)          %j代表供油罐
            if TKS(j, 2) == COTN1
                if TKS(j, 5) <= PET && TKS(j, 6) > PET
                    available(DSN, 1) = available(DSN, 1) + TKS(j, 3) - (PET - TKS(j, 5)) * DSFR(DSN);
                else
                    available(DSN, 1) = available(DSN, 1) + TKS(j, 3);
                end
            elseif TKS(j, 2) == COTN2
                if TKS(j, 5) <= PET && TKS(j, 6) > PET
                    available(DSN, 2) = available(DSN, 2) + TKS(j, 3) - (PET - TKS(j, 5)) * DSFR(DSN);   %供油状态
                else
                    available(DSN, 2) = available(DSN, 2) + TKS(j, 3);  %非供油状态
                end
            end
        end
    end

    total = zeros(1, size(FPORDER, 1));
    for i=1:size(available, 1)
        if FP(FPORDER(i, 1)) ~= 0
            total(i) = available(i,1);
        else
            total(i) = available(i,1) + available(i,2);
        end
    end

    %油量约束
    K = size(UDFR, 1);
    for i=1:K
        arfai = RT * UDFR(i,2);
        DSN = UDFR(i,1);

        %预留一个调度周期的原油
        if i==K
            total(DSN) = total(DSN) - 2 * K * arfai;
        else
            total(DSN) = total(DSN) - K * arfai;
        end

        tmp = total(DSN) / UDFR(i, 2);
        if tmp < t
            t = tmp;
        end
    end
end

%% 判断当前系统状态是否可调度
function [ flag ] = schedulable(FPORDER, DSFR, PIPEFR, RT, TKS, PET, FP, UD)
    %flag：可调度为1，否则为0
    UDFR = sortrows([UD;DSFR(UD)]',2);   %蒸馏塔的炼油速率(升序排列)
    available = zeros(size(DSFR, 2), 2);       %库存中各个蒸馏塔目前可用的油量
    for i = 1:size(UDFR, 1)
        DSN = UDFR(i,1);    %蒸馏塔
        COTN1 = FPORDER(DSN, 1);      %原油类型1
        COTN2 = FPORDER(DSN, 2);      %原油类型1
        for j = 1:size(TKS, 1)          %j代表供油罐
            if TKS(j, 2) == COTN1
                if TKS(j, 5) <= PET && TKS(j, 6) > PET
                    available(DSN, 1) = available(DSN, 1) + TKS(j, 3) - (PET - TKS(j, 5)) * DSFR(DSN);
                else
                    available(DSN, 1) = available(DSN, 1) + TKS(j, 3);
                end
            elseif TKS(j, 2) == COTN2
                if TKS(j, 5) <= PET && TKS(j, 6) > PET
                    available(DSN, 2) = available(DSN, 2) + TKS(j, 3) - (PET - TKS(j, 5)) * DSFR(DSN);   %供油状态
                else
                    available(DSN, 2) = available(DSN, 2) + TKS(j, 3);  %非供油状态
                end
            end
        end
    end

    total = zeros(1, size(FPORDER, 1));
    for i=1:size(available, 1)
        if FP(FPORDER(i, 1)) ~= 0
            total(i) = available(i,1);
        else
            total(i) = available(i,1) + available(i,2);
        end
    end

    %油量约束
    K = size(UDFR, 1);
    tag = zeros(1,K);
    for i=1:K
        arfai = RT * UDFR(i,2);
        DSN = UDFR(i,1);

        if (i==K && total(DSN) >= 4 * K * arfai) || total(DSN) >= 2 * K * arfai                 %统计下一个周期不需要的油罐数量
            tag(i) = 1;
        elseif (i==K && total(DSN) + 1 < 2 * K * arfai ) || total(DSN) + 1 < K * arfai      %判断是否满足最低油量约束
            flag=0;
            %disp('最低油量约束不满足！');
            return;
        end
    end

    %在调度周期内进行运油操作
    time = [];
    cur = PET;
    for i=1:K
        beta = available(i,1);
        arfai = RT * UDFR(i,2);
        %转运该蒸馏塔所需要炼的原油
        if tag(i) == 0
            %正常油品转运
            time = [time, cur];
            %油品切换需要供油罐
            if beta > 0 && beta < K * arfai && available(i,2) >= K * arfai - beta
                time = [time, cur + beta / PIPEFR];
            end
            cur = cur + K * arfai / PIPEFR;
        end
    end

    %供油罐可用时间ideltime
    TKSN = sortrows(TKS,6);
    ideltime = TKSN(:,6)';

    count = 0;%不可用供油罐的个数
    for i=1:size(time,2)
        if ideltime(i) > time(i)
            count = count + 1;
        end
    end

    %判断供油罐是否在需要时可用
    if count == 0
        flag = 1;
    else
        %disp('供油罐不足！');
        flag = 0;
    end
end