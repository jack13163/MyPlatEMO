%% ��������Ӧ��
function [ eff,scheduleplans ] = fat_2( pop )
    %����ȡֵ��Χ[0,1]
    %���߱���ά��w
    %Ŀ�꺯������M=4
    [popsize, w] = size(pop);   %������Ⱥ�Ĺ�ģ�;��߱�����ά��
    eff = pop(:,1:w);
    scheduleplans = cell(1,popsize);    %���ȼƻ�
    %����������Ⱥ
    for p = 1:popsize    
        %��������
        [RT,DSFR,PIPEFR,FPORDER,c1,c2,TKS,FP] = parameter_setting();
        
        ET = [];              %���й��͹޼���
        DSFET = [0,0,0];   %���������һ�ε����ͽ���ʱ��
        PET = 0;                %�ܵ����һ��ת�˵Ľ���ʱ��(ͣ��Ҳ��¼����)
        %������Ϣ��¼��
        scheduleplan = [];
        %��ȡ����
        x = pop(p,1:w);    
        %ͳ�Ƴ�ʼʱ���͹޵ļ���
        for i=1:size(TKS,1)             
            if TKS(i,3) == 0
                ET = [ET,i];
            end
        end
        %δʵ�����ͼƻ�������������
        UD = [1,2,3];
        %�����жϵ�ǰ״̬�Ƿ���Ե���
        if schedulable(FPORDER, DSFR, PIPEFR, RT, TKS, PET, FP, UD) 
            %��������˳�����ѡ��һ���ǿյĹ��͹�Ϊ��������������
            for i = 1:3     %i����������
                for j = 1:2     %j������ϰ�
                    for k = 1:size(TKS,1)       %k�����͹�
                        if FPORDER(i,j) == TKS(k,2)
                            Temp = DSFET(i);
                            %���㹩�͹޵Ĺ��ͽ���ʱ�䣬����ӵ����͹޶�Ӧ�Ĺ��ͼ�¼����
                            DSFET(i) = DSFET(i) + TKS(k,3) / DSFR(i);%DSFET(i)����������i�����һ�ι��ͽ���ʱ�䣬TKS(k,3)�����͹�k�ĵ�ǰ������DSFR(i)����������i����������
                            scheduleplan = [scheduleplan; i, k, Temp, DSFET(i), TKS(k,2)];
                            TKS(k, 4) = i;   %��¼������
                            TKS(k, 5) = Temp;   %��¼���Ϳ�ʼʱ��
                            TKS(k, 6) = DSFET(i);   %��¼���ͽ���ʱ��
                        end
                    end
                end
            end
            %���ݵ��ȣ������������н�
            [ Flag,x,~,~,~,~,~,~,~,scheduleplan ] = schedule( x,1,ET,UD,PET,DSFET,FP,TKS,scheduleplan,DSFR,PIPEFR,FPORDER,RT );
        else
            fprintf('��ʼ��治����ָ��������ϵͳ���ɵ��ȣ�\n');
        end
        %�������
        scheduleplans{p} = scheduleplan;
        scheduleplan=sortrows(scheduleplan,1);

        if Flag
            %������Ӧ�Ⱥ���
            f1 = gnum(scheduleplan);            %���͹޸���
            f2 = gchange(scheduleplan);        %�������Ĺ��͹��л�����
            f3 = gdmix(scheduleplan, c1);      %�ܵ���ϳɱ�
            f4 = gdimix(scheduleplan, c2);     %�޵׻�ϳɱ�
        else
            f1 = inf;
            f2 = inf;
            f3 = inf;
            f4 = inf;
        end

        eff(p,1:w) = x(1:w);   %����
        eff(p,w+1) = f1;
        eff(p,w+2) = f2;
        eff(p,w+3) = f3;
        eff(p,w+4) = f4;
    end
end

%% ��������
function [RT,DSFR,PIPEFR,FPORDER,c1,c2,TKS,FP] = parameter_setting()
    %% ��������(�������)    
   	RT = 6;                                                %פ��ʱ��
	DSFR = [333.3, 291.7, 625];                  %��������������(�ɼ�������������)
    PIPEFR = 1250;                                   %���鲻ͬ�Ĺܵ���������
    FPORDER = [5,1;6,2;4,3];                     %������������˳��
	c1 = [0 11 12 13 7 15;
            10 0 9 12 13 7;
            13 8 0 7 12 13;
            13 12 7 0 11 12;
            7 13 12 11 0 11;
            15 7 13 12 11 0];           %�ܵ���ϳɱ�
  	c2 = [0 11 12 13 10 15;
            11 0 11 12 13 10;
            12 11 0 10 12 13;
            13 12 10 0 11 12;
            10 13 12 11 0 11;
            15 10 13 12 11 0];          %�޵׻�ϳɱ�
    
	%% ��������(�ɱ����)
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
	];                                        %�͹޳�ʼ��棬������Զ�����������͹��Ϳ�ʼʱ��ͽ���ʱ��

	FP = [25992,49008,90000,0,0,0];        %���ϰ�(�ɼ���ԭ������)
end

%% ���ݵ���(ֻ��Ҫ�ҵ�һ�����н⼴��)
function [ Flag,X,step,ET,UD,PET,DSFET,FP,TKS,SchduleTable ] = schedule( X,step,ET,UD,PET,DSFET,FP,TKS,SchduleTable,DSFR,PIPEFR,FPORDER,RT )
% ֻ��Ҫ���α�����ȥ�������ֲ��ɵ��ȵĵط��ͻ���
% Ȼ����ű�����ֱ��FPΪ��
    Flag = false;
    %�����ֵ
    Old_step = step;
    Old_ET = ET;
    Old_UD = UD;
    Old_PET = PET;
    Old_DSFET = DSFET;
    Old_FP = FP;
    Old_TKS = TKS;
    Old_SchduleTable = SchduleTable;
    
    %��¼�������㼣
    footprint = zeros(1,length(UD)+1);
        
    %��������
    if sum(FP) == 0 || isempty(UD)
        Flag = true;
    elseif step>25   %�������ó���25
        Flag = false;
    else
        %������ȱ���
        while ~all(footprint)&&~Flag&&step<=25%�ҵ��⣬�ͷ���
            %���͹�
            TK_NO=getInt(X(2*step-1), size(ET, 2));
            %������
            DS = UD(getInt(X(2*step), size(UD, 2)-1)+1);
            %�������ж�
            f1 = false;
            if TK_NO==0
                PipeStoptime  = roundn(getPipeStoptime(FPORDER, DSFR, RT, TKS, PET, FP, UD),-6);%�������룬ȡ6λС��
                if PipeStoptime > 0 
                    f1 = true;
                    %ͣ��
                    [SchduleTable,TKS,ET,PET] = stop(ET,PET,TKS,SchduleTable,PipeStoptime,DSFR);
                end
            else
                %�Ե���
                [Flag,PET,TKS,FP,ET,UD,DSFET,SchduleTable] = tryschedule(ET(TK_NO),DS,DSFET,PET,PIPEFR,RT,ET,UD,DSFR,TKS,FP,FPORDER,SchduleTable);
                %�ж��Ե����Ƿ�ɹ�����һ״̬�Ƿ���Ե���
               	if Flag && schedulable(FPORDER, DSFR, PIPEFR, RT, TKS, PET, FP, UD)
                    f1 = true;
                end
            end
            %�ж��Ե����Ƿ�ɹ�����һ״̬�Ƿ���Ե���
            if f1
                step = step + 1;
                [ Flag,X,step,ET,UD,PET,DSFET,FP,TKS,SchduleTable ] = schedule( X,step,ET,UD,PET,DSFET,FP,TKS,SchduleTable,DSFR,PIPEFR,FPORDER,RT );
            else
                Flag = false;
                
                %���ݻع�
                [step,ET,UD,PET,DSFET,FP,TKS,SchduleTable] = rollback(Old_step,Old_ET,Old_UD,Old_PET,Old_DSFET,Old_FP,Old_TKS,Old_SchduleTable);
                
                %���Ĳ��ԣ��������Ƿ����ͣ�ˣ�
                for i = 1:length(footprint)
                    if footprint(i) == 0
                        %footprint�У���1λ����ܵ�ͣ��
                        if i == 1
                            while 0 ~= getInt(X(2*step-1),size(ET, 2))
                                X(2*step-1) = rand;
                            end
                            footprint(1) = 1;%ͣ�˱��
                        elseif isempty(ET)    %�Ѿ����Թ�ͣ�˲��У������޿չ�
                            footprint = ones(size(footprint));
                        else    %���Թ�ͣ�˲��У����пչ�
                            while 0 == getInt(X(2*step-1),size(ET, 2))
                                X(2*step-1) = rand;
                            end
                            while (i-1) ~= getInt(X(2*step),size(UD, 2)-1)+1
                                X(2*step) = rand;
                            end
                            footprint(i) = 1;%���������
                        end
                        break;
                    end
                end
            end
        end
    end
end

%% ���ݻع�
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

%% �������ԭ������
function [total] = gettotal(DSFR,TKS,PET)
	EO = TKS(TKS(:,2)~=inf,:);             %���еķǿ��͹�
   	total = zeros(1,length(DSFR));   
	for i=1:size(EO,1)
        %������
        ds = EO(i, 4);
        %ͳ�Ƹ����������Ŀ������
        if EO(i, 5) < PET && EO(i, 6) > PET
          	total(ds) = total(ds) + EO(i,3) - (PET - EO(i, 5)) * DSFR(ds);
        else
        	total(ds) = total(ds) + EO(i,3);
        end
	end
end

%% ͣ��
function [SchduleTable,TKS,ET,PET] = stop(ET,PET,TKS,SchduleTable,PipeStoptime,DSFR)
	PETOLD = PET;
    tk = 0;
    [feedendtimes,ind] = sort(TKS(:,6));
	%ͣ���ڼ䣬�����������Ƿ����ͽ���
	for i = 1:size(TKS,1)
        if feedendtimes(i) > PET && feedendtimes(i) <= PET + PipeStoptime
            tk = i;
            break;
        end
	end
	%���͹��ͷ�
	if tk ~= 0
        %����ת�˽���ʱ��Ϊ���͹޿��õ�ʱ��
        PET = feedendtimes(tk);
        %��ͣ���ڼ��ͷŵĹ��͹���ӵ�ET��
        ET = [ET, ind(tk)];      %����ǰ�͹���ӵ�ET��
        %���͹���Ϣ
        TKS(ind(tk), 2) = inf;
        TKS(ind(tk), 3) = 0;
        TKS(ind(tk), 4) = 0;
        TKS(ind(tk), 5) = 0;
        TKS(ind(tk), 6) = 0;
    else
        %����ͣ��
        PET = PET + PipeStoptime;
	end
	%����ת�˼�¼������Ҫ��ѡ���͹޹���
	SchduleTable = [SchduleTable; length(DSFR) + 1, 0, PETOLD, PET, 0];
end

%% ��ȫ����
function [V] = getVolume(DS,DSFET,PET,RT,UD,total,DSFR,PIPEFR)
	V = (DSFET(DS) - PET - RT) * PIPEFR;    %����פ��ʱ��Լ��
	VSec = inf;     %��������ԭ��ת�˵İ�ȫ���
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

%% �Ե���
function [Flag,PET,TKS,FP,ET,UD,DSFET,SchduleTable] = tryschedule(TK,DS,DSFET,PET,PIPEFR,RT,ET,UD,DSFR,TKS,FP,FPORDER,SchduleTable)
	%ԭ������
	if FP(FPORDER(DS, 1)) == 0
        COT = FPORDER(DS, 2);      %ԭ������2
	else
        COT = FPORDER(DS, 1);      %ԭ������1
	end
	%���������
	total = gettotal(DSFR,TKS,PET);    
    %��ȫ��ԭ�����
    V = roundn(getVolume(DS,DSFET,PET,RT,UD,total,DSFR,PIPEFR),-6);
	%��������Ե����ض�����ֵ
	if V >= 150
        if TKS(TK, 1) < V
            V = TKS(TK, 1);
        end
        if FP(COT) < V
            V = FP(COT);
        end

        %���ϰ�
        FP(COT) = FP(COT) - V;

        %���͹޼���
        ET(ET==TK) = [];       %ɾ��TK

        %ת���ڼ��ͷŵĹ��͹���ӵ����͹޼���
        for i = 1:size(TKS,1)
            if (TKS(i,6) > PET && TKS(i,6) <= PET + V / PIPEFR)
                %����ǰ�͹���ӵ�ET��
                ET = [ET, i];
                %���͹���Ϣ���
                TKS(i, 2) = inf;
                TKS(i, 3) = 0;
                TKS(i, 4) = 0;
                TKS(i, 5) = 0;
                TKS(i, 6) = 0;
            end
        end

        %�ܵ�ת�˽���ʱ��
        PETOLD = PET;
        PET = PET + V / PIPEFR;

        %��������ʱ��
        DSFETOLD = DSFET(DS);
        DSFET(DS) = DSFETOLD + V / DSFR(DS);

        %���͹޵�״̬��Ϣ
        TKS(TK, 2) = COT;
        TKS(TK, 3) = V;
        TKS(TK, 4) = DS;
        TKS(TK, 5) = DSFETOLD;
        TKS(TK, 6) = DSFET(DS);
        Flag = true;
        
    	%����ת�˼�¼�����ͼ�¼
        SchduleTable = [SchduleTable;  length(DSFR)+1, TK, PETOLD, PET, COT];
    	SchduleTable = [SchduleTable;  DS, TK, DSFETOLD, DSFET(DS), COT];
        %�ж�DS�Ƿ����ͳɹ�
     	if roundn(DSFET(DS),-6) == 240
            UD(UD==DS) = [];        %ɾ��DS
     	end
    else
        Flag = false;
	end
end

%% ���͹��л�����
function [ x ] = gchange( a )
    K = 0;
    
    %��ȡת�˼�¼
    for i=1:size(a,1)
        if a(i,1) >= max(a(:,1))
            K = i;
            break;
        end
    end
    x = K - 1;
end

%% �޵׻�ϴ���
function [ x ] = gdimix( a, c2 )
%c2���޵׻�ϳɱ�
    m2 = zeros(6,6);
    K = 0;
    
    %��ȡת�˼�¼
    for i=1:size(a,1)
        if a(i,1) >= max(a(:,1))
            K = i;
            break;
        end
    end
    K = K - 1;
    record = a(1:K,:);
    record = sortrows(record,2);
    
    %ͳ�Ƹ����͹޵�ʹ�ô���
    A = record(:,2)';
    B = unique(A);
    tks = [];
    for i=1:size(B,2)
        tks(i) = sum(A(:)==B(i));
    end
    
	%�������ͼ�¼������޵׻�ϴ���
    for i=1:size(tks,2)     %i�����͹�
        tmp = 0;
        for j=1:(size(record,1) - 1)      %j�������ͼ�¼
            if record(j,2) == i
                tmp =  tmp + 1;
                if tmp >= tks(i)     %tks(i)���͹�i�е�ʹ�ô���
                    break;
                end
                if record(j,5) ~= record(j+1,5)
                    m2(record(j,5), record(j+1,5)) = m2(record(j,5), record(j+1,5)) + 1;
                end
            end
        end
    end
    x = sum(sum(m2.*c2));          %�޵׻�ϴ���
end

%% ��ȡ�ܵ���ϴ���
function [ x ] = gdmix( a, c1 )
%c1���ܵ���ϳɱ�
    m1 = zeros(6,6);
    K = 0;
    
    %��ȡת�˼�¼
    for i=1:size(a,1)
        if a(i,1) >= max(a(:,1))
            K = i;
            break;
        end
    end
    PIPE = a(K:size(a,1),:);
    PIPE = PIPE(PIPE(:,5)~=0,:);   %ɾ�����е�ͣ�˼�¼
    
    for i=1:size(PIPE,1)-1
        if PIPE(i,5) ~= PIPE(i+1,5)
            m1(PIPE(i,5),PIPE(i+1,5)) = m1(PIPE(i,5),PIPE(i+1,5)) +1;
        end
    end
    x = sum(sum(m1.*c1));       %�ܵ���ϴ���
end

%% ���͹޸���
function [ x ] = gnum( a )
    TKSet = [];
    for i=1:size(a,1)
        if a(i,2) ~= 0 && ~ismember(a(i,2), TKSet)
            TKSet = [TKSet, a(i,2)];
        end
    end
    x = size(TKSet,2);      %���͹޵ĸ���
end

%% ����[0~1]֮��������������һ��[0~length]����
function [ r ] = getInt( a, length )
%% ������������
    persistent squence;
    x0 = 0.29583;
    len = 1000;
    %�ж��Ƿ��Ѿ����ɹ�
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
    %a��[0~1]
    if a == 1
        r = length;
    else
        r = fix((length + 1) * a);
    end
end

%% ��ȡ�ܵ�ͣ�˰�ȫʱ��
function [ t ] = getPipeStoptime(FPORDER, DSFR, RT, TKS, PET, FP, UD)
    t=inf;
    %1.0000  333.3000
    %3.0000  625.0000
    UDFR = sortrows([UD;DSFR(UD)]',2);   %����������������(��������)
    available = zeros(size(DSFR, 2), 2);       %����и���������Ŀǰ���õ�����
    for i = 1:size(UD, 2)
        DSN = UDFR(i,1);    %������
        COTN1 = FPORDER(DSN, 1);      %ԭ������1
        COTN2 = FPORDER(DSN, 2);      %ԭ������1
        for j = 1:size(TKS, 1)          %j�����͹�
            if TKS(j, 2) == COTN1
                if TKS(j, 5) <= PET && TKS(j, 6) > PET
                    available(DSN, 1) = available(DSN, 1) + TKS(j, 3) - (PET - TKS(j, 5)) * DSFR(DSN);
                else
                    available(DSN, 1) = available(DSN, 1) + TKS(j, 3);
                end
            elseif TKS(j, 2) == COTN2
                if TKS(j, 5) <= PET && TKS(j, 6) > PET
                    available(DSN, 2) = available(DSN, 2) + TKS(j, 3) - (PET - TKS(j, 5)) * DSFR(DSN);   %����״̬
                else
                    available(DSN, 2) = available(DSN, 2) + TKS(j, 3);  %�ǹ���״̬
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

    %����Լ��
    K = size(UDFR, 1);
    for i=1:K
        arfai = RT * UDFR(i,2);
        DSN = UDFR(i,1);

        %Ԥ��һ���������ڵ�ԭ��
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

%% �жϵ�ǰϵͳ״̬�Ƿ�ɵ���
function [ flag ] = schedulable(FPORDER, DSFR, PIPEFR, RT, TKS, PET, FP, UD)
    %flag���ɵ���Ϊ1������Ϊ0
    UDFR = sortrows([UD;DSFR(UD)]',2);   %����������������(��������)
    available = zeros(size(DSFR, 2), 2);       %����и���������Ŀǰ���õ�����
    for i = 1:size(UDFR, 1)
        DSN = UDFR(i,1);    %������
        COTN1 = FPORDER(DSN, 1);      %ԭ������1
        COTN2 = FPORDER(DSN, 2);      %ԭ������1
        for j = 1:size(TKS, 1)          %j�����͹�
            if TKS(j, 2) == COTN1
                if TKS(j, 5) <= PET && TKS(j, 6) > PET
                    available(DSN, 1) = available(DSN, 1) + TKS(j, 3) - (PET - TKS(j, 5)) * DSFR(DSN);
                else
                    available(DSN, 1) = available(DSN, 1) + TKS(j, 3);
                end
            elseif TKS(j, 2) == COTN2
                if TKS(j, 5) <= PET && TKS(j, 6) > PET
                    available(DSN, 2) = available(DSN, 2) + TKS(j, 3) - (PET - TKS(j, 5)) * DSFR(DSN);   %����״̬
                else
                    available(DSN, 2) = available(DSN, 2) + TKS(j, 3);  %�ǹ���״̬
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

    %����Լ��
    K = size(UDFR, 1);
    tag = zeros(1,K);
    for i=1:K
        arfai = RT * UDFR(i,2);
        DSN = UDFR(i,1);

        if (i==K && total(DSN) >= 4 * K * arfai) || total(DSN) >= 2 * K * arfai                 %ͳ����һ�����ڲ���Ҫ���͹�����
            tag(i) = 1;
        elseif (i==K && total(DSN) + 1 < 2 * K * arfai ) || total(DSN) + 1 < K * arfai      %�ж��Ƿ������������Լ��
            flag=0;
            %disp('�������Լ�������㣡');
            return;
        end
    end

    %�ڵ��������ڽ������Ͳ���
    time = [];
    cur = PET;
    for i=1:K
        beta = available(i,1);
        arfai = RT * UDFR(i,2);
        %ת�˸�����������Ҫ����ԭ��
        if tag(i) == 0
            %������Ʒת��
            time = [time, cur];
            %��Ʒ�л���Ҫ���͹�
            if beta > 0 && beta < K * arfai && available(i,2) >= K * arfai - beta
                time = [time, cur + beta / PIPEFR];
            end
            cur = cur + K * arfai / PIPEFR;
        end
    end

    %���͹޿���ʱ��ideltime
    TKSN = sortrows(TKS,6);
    ideltime = TKSN(:,6)';

    count = 0;%�����ù��͹޵ĸ���
    for i=1:size(time,2)
        if ideltime(i) > time(i)
            count = count + 1;
        end
    end

    %�жϹ��͹��Ƿ�����Ҫʱ����
    if count == 0
        flag = 1;
    else
        %disp('���͹޲��㣡');
        flag = 0;
    end
end