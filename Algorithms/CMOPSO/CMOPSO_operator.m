function Offspring = CMOPSO_operator(Global,Population)
% <operator> <real>
% Particle swarm optimization in CMOPSO
% proM ---  1 --- The expectation of number of bits doing mutation 
% disM --- 20 --- The distribution index of polynomial mutation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Begin
    [proM,disM] = Global.ParameterSet(1,20);
    P_Dec = Population.decs;     
    [N,D] = size(P_Dec); 
    P_Obj = Population.objs;
    V     = Population.adds(zeros(N,D));
    Off_P = zeros(N,D);
    Off_V = zeros(N,D);  
    
    %% Get leaders  
    Front     = NDSort(P_Obj,inf);    
    [~,rank]  = sortrows([Front',-CrowdingDistance(P_Obj,Front)']);
    LeaderSet = rank(1:10);
    
    %% Learning
    for i = 1 : N
        % Competition according to the angle 
        winner = LeaderSet(randperm(length(LeaderSet),2));
        c1     = dot(P_Obj(i,:),P_Obj(winner(1),:))/(norm(P_Obj(i,:))*norm(P_Obj(winner(1),:)));
        angle1 = rad2deg(acos(c1));
        c2     = dot(P_Obj(i,:),P_Obj(winner(2),:))/(norm(P_Obj(i,:))*norm(P_Obj(winner(2),:)));
        angle2 = rad2deg(acos(c2));
        mask   = (angle1 > angle2);
        winner = ~mask.*winner(1) + mask.*winner(2);
        % Learning
        r1 = rand(1,D);
        r2 = rand(1,D);
        Off_V(i,:) = r1.*V(i,:) + r2.*(P_Dec(winner,:)-P_Dec(i,:));
        Off_P(i,:) = P_Dec(i,:) + Off_V(i,:);
    end

    %% Polynomial mutation
    Lower = repmat(Global.lower,N,1);
    Upper = repmat(Global.upper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Off_P(temp) = Off_P(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                  (1-(Off_P(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Off_P(temp) = Off_P(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                  (1-(Upper(temp)-Off_P(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
 
   Offspring = INDIVIDUAL(Off_P,Off_V);
end