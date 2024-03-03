% Simulation parameter
%        Emax  e   emat   E0    tnode        c
%Map1    100  .8   1     .3   1 no renew    1000    Finished
%Map2    100  .8   1     .3   1 no renew    1000
%Map3    100  .8   1     .3   1 no renew  200




clear

for kk=1:1
    clearvars -except kk
    
    str='c'
    set=1;
    
    %%
    row=kk+1;  % row value update  in excel sheet
    sheet = 1; % no. of sheet in excel
    NS=1;
    NE=2*row+13;               % shor testpath(G,1,5)
    d_map=NE;   % dimension of road network
    v=2%  # of EV
    
    
    if str=='rc' | str=='c' |str=='rc' 
                load( ['Map1/',str,'10',num2str(set),'/Nodes_PickupTime_',str,num2str(d_map) ])  %  load a  DAG and pickup time

    elseif str=='be'
            load( ['RealSystem_Belgium/real/Nodes_PickupTime_',num2str(d_map) ])  %  load a  DAG and pickup time

    end
    
    
    tic
    %% data
    % matrix,  energy e, time t, distance  d
    d=full(adjacency(G,'weighted'));
    d=d(1:end,1:end);
    N=size(d,1);
    
    
    g1=1/22+0.299; g2=1/3+0.129;  % EV related parameter, i.e. charging rate
    g=[];
    for i=1:size(d,1)
        if rem(i,2)==0
            g=[g,g1];
        elseif rem(i,2)==1
            g=[g,g2];
        end
    end
    g(end)=0;
    
    
    tnode=G.Nodes.PickupTime';%%%%%%% Pickup time
    t=d/90;  % 1/18   3/40
    e=0.24*d;
    Emax=90;
    gpbar=g;
    emat=e;
    tmat=t;
    w1=1;
    c=19;

    cc=9.05*ones(1,N);
    for i=1:N
        if i==NS
            for j=1:N
                Cost(i,j)=c;
            end
        else
            for j=1:N
                Cost(i,j)=-cc(i);
            end
        end
    end
    
       %%  Objective Function
 
    
    
    X = sdpvar(N,N,'full');  % n^2 d.v.
    Z= sum(sum(X.*Cost))+sum(sum(X.*tmat));   %  Travel time+Travel cost
    
    %%  Constraints
    C=[];
    C=[C,0<=X<=1];
    for i=1:N;   %Remove unlinked arcs
        for j=1:N
            if emat(i,j)==0
                C=[C,X(i,j)==0];
            end
        end
    end
    
    %% Flow conservation
    
    for i=1:N
        if i==NS
            C=[C, sum(X(i,:))-sum(X(:,i))<=v];
        elseif i==NE
            C=[C,sum(X(i,:))-sum(X(:,i))<=-v];
        else
            C=[C,sum(X(i,:))-sum(X(:,i))==0];
            C=[C,sum(X(i,:))<=1];    % Prevent the occurence of one node having two flows
        end
    end
    
    %% Time constraints, Fixed pickup time constraints
    for i=1:N
        for j=1:N
            if tmat(i,j)>0
                C=[C,0>=(tnode(i)+tmat(i,j)+g(i)*emat(i,j)-tnode(j))*X(i,j)];  %  Emax
            end
        end
    end
    
    
    %%
    tic
ops = sdpsettings('solver','gurobi');    %cplex
    result  = optimize(C,Z,ops)
    toc
    %%
    
    if result.problem == 0
        dxx=value(X);
        obj=value(Z)
    else
        disp('Something wrong');
    end
    [i,j]=find(dxx==1);
    XY=[i j];
    % path=unique(sort([XY(:,1)' XY(:,2)']))
    %% Visulization
    final=digraph(dxx);
    ServedCustoms=size(unique(final.Edges.EndNodes),1)-2;
    
    % From endnodes to path
    xy=XY;
    for k=1:v
        q=NS;
        route(k,1)=q;
        i=1;
        while i<=size(xy,1)
            index=find(xy(:,1)==q);
            
            if i==1
                k==1:size(index,1);
                k=find(ans==1);
                index=index(k);
            end
            route(k,i+1)=xy(index(1),2);
            q=route(k,i+1);
            if q==NE
                break
            end
            i=i+1;
        end
    end
    
    
    
    %% New road network for individual EV
    % Given traversed path for individual EV, route(k,:)
    for q=1:v
        net=route(q,:);
        net=net(net>0);
        N=size(net,2);
        
        for i=1:N
            for j=1:N
                eMat(i,j)=emat(net(i),net(j));
            end
        end
        %% Chagring problem for Multiple EVs
        
        p=zeros(1,N);
        w2=1;   w3=1;
        r=sdpvar(1,N,'full')
        E=sdpvar(1,N,'full')
        
        Z2=0;
        for i=1:N
            Z2=Z2+(w2*gpbar(i))*r(i);
        end
        
        
        C2=[];
        for  i=1:N-1
            
            C2=[C2, 0<=E(i)<=Emax];
            C2=[C2, E(i+1)==E(i)+r(i+1)-eMat(i,i+1)];
            C2=[C2,0<=r(i)<=eMat(i,1+i) ]; %%%%   Should pay attention in this eq.
            
        end
        C2=[C2, 0<=E(N)<=Emax];
        C2=[C2, E(1)==.7*Emax];
        C2=[C2,0==r(N) ];
        
        
        % C2=[C2,0==r(1) ];
        % C2=[C2,0==r(N) ];
        
        result2(q)  = optimize(C2,Z2);
        if result2(q).problem == 0
            R{q}=value(r);
            En{q}=value(E);
            
            objective(q)=value(Z2);
        else
            break
            disp('Something wrong');
        end
    end

    TotalCostLP=sum(objective)+obj
    
    disp('Total calculation time: ')
    toc

    
    figure(kk)
    plot(final)
    % title(['LP-Served nodes-',num2str(ServedCustoms),'--Solver time-',num2str(sum([result.solvertime result2.solvertime])),'obj-',num2str(TotalCostLP)])%,title(num2str(ServedCustoms))
    % saveas(gcf,['C:\Users\vv\Desktop\',num2str(NS),'-',num2str(NE),'(',num2str(v),')','LP.png'])
    % saveas(gcf,['/Users/vulcanyao/Desktop/',num2str(NS),'-',num2str(NE),'LP.png'])
    
    %%
    % A={'NS','NE','number of vehicle','served customers','solver time','total cost'};
%     A={NS, NE,v,ServedCustoms,sum([result.solvertime result2.solvertime]),TotalCostLP }
%     range = ['H',num2str(row)];
% %         xlswrite(['LP_MIP_version','Belgium.xls'],A,sheet,range);
%        dlmwrite('Test_TLP.csv',A,'-append');

end