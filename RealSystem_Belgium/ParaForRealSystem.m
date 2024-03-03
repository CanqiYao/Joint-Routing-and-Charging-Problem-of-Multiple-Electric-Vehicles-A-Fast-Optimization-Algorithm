

clear
kmax=10;
k=1;
while k<=kmax
    clearvars -except k*

load('RawData.mat')
tmin=tpick;



nn=15+k*2;   % Randomly select nn.  nodes w.r.t.  time constraints
nodes=randperm(100,nn);
tmin=sort(tmin(nodes));

for i=1:nn
    for j=i:nn
        dnew(i,j)=d(nodes(i),nodes(j));
    end
end


G=digraph(dnew);

G.Nodes.PickupTime=tmin';   % Given set,  B


save(['real','/Nodes_PickupTime_',num2str(nn)],'G')

% save([num2str(str),num2str(set),'\Nodes_PickupTime_',str,num2str(nn)],'G')

k=k+1;
end


%% Insights



% figure(1)
% plot(G)
% figure(2)
% plot(x,y,'o')
% xlabel('x')
% ylabel('y')
% title('City map with 50 nodes and 1225 edges ')
% saveas(gcf,'map.jpg')




%%  Insights
% nn=size(G.Nodes,1);
% tmin=G.Nodes.PickupTime;
% x=G.Nodes.x;
% y=G.Nodes.y;
% 
% 
% for i=1:nn
%     for j=i:nn
%         d(i,j)=hypot(x(i)-x(j),y(i)-y(j));
%     end
% end
% d=triu(d,1);
% tmat=d/18;  %
% emat=.23*d;   %  23kmh/km 
% 
% 
% var(var(tmat))



% for i=1:nn-1
% deltaT(i)=tmin(i+1)-tmin(i);    
% end
% figure(3)
% bar(deltaT)
% hold on
% mean(deltaT)
% % title('')
% for j=1:nn
% vv(j)=mean(tmat(j,:));
% end
% figure
% bar(vv)
% legend('Given ','actual_{max}')