clear all;
dataObjs = findobj(gcf,'-property','YData');
for i = 1:length(dataObjs)
    x{i} = dataObjs(i).XData;
    y{i} = dataObjs(i).YData;
end

%% Change nondimenional parameters (Morland->Kuo)
lsy = {'-ko',':ro','--bo','-.go'};
hold on;
u = [30 30 80 80]; 
lam = [4 20 4 20];
for i = 1:length(y)
    [ud_nd,delta_nd,lam_nd,Re] = wtoMor(u(i),lam(i));
    ynd{i} = Mortow(ud_nd,delta_nd,lam_nd,Re,y{i});
    dataObjs(i).YData = ynd{i};
end