clear;
load('h02050.mat');
[X,Y] = meshgrid(k_list,fr);
% [X,Y] = meshgrid(k_list,Re_list);
% fr = Re_list;
oi1 = squeeze(imag(cmodes(1,:,:))).'.*k_list;
oi2 = squeeze(imag(cmodes(2,:,:))).'.*k_list;
oi3 = squeeze(imag(cmodes(3,:,:))).'.*k_list;

%%
mf = []; mk = [];
for i = 1:length(fr)
    [~, ind] = findpeaks(oi2(i,1:45));
    mk = [mk ind];
    mf = [mf i*ones(1,length(ind))];
end
mf1 = nan(1,length(fr)); mk1 = nan(1,length(fr)); mo1 = nan(1,length(fr));
for i = 1:1:length(fr)
    ind = mf==i;
    if sum(ind)>1
        modeb = i;
        break;
    end
    mf1(i) = mf(ind);
    mk1(i) = mk(ind);
    mo1(i) = oi2(mf(ind),mk(ind));
end
mf3 = nan(1,length(fr)); mk3 = nan(1,length(fr)); mo3 = nan(1,length(fr));
for i = length(fr):-1:1
    ind = mf==i;
    if sum(ind)>1
        modee = i;
        break;
    end
    mf3(i) = mf(ind);
    mk3(i) = mk(ind);
    mo3(i) = oi2(mf(ind),mk(ind));
end

for j = modeb:modee
    ind = mf==j;
    a = 1:length(mf);
    ind = a(ind);
    mf1(j) = mf(ind(1));
    mk1(j) = mk(ind(1));
    mo1(j) = oi2(mf(ind(1)),mk(ind(1)));
    mf3(j) = mf(ind(2));
    mk3(j) = mk(ind(2));
    mo3(j) = oi2(mf(ind(2)),mk(ind(2)));
end

%%
o1 = -100*ones(length(fr),length(k_list));
o3 = -100*ones(length(fr),length(k_list));
for i = 1:length(fr)
    if isnan(mk1(i))
        o3(i,:) = oi2(i,:);
    elseif isnan(mk3(i))
        o1(i,:) = oi2(i,:);
    else
        [~, mink] = findpeaks(-oi2(i,1:45));
        o1(i,1:mink) = oi2(i,1:mink);
        o1(i,mink+1:end) = oi3(i,mink+1:end);
        o3(i,mink:end) = oi2(i,mink:end);
        o3(i,1:mink-1) = oi3(i,1:mink-1);
    end
end
%%
o4 = oi1;
o4(o4<=0) = -1;
o1(o1<=0) = -1;
o3(o3<=0) = -1;
figure;
surf(X,Y,o1); hold on;
surf(X,Y,o3);
surf(X,Y,o4);
hold off;
zlim([0 1.1*(max(oi1,[],'all'))]);
caxis([0 max(oi1,[],'all')]);
colormap(parula);
colorbar('location','eastoutside');
view([1 -2 1.5]);
xlabel('$k$');
ylabel('$Fr$');
zlabel('$\omega_i$');
%%
omodes = {o1,o3,o4};
for i = 1:length(omodes)
    figure('position',[50 50 1080 720]);
    surf(X,Y,omodes{i});
    zlim([0 1.1*(max(oi1,[],'all'))]);
    caxis([0 max(oi1,[],'all')]);
    colormap(parula);
    colorbar('location','eastoutside');
    view([1 -2 1.5]);
    xlabel('$k$');
    ylabel('$Fr$');
    zlabel('$\omega_i$');
end
%%
omodes = {o1,o3,o4};
for i = 1:length(omodes)
    figure('position',[50 50 1080 720]);
    surf(X,Y,omodes{i});
    zlim([0 1.1*(max(oi1,[],'all'))]);
    caxis([0 max(oi1,[],'all')]);
    colormap(parula);
    colorbar('location','eastoutside');
    view(2);
    xlabel('$k$');
    ylabel('$Fr$');
    zlabel('$\omega_i$');
end
%%
omodes = {o1,o3,o4};
fcolor = {'r','b','#7E2F8E'};
fa = [0.3,0.3,0.5];
figure; hold on;
zlim([0 1.1*(max(oi1,[],'all'))]);
for i = 1:length(omodes)
    surf(X,Y,omodes{i},'FaceColor',fcolor{i}, 'FaceAlpha',fa(i), 'EdgeColor','none');
end
view(2);
xlabel('$k$');
ylabel('$Fr$');
zlabel('$\omega_i$');
hold off; box on;