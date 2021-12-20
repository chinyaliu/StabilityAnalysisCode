clear;
load('Hfr253080.mat')
oi{1} = squeeze(imag(cmodes(1,:,:))).'.*k_list;
oi{2} = squeeze(imag(cmodes(2,:,:))).'.*k_list;
oi{3} = squeeze(imag(cmodes(3,:,:))).'.*k_list;

%% plot H
figure;
[X,Y] = meshgrid(k_list,H_list);
surf(X,Y,oi{1}); hold on;
surf(X,Y,oi{2});
surf(X,Y,oi{3}); hold off;
zlim([0 0.15]);

%%
hold on;
maxh = cell(length(oi),1); minh = cell(length(oi),1);
maxk = cell(length(oi),1); mink = cell(length(oi),1);
maxo = cell(length(oi),1); mino = cell(length(oi),1);
for j = 1:length(oi)
    mh = []; mk = []; mo = [];
    for i = 1:length(H_list)
        [~, ind] = findpeaks(oi{j}(i,:));
        mk = [mk ind];
        mh = [mh i*ones(1,length(ind))];
        mo = [mo (ind-1)*length(H_list)+i];
    end
    maxh{j} = mh;
    maxk{j} = mk;
    maxo{j} = mo;
    
    plot3(k_list(mk),H_list(mh),oi{j}(mo),'o');
    
    mh = []; mk = []; mo = [];
    for i = 1:length(H_list)
        [~, ind] = findpeaks(-oi{j}(i,:));
        mk = [mk ind];
        mh = [mh i*ones(1,length(ind))];
        mo = [mo (ind-1)*length(H_list)+i];
    end
    minh{j} = mh;
    mink{j} = mk;
    mino{j} = mo;
    
    plot3(k_list(mk),H_list(mh),oi{j}(mo),'s');
end