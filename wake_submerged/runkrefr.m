clear;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end
%%
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,~,~,Fr2,Re,eps,c0,h,f] = pars_wake;

H = 0; % compare with half-submerged
% Re_list = logspace(2,5,20);
% k_list = [linspace(0.04,1.5,40) linspace(1.6,2.5,10)];
% fr2 = 0.5^2;
% cmodes = nan(4,length(k_list),length(Re_list));
fr2_list = linspace(0.5,8,50).^2;
% k_list = [linspace(0.04,1,50) linspace(1.02,1.8,20) linspace(1.9,2.5,5)];
k_list = linspace(0.04,2.5,120);
cmodes = nan(4,length(k_list),length(fr2_list));
zc = 0.74708299+H;

%%
tic;
parfor ii = 1:length(fr2_list)
    p1 = wSubmerged(N,H,1,1,Re,fr2_list(ii),bflow);
% parfor ii = 1:length(Re_list)
%     p1 = wSubmerged(N,H,1,1,Re_list(ii),fr2,bflow);
    p1.numMeth(method);
    addvar = struct('zL1',zc,'eps',eps);
    csall = findmode(k_list,addvar,H,p1,alg, de_singularize, do_balancing, eig_spectrum, f);
    cmodes(:,:,ii) = csall;
end
toc;

% %% plot Fr
% fr = fr2_list.^(0.5);
% figure;
% [X,Y] = meshgrid(k_list,fr);
% oi1 = squeeze(imag(cmodes(1,:,:))).'.*k_list;
% oi2 = squeeze(imag(cmodes(2,:,:))).'.*k_list;
% oi3 = squeeze(imag(cmodes(3,:,:))).'.*k_list;
% surf(X,Y,oi1);

%% plot Re
figure;
[X,Y] = meshgrid(k_list,sqrt(fr2_list));
% [X,Y] = meshgrid(k_list,Re_list);
oi1 = squeeze(imag(cmodes(1,:,:))).'.*k_list;
oi2 = squeeze(imag(cmodes(2,:,:))).'.*k_list;
oi3 = squeeze(imag(cmodes(3,:,:))).'.*k_list;
surf(X,Y,oi1); hold on;
surf(X,Y,oi2);
surf(X,Y,oi3);
% set(gca,'yscale','log');

%% plot Fr (max eig)
figure;
[X,Y] = meshgrid(k_list,sqrt(fr2_list));
% [X,Y] = meshgrid(k_list,Re_list);
oi1 = squeeze(imag(cmodes(1,:,:))).'.*k_list;
oi2 = squeeze(imag(cmodes(2,:,:))).'.*k_list;
oi3 = squeeze(imag(cmodes(3,:,:))).'.*k_list;
omax = max(oi1,oi2);
omax = max(omax,oi3);
omax(omax==0) = NaN;
pcolor(k_list,sqrt(fr2_list),omax); 
% shading flat;
colorbar;

function csall = findmode(k_list,addvar,H,p1,alg, de_singularize, do_balancing, eig_spectrum, f)
    h = 3*2*pi./real(k_list)+H;
    csall = nan(4,length(k_list));
    for i = 1:length(k_list)
        p1.k = k_list(i); p1.h = h(i);
        oall = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        oall = oall(real(oall)>-50); % Remove the eigenvalues assigned by de-singularizing
        [o, osort] = maxeig(oall);
        if real(o) > 0 && imag(o)>1e-5
            addvar.zL1 = p1.invbf(real(o)/k_list(i));
        end

        c = osort/k_list(i);
        if length(c)<4
            c = [c;zeros(4-length(c),1)];
        end
        csall(:,i) = c(1:4);
    end
end

function [o, w] = maxeig(ev)
    a = 1:length(ev);
    ev_ind = a(abs(ev)<1e+3 & abs(ev)>1e-5 & abs(imag(ev)) > 1e-10);
    [~, ind] = sort(imag(ev(ev_ind)),'descend');
    ev_ind = ev_ind(ind);
    w = ev(ev_ind);
    if isempty(w)
        o = NaN;
    else
        o = w(1);
    end
end