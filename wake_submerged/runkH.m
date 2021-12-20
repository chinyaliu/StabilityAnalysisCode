clear;
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,~,~,Fr2,Re,eps,c0,h,f] = pars_wake;

H_list = linspace(0,4,30);
% H_list = [linspace(0,5,30) linspace(6,10,5)];
% k_list = [linspace(0.04,1.5,60) linspace(1.53,2.5,20)];
k_list = linspace(0.04,2.5,80);
cmodes = nan(4,length(k_list),length(H_list));
zc = 0.74708299+H_list;

%%
tic;
parfor ii = 1:length(H_list)
    p1 = wSubmerged(N,H_list(ii),1,1,Re,Fr2);
    p1.numMeth(method);
    csall = findmode(k_list,zc(ii),p1,alg, de_singularize, do_balancing, eig_spectrum, f);
    cmodes(:,:,ii) = csall;
end
toc;

%% plot H
figure;
[X,Y] = meshgrid(k_list,H_list);
oi1 = squeeze(imag(cmodes(1,:,:))).'.*k_list;
oi2 = squeeze(imag(cmodes(2,:,:))).'.*k_list;
oi3 = squeeze(imag(cmodes(3,:,:))).'.*k_list;
surf(X,Y,oi1); hold on;
surf(X,Y,oi2);
surf(X,Y,oi3); hold off;
zlim([0 0.2]);

function csall = findmode(k_list,zc,p1,alg, de_singularize, do_balancing, eig_spectrum, f)
    addvar = struct('zL1',zc);
    h = 3*2*pi./real(k_list)+p1.H;
    csall = nan(4,length(k_list));
    for i = 1:length(k_list)
        p1.k = k_list(i); p1.h = h(i);
        oall = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        oall = oall(real(oall)>-50); % Remove the eigenvalues assigned by de-singularizing
        [o, osort] = maxeig(oall);
        if real(o) > 0 && imag(o)>1e-5
            addvar.zL1 = p1.criticalH(real(o)/k_list(i));
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