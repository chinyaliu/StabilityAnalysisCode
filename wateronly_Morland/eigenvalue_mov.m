close all; clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig"];
alg = solveGEPmeth(1);
baseflowlist = ["exponential", "error function"];
bflow = baseflowlist(2);
% Inputs
do_balancing = 'y';
eig_spectrum = 'all';
N = 600;
ud_nd = 2;
delta_nd = 0.291;
% h = @(x) 3*x;
h = 3*delta_nd;
ddm_number = 2;
f = wMorland.ddmtype(ddm_number);

%% Select mode to observe
col = repmat(get(gca,'colororder'),5,1);
cs(1,1) = 4.1457;
cs(2,1) = -0.35307;
cs_col(1:2,:) = col(1:2,:);
count = 3;

%% Run solver
lambda_list = linspace(0.1,4,400);
ki = 0.1i;
c_list = cell(length(lambda_list),1);
% c_list = NaN(2,length(lambda_list));
addvar = struct('zL1',delta_nd,'eps',0.2);
yline(0,'linewidth',1.5,'color','#898989');
hold on;
xline(0,'linewidth',1.5,'color','#898989');
xlabel('$\omega_r$');
ylabel('$\omega_i$','rotation',0, 'HorizontalAlignment','right');
ylim([-0.3 0.15]);
xlim([-5 20]);
t1 = tic;
flow1 = wMorland(N,h,ud_nd,delta_nd,1,method,bflow);
for i = 1:length(lambda_list)
    lam = lambda_list(i);
    fprintf('wavelength = %.2f\n', lam);
    flow1.k = 2*pi/lam-ki;
%     flow1.h = h;
    c_all = flow1.solver(alg, do_balancing, eig_spectrum, f, addvar);
    if isnan(flow1.zc)
        addvar.zL1 = delta_nd;
    else
        addvar.zL1 = -flow1.zc;
    end
    
%     % Choose selected eigenvalues
%     if (i~=1)
%         cmatn = repmat(c_alln,1,length(c_all));
%         cmat = repmat(c_all.',length(c_alln),1);
%         [~,ind] = min(abs(cmatn-cmat),[],2);
%         cnind = c_all(setdiff(1:end,ind));
%         c_alln = c_all(ind);
%     else
%         cnind = c_all;
%         c1 = repmat(cs.',1,length(cnind));
%         c2 = repmat(cnind.',length(cs),1);
%         [~,ind] = min(abs(c1-c2),[],2);
%         c_alln = cnind(setdiff(1:end,ind));
%     end
%     if ~isempty(cnind)
%         % mode 1
%         [mc,ind] = min(abs(cnind-cs(1)));
%         c = cnind(ind);
%         cs(1) = c;
%         c_list(1,i) = c;
%         % mode 2
%         c_max = maxeig(cnind);
%         if ~isnan(c_max)
%             cs(2) = c_max;
%             c_list(2,i) = c_max;
%         elseif sum(real(cnind)<0)
%             [~, ind] = min(real(cnind));
%             cs(2) = cnind(ind);
%             c_list(2,i) = cnind(ind);
%         end
%     end

    % Choose selected eigenvalues
    c_new = finddiscrete(c_all,ud_nd);
    c0 = repmat(cs.',length(c_new),1);
    c1 = repmat(c_new,1,length(cs));
    [~,min_ind] = min(abs(c0-c1),[],1);
    [~,min2] = min(abs(c0-c1),[],2);
%     if length(min_ind) == length(unique(min_ind))
%         cs = c_new(min_ind);
%         c_list{i} = cs;
%     else
%     end
    for j = 1:length(cs)
        if min2(min_ind(j)) == j
            cs(j) = c_new(min_ind(j));
        else
            cs(j) = NaN;
            cs_col(j,:) = NaN;
        end
    end
    cs_col(isnan(cs),:) = [];
    cs(isnan(cs)) = [];
    if length(cs)<length(c_new)
        num = length(c_new)-length(cs);
        cs = [cs;c_new(setdiff(1:end,min_ind))];
        cs_col = [cs_col;col(count:count+num-1,:)];
        count = count+num;
    end
    c_list{i} = cs;
    
    o = c_all*flow1.k;
    g(1) = scatter(real(o),imag(o),'ok');
    for j = 1:length(cs)
        plot(real(cs(j)*flow1.k),imag(cs(j)*flow1.k),'.','Markersize',12,'color',cs_col(j,:));
    end
%     g(1) = scatter(real(c_all),imag(c_all),'ok');
%     plot(real(cs(1)),imag(cs(1)),'.','Markersize',12,'color',col(1,:));
%     plot(real(cs(2)),imag(cs(2)),'.','Markersize',12,'color',col(2,:));
    titext = sprintf('$\\lambda =%.2f$',lam);
    title(titext);
    F(i) = getframe(gcf);
    delete(g);
end
hold off;
close(gcf);
toc(t1);

%%
writerObj = VideoWriter('dt0291o2_mode01.mp4','MPEG-4');
writerObj.FrameRate = 20;
writerObj.Quality = 100;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);

function c = maxeig(ev)
    a = 1:length(ev);
    ev_ind = a(abs(ev)<1e+3 & abs(ev)>1e-6 & abs(imag(ev)) > 1e-5);
    [~, ind] = sort(imag(ev(ev_ind)),'descend');
    ev_ind = ev_ind(ind);
    w = ev(ev_ind);
    if isempty(w)
        c = NaN;
    else
        c = w(1);
    end
end

function c = finddiscrete(ev,ulim)
    a = 1:length(ev); 
    crange = (real(ev)>0) & (real(ev)<ulim);
    aa = a(crange);
    abch = isoutlier(imag(ev(aa)),'movmedian',5);
    newa = aa(abch);
    newa(abs(imag(ev(newa)))<1e-4) = [];
    aa = [a(real(ev)<0) newa a(real(ev)>ulim)];
    c = ev(aa);
end
