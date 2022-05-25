close all; clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,ud_nd,delta_nd,lamtemp,c0,htemp,f,epss,Re] = pars_Morland;
lmaxmin = findneutral(ud_nd,delta_nd);
lambda_list = linspace(lmaxmin{1}(2),lmaxmin{1}{1},50);
h = lambda_nd*htemp/lamtemp;
c0 = min(sqrt(0.5*lambda_list), 0.99*ud_nd);
c_list = cell(length(lambda_list),1);
in_init = {N,h(1),ud_nd,delta_nd,lambda_list(1),method,bflow,Re};

%% Select mode to observe
col = repmat(get(gca,'colororder'),5,1);
cs(1,1) = 4.1457;
cs(2,1) = -0.35307;
cs_col(1:2,:) = col(1:2,:);
count = 3;

%% Run solver
figure;
yline(0,'linewidth',1.5,'color','#898989');
hold on;
xline(0,'linewidth',1.5,'color','#898989');
xlabel('$\omega_r$');
ylabel('$\omega_i$','rotation',0, 'HorizontalAlignment','right');
%ylim([-0.3 0.15]);
%xlim([-5 20]);
t1 = tic;
p1 = wMorland(in_init{:});
zL = p1.invbf(c0);
cutz = zL(1);
addvar = struct('zL1',zL(1),'eps',epss);
for i = 1:length(lambda_list)
    p1.setprop('lambda',lambda_list(i),'h',h(i));
    fprintf('wavelength = %.2f\n', lambda_list(i));
    oall = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    if ~isnan(p1.zc)
        addvar.zL1 = p1.zc;
    end
    
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

