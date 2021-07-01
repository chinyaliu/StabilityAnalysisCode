% close all; clear all;% clc
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig"];
alg = solveGEPmeth(1);
do_balancing = 'n';
eig_spectrum = 'all';
Re = 1000;
Fr2 = 2.25;
N = 600;
k = linspace(0.01,4,200);
% h = 6*ones(1,length(k));
h = 2*pi./k;
inflec_pt = -0.74708299;
zL = 0.74708299;
cutz = NaN(1,length(k));
cutz(1) = -inflec_pt;
numberofDDM = 4;
eps = 0.01;
f = wZhang_ddm.ddmtype(numberofDDM);
in_init = {N,k(1),h(1),Re,Fr2};
%% Select specific modes to observe
% 100
cs(1,1) = -5.69935554321657 + 0.155435654161913i;
cs(1,2) = 0.0943324372553924 - 0.359774643534542i;
cs(1,3) = 7.64805342765973 - 0.0885146263928389i;
% 1000
cs(2,1) = -5.69042132200112 + 0.0157451056762317i;
cs(2,2) = -0.0903635749346843 - 0.0288947409210721i;
cs(2,3) = 7.64522803387940 - 0.00878703878220784i;
% inf
cs(3,1) = -5.69036042186261;
cs(3,2) = 0.0246719832072093 + 0.0344079342992104i;
cs(3,3) = 7.64521877450795;
%% Run solver
tic;
p1 = wZhang_ddm(in_init{:});
p1.numMeth(method);
oall = cell(1,3);
R = [1e2,1e3,inf];
css = cell(3,1);
oss = cell(3,1);
for j = 1:length(R)
    fprintf('Re = %4d\n', R(j));
    p1.chgRe(R(j));
    addvar = struct('zL1',zL,'eps',eps);
    css{j} = NaN(3,length(k));
    oss{j} = NaN(3,length(k));
    for i = 1:length(k)
        fprintf('k = %.2f\n', k(i));
        p1.k = k(i); p1.h = h(i);
        oall = p1.solver(alg, do_balancing, eig_spectrum, f, addvar);
        o = maxeig(oall);
        if real(o) > 0
            cutz(i)=-p1.criticalH(real(o)/k(i));
        else
            cutz(i)=cutz(1);
        end
        addvar.zL1 = cutz(i);
        
        call = oall/k(i);
        if isinf(R(j))
            a = 1:length(call); 
            crange = ((real(call)-0.0012>-1e-5) & (real(call)-1<=1e-5));
            dis = ((real(call)-1).^2 +imag(call).^2)>1e-5;
            aa = a(crange&dis);
            abch = isoutlier(imag(call(aa)),'movmedian',5);
            aa = [a(dis&~crange) aa(abch)];
            cnind = call(aa);
        else
            if (i~=1) && (abs(k(i)-0.8)>1e-5)
                cmatn = repmat(calln,1,length(call));
                cmat = repmat(call.',length(calln),1);
                [~,ind] = min(abs(cmatn-cmat),[],2);
                cnind = call(setdiff(1:end,ind));
                calln = call(ind);
            else
                cnind = call(imag(call)>-2);
                c1 = repmat(cs(j,:).',1,length(cnind));
                c2 = repmat(cnind.',length(cs),1);
                [~,ind] = min(abs(c1-c2),[],2);
                calln = cnind(setdiff(1:end,ind));
            end
        end

        % Choose selected eigenvalues
        if ~isempty(cnind)
            for m = 1:size(cs,1)
                [mc,ind] = min(abs(cnind-cs(j,m)));
                c = cnind(ind);
                [~,ind2] = min(abs(c-cs(j,:)));
                if ind2 == j
                    cs(j,m) = c;
                    css{j}(m,i) = c;
                    oss{j}(m,i) = c*k(i);
                end
            end
        end
    end
end
toc;

%% Read Zhang's results & Plot oi vs k
mk = {'o','+','x'};
ln = {'--','-.','-'};
col = get(gca,'colororder');
dn = {'100','1000','$\infty$'};
imagedata = imread('Zhang.bmp');
im2 = imbinarize(imagedata).*255;
im2(im2==0) = 100;
im2 = cast(im2,'uint8');
fig = figure('position',[50,0,1000,720]);
imagesc([0,4],[0.04,0],im2);
set(gca,'YDir','normal');
xticks(0:1:4);
yticks(0:0.01:0.04);
hold on;
for i = 1:3
    [~,ind] = max(imag(oss{i}),[],1,'linear');
    osmax = oss{i}(ind);
    plot(k,imag(osmax),'.','markersize',6,'color',col(i,:));
end
hold off;
xlabel('$k$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot oi vs k
fig = figure('position',[50,0,1000,720]);
hold on;
for i = 1:3
    for j = 1:3
        plot(k,imag(oss{i}(j,:)),ln{i},'linewidth',1.5,'color',col(j,:));
    end
end
hold off;
xlabel('$k$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
xlim([0 4]);
ylim([-0.03 0.03]);
xticks(0:1:4);
yticks(-0.03:0.01:0.03);

%% Plot c real vs imag
figure; 
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for i = 1:3
    for j = 1:3
        plot(real(css{i}(j,:)),imag(css{i}(j,:)),'.','Markersize',12,'color',col(j,:));
    end
end
hold off; box on;
xlabel('$c_r$','fontsize',30);
xticks(0:1:4);
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot real omega
figure; 
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for i = 1:3
    for j = 1:3
        plot(k,real(oss{i}(j,:)),ln{i},'linewidth',1.5,'color',col(j,:));
    end
end
hold off; box on;
xlabel('$k$','fontsize',30);
xticks(0:1:4);
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('$\omega_r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot real c
figure; 
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for i = 1:3
    for j = 1:3
        plot(k,real(css{i}(j,:)),ln{i},'linewidth',1.5,'color',col(j,:));
    end
end
hold off;
xlabel('$k$','fontsize',30);
xticks(0:1:4);
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('$c_r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

function o = maxeig(ev)
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