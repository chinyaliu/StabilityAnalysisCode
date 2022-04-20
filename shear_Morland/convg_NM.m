clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end 
%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,~,~,ud_nd,delta_nd,lambda_nd,c0,h,f,epss] = pars_Morland(2);
N = 100:100:1500;
eig_spectrum = 'max';
lnk = length(lambda_nd);
lnn = length(N);
hh = h(1)/lambda_nd(1);
Re = inf;

%% Convergence test folder
subFolder = ['fig_convergence/' datestr(datetime('now','TimeZone','local'),'mm-dd-HHMM')];
if ~exist(subFolder, 'dir')
    mkdir(subFolder);
else
    error('Folder already exists.');
end

%% Flow property
fid = fopen([subFolder '/_caseinfo.txt'],'w');
fprintf(fid, '%s\n','Initial:');
fprintf(fid,'%15s %10s\n','Base flow',bflow);
fprintf(fid,'%15s %.1e\n','Re',Re);
fprintf(fid,'%15s %15s\n','B.B.C.','Exponential decay');
fprintf(fid,'%15s %1.1f wavelength\n','h',hh);
fprintf(fid,'%15s %10s\n','DDM method',func2str(f));
fprintf(fid,'%15s %4s\n','GEP solver',alg);
fprintf(fid,'%15s %1s\n','Balancing',do_balancing);
fprintf(fid,'%15s %1s\n','De-singularize',de_singularize);
fprintf(fid, '%s\n','Test Case:');
fprintf(fid,'%7s %7s %7s\n','ud','delta','lambda');
for i = 1:lnk
    fprintf(fid,'%7.3f %7.3f %7.3f\n',ud_nd(i),delta_nd(i),lambda_nd(i));
end
fprintf(fid,'%s\n',repelem('=',40));
fclose(fid);

%% Infinite boundary & truncation height
nh = linspace(0.5,5,10);
lnh = length(nh);
cih1 = nan(lnk,lnh);
Nh1 = nan(lnk,lnh);
cih2 = cih1; Nh2 = Nh1;
tic;
parfor i = 1:lnk
    nhloop = nh;
    for j = 1:lnh
        htemp = nhloop(j)*lambda_nd(i);
        case1 = wMorland(100,htemp,ud_nd(i),delta_nd(i),lambda_nd(i),method,bflow);
        addvar = struct('zL1',case1.invbf(c0(i)),'eps',epss);
        [cih1(i,j), Nh1(i,j)] = convgmode(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        [cih2(i,j), Nh2(i,j)] = convgmode2(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    end
end
toc;
% Analyze convergence
dd1 = []; dd2 = [];
for i = 1:lnk
    dtemp = diff(abs(cih1(i,:)-cih1(i,end)));
    dd1(i) = find(dtemp<5*dtemp(end-1),1,'first');
    dtemp = diff(abs(cih2(i,:)-cih2(i,end)));
    dd2(i) = find(dtemp<5*dtemp(end-1),1,'first');
end
ind1 = max(dd1); ind2 = max(dd2);
bbc = 1; % exponential decay
nh0 = ceil(nh(ind1+1))+1;
if ind1>ind2
    bbc = 2; % free-slip
    nh0 = ceil(nh(ind2+1))+1;
end    
h = h/hh*nh0;
%% Plot
for i = 1:lnk
    nam(i) = "$\"+sprintf("lambda = %.3f$",lambda_nd(i));
end
f1(1) = plotcon(nh,abs(cih1-cih1(:,end)),'$h/\lambda$','$\ | \ c_i - c_{0,i}\ |$',nam);
xlim([0.5 5]);
xticks(1:1:5);
ylim([1e-16 1e-2]);
expfig([subFolder '/exp']);
f1(2) = plotcon(nh,abs(cih2-cih2(:,end)),'$h/\lambda$','$\ | \ c_i - c_{0,i}\ |$',nam);
xlim([0.5 5]);
xticks(1:1:5);
ylim([1e-16 1e-2]);
expfig([subFolder '/fs']);

%% Domain decomposition
ci1 = nan(lnk,lnn);
ci2 = ci1; ci4 = ci1;
tic;
parfor i = 1:lnk
    nloop = N;
    for j = 1:lnn
        case1 = wMorland(100,h(i),ud_nd(i),delta_nd(i),lambda_nd(i),method,bflow);
        addvar = struct('zL1',case1.invbf(c0(i)),'eps',epss);
        case1.N = nloop(j);
        c = solver(bbc,case1,alg,de_singularize, do_balancing, eig_spectrum, wMorland.ddmtype(1), addvar);
        ci1(i,j) = imag(c(1));
        c2 = solver(bbc,case1,alg,de_singularize, do_balancing, eig_spectrum, wMorland.ddmtype(2), addvar);
        ci2(i,j) = imag(c2(1));
        c4 = solver(bbc,case1,alg,de_singularize, do_balancing, eig_spectrum, wMorland.ddmtype(44), addvar);
        ci4(i,j) = imag(c4(1));
    end
end
toc;
% Analyze convergence
dd1 = []; dd2 = []; dd4 = [];
for i = 1:lnk
    dtemp = abs(ci1(i,:)-ci1(i,end));
    dd1(i) = sum(dtemp(end-4:end));
    dtemp = abs(ci2(i,:)-ci2(i,end));
    dd2(i) = sum(dtemp(end-4:end));
    dtemp = abs(ci4(i,:)-ci4(i,end));
    dd4(i) = sum(dtemp(end-4:end));
end
ind = [max(dd1) max(dd2) max(dd4)];
[~,ind] = min(ind);
switch ind
case 1
    f = wMorland.ddmtype(1);
case 2
    f = wMorland.ddmtype(2);
case 3
    f = wMorland.ddmtype(44);
end
%% Plot
f1(3) = plotcon(N,abs(ci1-ci1(:,end)),'$N$','$\ | \ c_i - c_{0,i}\ |$',nam);
xlim([N(1) N(end)]);
ylim([1e-17 1e-1]);
expfig([subFolder '/d1']);
f1(4) = plotcon(N,abs(ci2-ci2(:,end)),'$N$','$\ | \ c_i - c_{0,i}\ |$',nam);
xlim([N(1) N(end)]);
ylim([1e-17 1e-1]);
expfig([subFolder '/d2']);
f1(5) = plotcon(N,abs(ci4-ci4(:,end)),'$N$','$\ | \ c_i - c_{0,i}\ |$',nam);
xlim([N(1) N(end)]);
ylim([1e-17 1e-1]);
expfig([subFolder '/d4']);

%% GEP solving algorithm
alglist = ["eig", "qr", "invB"];
tt1 = nan(lnk,length(nh));
tt2 = tt1; tt3 = tt1;
ci1 = tt1; ci2 = tt1; ci3 = tt1;
tic;
parfor i = 1:lnk
    algloop = ["eig", "qr", "invB"];
    nloop = N;
    for j = 1:lnn
        case1 = wMorland(100,h(i),ud_nd(i),delta_nd(i),lambda_nd(i),method,bflow);
        addvar = struct('zL1',case1.invbf(c0(i)),'eps',epss);
        case1.N = nloop(j);
        t1 = tic;
        c = solver(bbc,case1,algloop(1), de_singularize, do_balancing, eig_spectrum, f, addvar);
        tt1(i,j) = toc(t1);
        ci1(i,j) = imag(c(1));
        t1 = tic;
        c = solver(bbc,case1,algloop(2), de_singularize, do_balancing, eig_spectrum, f, addvar);
        tt2(i,j) = toc(t1);
        ci2(i,j) = imag(c(1));
        t1 = tic;
        c = solver(bbc,case1,algloop(3), de_singularize, do_balancing, eig_spectrum, f, addvar);
        tt3(i,j) = toc(t1);
        ci3(i,j) = imag(c(1));
    end
end
toc;
% Analyze convergence
dd1 = []; dd2 = []; dd3 = [];
for i = 1:lnk
    dtemp = abs(ci1(i,:)-ci1(i,end));
    dd1(i) = sum(dtemp(end-4:end));
    dtemp = abs(ci2(i,:)-ci2(i,end));
    dd2(i) = sum(dtemp(end-4:end));
    dtemp = abs(ci3(i,:)-ci3(i,end));
    dd3(i) = sum(dtemp(end-4:end));
end
ind = [max(dd1) max(dd2) max(dd3)];
[~,ind] = min(ind);
alg = alglist(ind);
%% Plot
nam = ["$\mathrm{qz}$" "$\mathrm{inv}(A)$" "$\mathrm{inv}(B)$"];
ct = length(f1);
for i = 1:lnk
    y1 = [ci1(i,:); ci2(i,:); ci3(i,:)];
    y2 = [tt1(i,:); tt2(i,:); tt3(i,:)];
    f1(ct+1) = plotcon(N,abs(y1-y1(:,end)),'$N$','$\ | \ c_i - c_{0,i}\ |$',nam);
    xlim([N(1) N(end)]);
    ylim([1e-17 1e-3]);
    expfig([subFolder '/' sprintf('lam%03.0f',100*lambda_nd(i))]);
    f1(ct+2) = plotcon(N,y2,'$N$','computational time (s)',nam);
    legend('location','northwest');
    xlim([N(1) N(end)]);
    set(gca,'YScale','linear');
    expfig([subFolder '/' sprintf('lam%03.0f_time',100*lambda_nd(i))]);
    ct = ct+2;
end

%% De-singularizing
if strcmp(alg,'invB')
    de_singularize = 'y';
else
    ci1 = nan(lnk,length(nh));
    ci2 = ci1;
    tic;
    parfor i = 1:lnk
        nloop = N;
        for j = 1:lnn
            case1 = wMorland(100,h(i),ud_nd(i),delta_nd(i),lambda_nd(i),method,bflow);
            addvar = struct('zL1',case1.invbf(c0(i)),'eps',epss);
            case1.N = nloop(j);
            c = solver(bbc,case1,alg, 'n', do_balancing, eig_spectrum, f, addvar);
            ci1(i,j) = imag(c(1));
            c = solver(bbc,case1,alg, 'y', do_balancing, eig_spectrum, f, addvar);
            ci2(i,j) = imag(c(1));
        end
    end
    toc;
    % Analyze convergence
    dd1 = []; dd2 = [];
    for i = 1:lnk
        dtemp = abs(ci1(i,:)-ci1(i,end));
        dd1(i) = sum(dtemp(end-4:end));
        dtemp = abs(ci2(i,:)-ci2(i,end));
        dd2(i) = sum(dtemp(end-4:end));
    end
    ind = [max(dd1) max(dd2)];
    [~,ind] = min(ind);
    if ind == 1
        de_singularize = 'n';
    else
        de_singularize = 'y';
    end

    % Plot
    nam = ["original" "de-singularized"];
    ct = length(f1);
    for i = 1:lnk
        ct = ct+1;
        y1 = [ci1(i,:); ci2(i,:);];
        y2 = [tt1(i,:); tt2(i,:);];
        f1(ct) = plotcon(N,abs(y1-y1(:,end)),'$N$','$\ | \ c_i - c_{0,i}\ |$',nam);
        xlim([N(1) N(end)]);
        ylim([1e-17 1e-3]);
        expfig([subFolder '/' sprintf('lam%03.0f_sing',100*lambda_nd(i))]);
    end
end

%% Balancing
ci1 = nan(lnk,length(nh));
ci2 = ci1;
tic;
parfor i = 1:lnk
    nloop = N;
    for j = 1:lnn
        case1 = wMorland(100,h(i),ud_nd(i),delta_nd(i),lambda_nd(i),method,bflow);
        addvar = struct('zL1',case1.invbf(c0(i)),'eps',epss);
        case1.N = nloop(j);
        c = solver(bbc,case1,alg, de_singularize, 'n', eig_spectrum, f, addvar);
        ci1(i,j) = imag(c(1));
        c = solver(bbc,case1,alg, de_singularize, 'y', eig_spectrum, f, addvar);
        ci2(i,j) = imag(c(1));
    end
end
toc;
% Analyze convergence
dd1 = []; dd2 = [];
for i = 1:lnk
    dtemp = abs(ci1(i,:)-ci1(i,end));
    dd1(i) = sum(dtemp(end-4:end));
    dtemp = abs(ci2(i,:)-ci2(i,end));
    dd2(i) = sum(dtemp(end-4:end));
end
ind = [max(dd1) max(dd2)];
[~,ind] = min(ind);
if ind == 1
    do_balancing = 'n';
else
    do_balancing = 'y';
end

% Plot
nam = ["original" "balanced"];
ct = length(f1);
for i = 1:lnk
    ct = ct+1;
    y1 = [ci1(i,:); ci2(i,:);];
    y2 = [tt1(i,:); tt2(i,:);];
    f1(ct) = plotcon(N,abs(y1-y1(:,end)),'$N$','$\ | \ c_i - c_{0,i}\ |$',nam);
    xlim([N(1) N(end)]);
    ylim([1e-17 1e-3]);
    expfig([subFolder '/' sprintf('lam%03.0f_bal',100*lambda_nd(i))]);
end

%% Flow property (Updated)
if bbc == 1
    Bbc = 'Exponential decay';
else
    Bbc = 'Free-slip';
end
fid = fopen([subFolder '/_caseinfo.txt'],'a');
fprintf(fid, '\n%s\n','New:');
fprintf(fid,'%15s %10s\n','Base flow',bflow);
fprintf(fid,'%15s %.1e\n','Re',Re);
fprintf(fid,'%15s %15s\n','B.B.C.','Exponential decay');
fprintf(fid,'%15s %1.1f wavelength\n','h',nh0);
fprintf(fid,'%15s %10s\n','DDM method',func2str(f));
fprintf(fid,'%15s %4s\n','GEP solver',alg);
fprintf(fid,'%15s %1s\n','Balancing',do_balancing);
fprintf(fid,'%15s %1s\n','De-singularize',de_singularize);
fprintf(fid,'%s\n',repelem('=',40));
fclose(fid);

function f1 = plotcon(x,y,xl,yl,nam)
lsy = {'-ko',':ro','--bo','-.g'};
f1 = figure; hold on;
for i = 1:length(nam)
    f(i) = plot(x, y(i,:), lsy{i}, 'Displayname', nam(i));
end
hold off; grid on; box on;
set(f(:),'linewidth',3,'markersize',8);
set(gca,'fontsize',30, 'YScale', 'log');
xticks(500:500:1500);
xlabel(xl,'fontsize',36);
ylabel(yl,'fontsize',36);
legend('location','northeast');
end

function [ci, Nc] = convgmode(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar)
    Nc = N(end);
    ctemp = 0;
    for i = 1:length(N)
        case1.N = N(i);
        c = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        ci = imag(c);
        if ~isnan(case1.zc)
            addvar.zL1 = case1.zc;
            if i~=1
                if (abs((ci-ctemp)/ci)<1e-8 || i == length(N))
                    Nc = N(i);
                    break;
                else
                    ctemp = ci;
                end
            else
                ctemp = imag(c);
            end
        end
    end
end

function [ci, Nc] = convgmode2(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar)
    Nc = N(end);
    ctemp = 0;
    for i = 1:length(N)
        case1.N = N(i);
        c = case1.solver2(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        ci = imag(c);
        if ~isnan(case1.zc)
            addvar.zL1 = case1.zc;
            if i~=1
                if (abs((ci-ctemp)/ci)<1e-8 || i == length(N))
                    Nc = N(i);
                    break;
                else
                    ctemp = ci;
                end
            else
                ctemp = imag(c);
            end
        end
    end
end

function o = solver(num, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar)
    if num == 1
        o = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    else
        o = case1.solver2(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    end
end