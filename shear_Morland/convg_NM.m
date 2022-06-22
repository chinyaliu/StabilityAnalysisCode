clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end 
%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,~,~,ud_nd,delta_nd,lambda_nd,c0,h,f,epss,Re] = pars_Morland(6);
N = 100:100:1500;
eig_spectrum = 'max';
lnk = length(lambda_nd);
lnn = length(N);
hh = h(1)/lambda_nd(1);
if length(ud_nd)==1
    ud_nd = ud_nd*ones(1,lnk);
    delta_nd = delta_nd*ones(1,lnk);
end
if length(Re) == 1
    Re = Re*ones(1,lnk);
end
pltylab = '$\ | \ \omega_i - \omega_{0,i}\ |$';

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
fprintf(fid,'%15s %15s\n','B.B.C.','Exponential decay');
fprintf(fid,'%15s %1.1f wavelength\n','h',hh);
fprintf(fid,'%15s %10s\n','DDM method',func2str(f));
fprintf(fid,'%15s %4s\n','GEP solver',alg);
fprintf(fid,'%15s %1s\n','Balancing',do_balancing);
fprintf(fid,'%15s %1s\n','De-singularize',de_singularize);
fprintf(fid, '%s\n','Test Case:');
fprintf(fid,'%7s %7s %7s %7s\n','ud','delta','lambda','Re');
for i = 1:lnk
    fprintf(fid,'%7.3f %7.3f %7.3f %.1e\n',ud_nd(i),delta_nd(i),lambda_nd(i),Re(i));
end
fprintf(fid,'%s\n',repelem('=',40));
fclose(fid);

%% Infinite boundary & truncation height
nh = linspace(0.5,5,10);
lnh = length(nh);
oih1 = nan(lnk,lnh);
Nh1 = nan(lnk,lnh);
oih2 = oih1; Nh2 = Nh1;
tic;
parfor i = 1:lnk
    nhloop = nh;
    for j = 1:lnh
        htemp = nhloop(j)*lambda_nd(i);
        case1 = wMorland(100,htemp,ud_nd(i),delta_nd(i),lambda_nd(i),method,bflow,Re(i));
        addvar = struct('zL1',case1.invbf(c0(i)),'eps',epss);
        [oih1(i,j), Nh1(i,j)] = convgmode(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        [oih2(i,j), Nh2(i,j)] = convgmode2(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    end
end
toc;
% Analyze convergence
dd1 = []; dd2 = [];
for i = 1:lnk
    dtemp = abs(diff(abs(oih1(i,:)-oih1(i,end))));
    dd1(i) = find(dtemp<5*dtemp(end-1),1,'first');
    dtemp = diff(abs(oih2(i,:)-oih2(i,end)));
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
f1(1) = plotcon(nh,abs(oih1-oih1(:,end)),'$h/\lambda$',pltylab,nam);
xlim([0.5 5]);
xticks(1:1:5);
ylim([1e-16 1e-2]);
expfig([subFolder '/exp']);
f1(2) = plotcon(nh,abs(oih2-oih2(:,end)),'$h/\lambda$',pltylab,nam);
xlim([0.5 5]);
xticks(1:1:5);
ylim([1e-16 1e-2]);
expfig([subFolder '/fs']);

%% Domain decomposition
oi1 = nan(lnk,lnn);
oi2 = oi1; oi4 = oi1;
tic;
parfor i = 1:lnk
    nloop = N;
    for j = 1:lnn
        case1 = wMorland(100,h(i),ud_nd(i),delta_nd(i),lambda_nd(i),method,bflow,Re(i));
        addvar = struct('zL1',case1.invbf(c0(i)),'eps',epss);
        case1.N = nloop(j);
        o = solver(bbc,case1,alg,de_singularize, do_balancing, eig_spectrum, wMorland.ddmtype(1), addvar);
        oi1(i,j) = imag(o(1));
        o2 = solver(bbc,case1,alg,de_singularize, do_balancing, eig_spectrum, wMorland.ddmtype(2), addvar);
        oi2(i,j) = imag(o2(1));
        o4 = solver(bbc,case1,alg,de_singularize, do_balancing, eig_spectrum, wMorland.ddmtype(44), addvar);
        oi4(i,j) = imag(o4(1));
    end
end
toc;
% Analyze convergence
dd1 = []; dd2 = []; dd4 = [];
for i = 1:lnk
    dtemp = abs(oi1(i,:)-oi1(i,end));
    dd1(i) = sum(dtemp(end-4:end));
    dtemp = abs(oi2(i,:)-oi2(i,end));
    dd2(i) = sum(dtemp(end-4:end));
    dtemp = abs(oi4(i,:)-oi4(i,end));
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
f1(3) = plotcon(N,abs(oi1-oi1(:,end)),'$N$',pltylab,nam);
xlim([N(1) N(end)]);
ylim([1e-17 1e-1]);
expfig([subFolder '/d1']);
f1(4) = plotcon(N,abs(oi2-oi2(:,end)),'$N$',pltylab,nam);
xlim([N(1) N(end)]);
ylim([1e-17 1e-1]);
expfig([subFolder '/d2']);
f1(5) = plotcon(N,abs(oi4-oi4(:,end)),'$N$',pltylab,nam);
xlim([N(1) N(end)]);
ylim([1e-17 1e-1]);
expfig([subFolder '/d4']);

%% GEP solving algorithm
alglist = ["eig", "qr", "invB"];
tt1 = nan(lnk,length(nh));
tt2 = tt1; tt3 = tt1;
oi1 = tt1; oi2 = tt1; oi3 = tt1;
tic;
parfor i = 1:lnk
    algloop = ["eig", "qr", "invB"];
    nloop = N;
    for j = 1:lnn
        case1 = wMorland(100,h(i),ud_nd(i),delta_nd(i),lambda_nd(i),method,bflow,Re(i));
        addvar = struct('zL1',case1.invbf(c0(i)),'eps',epss);
        case1.N = nloop(j);
        t1 = tic;
        o = solver(bbc,case1,algloop(1), de_singularize, do_balancing, eig_spectrum, f, addvar);
        tt1(i,j) = toc(t1);
        oi1(i,j) = imag(o(1));
        t1 = tic;
        o = solver(bbc,case1,algloop(2), de_singularize, do_balancing, eig_spectrum, f, addvar);
        tt2(i,j) = toc(t1);
        oi2(i,j) = imag(o(1));
        t1 = tic;
        o = solver(bbc,case1,algloop(3), de_singularize, do_balancing, eig_spectrum, f, addvar);
        tt3(i,j) = toc(t1);
        oi3(i,j) = imag(o(1));
    end
end
toc;
% Analyze convergence
dd1 = []; dd2 = []; dd3 = [];
for i = 1:lnk
    dtemp = abs(oi1(i,:)-oi1(i,end));
    dd1(i) = sum(dtemp(end-4:end));
    dtemp = abs(oi2(i,:)-oi2(i,end));
    dd2(i) = sum(dtemp(end-4:end));
    dtemp = abs(oi3(i,:)-oi3(i,end));
    dd3(i) = sum(dtemp(end-4:end));
end
ind = [max(dd1) max(dd2) max(dd3)];
[~,ind] = min(ind);
alg = alglist(ind);
%% Plot
nam = ["$\mathrm{qz}$" "$\mathrm{inv}(A)$" "$\mathrm{inv}(B)$"];
ct = length(f1);
for i = 1:lnk
    y1 = [oi1(i,:); oi2(i,:); oi3(i,:)];
    y2 = [tt1(i,:); tt2(i,:); tt3(i,:)];
    f1(ct+1) = plotcon(N,abs(y1-y1(:,end)),'$N$',pltylab,nam);
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
    oi1 = nan(lnk,length(nh));
    oi2 = oi1;
    tic;
    parfor i = 1:lnk
        nloop = N;
        for j = 1:lnn
            case1 = wMorland(100,h(i),ud_nd(i),delta_nd(i),lambda_nd(i),method,bflow,Re(i));
            addvar = struct('zL1',case1.invbf(c0(i)),'eps',epss);
            case1.N = nloop(j);
            o = solver(bbc,case1,alg, 'n', do_balancing, eig_spectrum, f, addvar);
            oi1(i,j) = imag(o(1));
            o = solver(bbc,case1,alg, 'y', do_balancing, eig_spectrum, f, addvar);
            oi2(i,j) = imag(o(1));
        end
    end
    toc;
    % Analyze convergence
    dd1 = []; dd2 = [];
    for i = 1:lnk
        dtemp = abs(oi1(i,:)-oi1(i,end));
        dd1(i) = sum(dtemp(end-4:end));
        dtemp = abs(oi2(i,:)-oi2(i,end));
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
        y1 = [oi1(i,:); oi2(i,:);];
        y2 = [tt1(i,:); tt2(i,:);];
        f1(ct) = plotcon(N,abs(y1-y1(:,end)),'$N$',pltylab,nam);
        xlim([N(1) N(end)]);
        ylim([1e-17 1e-3]);
        expfig([subFolder '/' sprintf('lam%03.0f_sing',100*lambda_nd(i))]);
    end
end

%% Balancing
oi1 = nan(lnk,length(nh));
oi2 = oi1;
tic;
parfor i = 1:lnk
    nloop = N;
    for j = 1:lnn
        case1 = wMorland(100,h(i),ud_nd(i),delta_nd(i),lambda_nd(i),method,bflow,Re(i));
        addvar = struct('zL1',case1.invbf(c0(i)),'eps',epss);
        case1.N = nloop(j);
        o = solver(bbc,case1,alg, de_singularize, 'n', eig_spectrum, f, addvar);
        oi1(i,j) = imag(o(1));
        o = solver(bbc,case1,alg, de_singularize, 'y', eig_spectrum, f, addvar);
        oi2(i,j) = imag(o(1));
    end
end
toc;
% Analyze convergence
dd1 = []; dd2 = [];
for i = 1:lnk
    dtemp = abs(oi1(i,:)-oi1(i,end));
    dd1(i) = sum(dtemp(end-4:end));
    dtemp = abs(oi2(i,:)-oi2(i,end));
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
    y1 = [oi1(i,:); oi2(i,:);];
    y2 = [tt1(i,:); tt2(i,:);];
    f1(ct) = plotcon(N,abs(y1-y1(:,end)),'$N$',pltylab,nam);
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
fprintf(fid,'%15s %15s\n','B.B.C.',Bbc);
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

function [oi, Nc] = convgmode(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar)
    Nc = N(end);
    otemp = 0;
    for i = 1:length(N)
        case1.N = N(i);
        o = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        oi = imag(o);
        if ~isnan(case1.zc)
            addvar.zL1 = case1.zc;
        end
        if i~=1
            if (abs(oi-otemp)<1e-8 || i == length(N))
                Nc = N(i);
                break;
            else
                otemp = oi;
            end
        else
            otemp = oi;
        end
    end
end

function [oi, Nc] = convgmode2(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar)
    Nc = N(end);
    otemp = 0;
    for i = 1:length(N)
        case1.N = N(i);
        o = case1.solver2(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        oi = imag(o);
        if ~isnan(case1.zc)
            addvar.zL1 = case1.zc;
        end
        if i~=1
            if (abs(oi-otemp)<1e-8 || i == length(N))
                Nc = N(i);
                break;
            else
                otemp = oi;
            end
        else
            otemp = oi;
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