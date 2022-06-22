clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end
%% Solver & Algorithm list
[method,~,bflow,de_singularize,do_balancing,~,~,ud_nd,delta_nd,lambda_nd,c0,h,f,epss,Re] = pars_Morland;
N = 100:100:1500;
eig_spectrum = 'max';
alg = ["qr", "qz", "invB"];
lnk = length(alg);

%% Run solver
nh = linspace(0.5,5,20);
lnh = length(nh);
oi = nan(lnk,lnh);
Ni = oi;
tic;
parfor i = 1:lnk
    nhloop = nh;
    for j = 1:lnh
        htemp = nhloop(j)*lambda_nd;
        case1 = wMorland(100,htemp,ud_nd,delta_nd,lambda_nd,method,bflow,Re);
        addvar = struct('zL1',case1.invbf(c0),'eps',epss);
        [oi(i,j), Ni(i,j)] = convgmode(N, case1, alg(i), de_singularize, do_balancing, eig_spectrum, f, addvar);
    end
end
doi = abs(diff(oi,1,2));
toc;

%% Plot figure
lsy = {'-ko',':ro','--bo','-.go'};
figure;
nam = ["inv(A)","qz","inv(B)"];
yf(1) = semilogy(nh(1:end-1), doi(1,:), lsy{1}, 'Displayname',nam(1));
hold on;
for i = 2:lnk
    yf(i) = semilogy(nh(1:end-1), doi(i,:), lsy{i}, 'Displayname', nam(i));
end
hold off;
set(yf(:),'linewidth',3,'markersize',8);
set(gca,'fontsize',30);
xlim([0.5 5])
xlabel('$h/\lambda$','fontsize',36);
ylabel('$\ | \ \omega_{i,m} - \omega_{i,m+1}\ |$','fontsize',36);
legend('location','east');
grid on;

%%
function [oi, Nc] = convgmode(N, case1, alg, de_singularize, do_balancing, eig_spectrum, f, addvar)
    Nc = N(end);
    otemp = 0;
    for i = 1:length(N)
        case1.N = N(i);
        o = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        oi = imag(o);
        if ~isnan(case1.zc)
            addvar.zL1 = case1.zc;
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
end