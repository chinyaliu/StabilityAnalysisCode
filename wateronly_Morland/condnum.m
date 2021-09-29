clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,ud_nd,delta_nd,lambda_nd,c0,h,f] = pars_Morland(2);
testpar = 'h';

%% Find the condition numbers
t1 = tic;
case1 = wMorland(N,h,ud_nd,delta_nd,lambda_nd,method,bflow);
addvar = struct('zL1',-case1.criticalH(c0),'eps',0.1);
switch testpar
    case('h')
        test_list = linspace(0.5,5,30);
        n_list = N*ones(1,length(test_list));
        h_list = test_list.*lambda_nd;
    case('n')
        test_list = 300:50:1500;
        n_list = test_list;
        h_list = h*ones(1,length(test_list));
end
cA = NaN(1,length(test_list));
cAb = NaN(1,length(test_list));
cB = NaN(1,length(test_list));
cBb = NaN(1,length(test_list));
for ii = 1:length(test_list)
    fprintf('N = %4d, h = %.4f\n', n_list(ii), h_list(ii));
    case1.N = n_list(ii);
    case1.h = h_list(ii);
    [cA(ii),cB(ii),cAb(ii),cBb(ii)] = case1.calcond(de_singularize, do_balancing, f, addvar);
end
toc(t1);
%% Plot results
figure;
hold on;
plot(test_list, cA, '--bo', 'Displayname', 'A');
plot(test_list, cAb, '-bo', 'Displayname', 'balanced A');
plot(test_list, cB, '--ro', 'Displayname', 'B');
plot(test_list, cBb, '-ro', 'Displayname', 'balanced B');
hold off; box on; grid on;
set(gca,'YScale','log');
legend('location', 'southeast');
ylabel('Condition number of inverse');
switch testpar
    case('h')
        xlabel('$n\lambda $');
    case('n')
        xlabel('$N$');
end