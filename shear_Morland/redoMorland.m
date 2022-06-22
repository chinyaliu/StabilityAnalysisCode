clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    

%% Inputs
[method,alg,~,de_singularize,do_balancing,~,N,ud_nd,~,~,~,~,f,epss,Re] = pars_Morland(1,'inv');
eig_spectrum = 'max';
ptnum = 30;

%% Find maximum
baseflowlist = ["exponential", "error function"];
t1 = tic;
parfor j = 1:length(baseflowlist)
    bflow = baseflowlist(j);
    delta_list = linspace(0,0.8,ptnum);
    lambda_list = NaN(1,length(delta_list));
    o_list = NaN(1,length(delta_list));
    for i = 1:length(delta_list)
        flow1 = wMorland(N,1,ud_nd,delta_list(i),1,method,bflow,Re);
        [xL, xU] = findneutral(ud_nd,delta_list(i));
        if ~isnan(xL)
            advar = {flow1,alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',-0.1,'eps',epss)};
            lambda = goldensearchmax(xL,xU,@solveflow,advar{:});
            if ~isnan(lambda)
                lambda_list(i) = lambda;
                flow1.setprop('lambda',lambda,'h',2*lambda);
                o = flow1.solver(advar{2:end});
                o_list(i) = o(1);
            end
        end
    end
    data(j) = struct('delta',delta_list,'o',o_list,'c',0.5*o_list.*lambda_list/pi,'lambda',lambda_list);
end
toc(t1);
global expdata erfdata
expdata = data(1); erfdata = data(2);

%% plot growth rate
ff(1) = plotresult('Morland_growth.bmp',[0,0.8],[1.5,0],(@(x) x.('delta')),(@(x) imag(x.('o'))));
xlabel('$\Delta$');
ylabel('$\omega_i$','rotation',0);
legend('fontsize',18);

%% plot phase speed
ff(2) = plotresult('Morland_cr.bmp',[0,0.8],[0.5,0],(@(x) x.('delta')),(@(x) real(x.('c'))));
xlabel('$\Delta$');
ylabel('$c_r$','rotation',0);
legend('location','southeast','fontsize',18);

%% plot lambda
ff(3) = plotresult('Morland_lam.bmp',[0,0.8],[2.5,0],(@(x) x.('delta')),(@(x) x.('lambda')));
xlabel('$\Delta$');
ylabel('$\lambda$','rotation',0);
legend('location','southeast','fontsize',18);

%%
function f = plotresult(fignam,xL,yL,x,y)
% Read Morland's results
f = figure('position',[150 100 720 640]);
im1 = imread(fignam);
imagesc(xL,yL,im1);
set(gca,'YDir','normal');
% Plot results
global expdata erfdata
expdisp = {'bx','Displayname','exponential','Markersize',8};
erfdisp = {'ro','Displayname','error function','Markersize',8};
hold on;
plot(-1,-1,'-k','Displayname','Morland piecewise-linear');
plot(-1,-1,'--k','Displayname','Morland error function');
plot(x(erfdata),y(erfdata),erfdisp{:},'linewidth',2);
plot(x(expdata),y(expdata),expdisp{:},'linewidth',2);
hold off;
grid on; axis square;
end

function oi = solveflow(lam,flow1,alg, de_singularize, do_balancing, eig_spectrum, f, addvar)
    flow1.setprop('lambda',lam,'h',2*lam);
    o = flow1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    oi = imag(o(1));
end