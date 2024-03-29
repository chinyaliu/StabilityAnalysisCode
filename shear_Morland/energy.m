clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,ud_nd,delta_nd,lambda_nd,c0,h,f,epss,Re] = pars_Morland(3,'inv');

%% Run
t1 = tic;
for i = 1:length(lambda_nd)
    % Solve for the eigenvalues
    case1 = wMorland(N, h(i), ud_nd, delta_nd, lambda_nd(i), method, bflow, Re);
    addvar = struct('zL1',case1.invbf(c0(i)),'eps',epss);
    [o, an] = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    
    % Choose discrete eigenvalues
    [o_chosen, an_c] = choose_eiv(o, an, case1.k) ;
%     [~,ind] = max(imag(o));
%     o_chosen = o(ind);
%     an_c = an(:,ind);

    % Calculate energy
    chm = 1;
    if imag(o_chosen(chm))>0% && real(o(chm))>0
        if i == 1
            [ENG(i),f1] = plotenergy(case1, "res", o_chosen(chm), an_c(:,chm));
            [~,f2] = plotenergy(case1, "resz", o_chosen(chm), an_c(:,chm));
            axf2 = get(f2,'CurrentAxes');
            title(axf2,'$\frac{d\tau(z)}{dz}$');
            [~,f3] = plotenergy(case1, "invu", o_chosen(chm), an_c(:,chm));
            axf3 = get(f3,'CurrentAxes');
            title(axf3,'$\frac{1}{|U-c|^2}$');
        else
            [ENG(i),f1] = plotenergy(case1, "res", o_chosen(chm), an_c(:,chm), f1);
            [~,f2] = plotenergy(case1, "resz", o_chosen(chm), an_c(:,chm), f2);
            [~,f3] = plotenergy(case1, "invu", o_chosen(chm), an_c(:,chm), f3);
        end
    else
        if i > 1 
            close([f1 f2 f3]); 
        end
        error('No unstable mode found at lambda = %.3f.',lambda_nd(i));
    end
end
toc(t1);

%% Function defs
% Choose eigenvalue
function [o_chosen, an_c] = choose_eiv(o, an, k)
c = o/k;
a = 1:length(c);
crange = (real(c)>0) & (real(c)<0.5) & (imag(c)<0);
aa = a(crange);
abch = isoutlier(imag(c(aa)), 'movmedian',5);
aa = [a(imag(c)>0) a(real(c)<0) aa(abch)];
o_chosen = o(aa);
an_c = an(:, aa);
end

function [ENG, f1] = plotenergy(case1, engname, o, an, varargin)
[ENG,z] = energy_boomkamp(case1, o, an);
if isempty(varargin)
    f1 = figure('position', [0 0 360 720]);
    plot(ENG.(engname), z, 'color', [0.93,0.00,0.00], 'linewidth', 3, 'DisplayName', num2str(case1.getprop('lambda'),'$%.3f$'));
%     plot(ENG.(engname), z, 'k', 'linewidth', 2, 'DisplayName', num2str(case1.k,'$%.2f$'));
    hold on;
    zc = case1.zc;
    yf(1) = xline(0, '--', 'color', '#bdbdbd');
    yf(2) = yline(zc, '--', 'color', [0.93,0.00,0.00]);
    hold off; box on;
    set(yf(:),'linewidth',2,'HandleVisibility','off');
    if (case1.h > 6)
        blim = -6;
    else
        blim = -case1.h;
    end
    ylabel('$z$', 'rotation', 0, 'HorizontalAlignment', 'right');
    ylim([blim 0]);
    title('$\tau(z)$');
    grid on; box on;
    leg = legend('location', 'southwest');
    title(leg,'$\lambda$','interpreter','latex');
else
    ls = {'--',':','-.'};
    f1 = varargin{1};
    axf1 = get(f1,'CurrentAxes');
    hnum = length(findobj(axf1, 'Type','Line'));
    hold(axf1,'on');
    plot(axf1, ENG.(engname), z, 'color', [0.00,0.57,0.00], 'linewidth', 3, 'DisplayName', num2str(case1.getprop('lambda'),'$%.3f$'));
    yline(axf1,case1.zc, '--', 'color', [0.00,0.57,0.00],'linewidth',2,'HandleVisibility','off');
%     plot(axf1, ENG.(engname), z, 'k', 'linewidth', 2, 'linestyle', ls{hnum}, 'DisplayName', num2str(case1.k,'$%.2f$'));
%     yline(axf1,-case1.zc, 'color', [0.00,0.57,0.00],'linewidth',2.5,'HandleVisibility','off');
    hold(axf1,'off');
end
end