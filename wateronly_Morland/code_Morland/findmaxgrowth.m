function [lambda,maximum_c] = findmaxgrowth(init_var, sol_vars, cut, varargin)
    flow1 = wMorland(init_var{:});
    clim = init_var{3}/2;
    crange = @(c) sqrt((clim-real(c))^2+imag(c)^2);
    [inimin, inimax] = findneutral(init_var{3},init_var{4});
    if isnan(inimin)
        lambda = NaN;
        maximum_c = NaN;
    else
        inrange = @(lam) lam<=inimax & lam>=inimin;
        if length(varargin) == 1
            lamR = 0.5*min(inimax-varargin{1},varargin{1}-inimin);
            lamguess = varargin{1};
        else
            lamR = 0.5*(inimin+inimax);
            lamguess = 0.5*(inimax-inimin);
        end
        sol_vars{end}.zL1 = -case1.criticalH(sqrt(0.5*(lamguess+1./lamguess)));
        [lambda, maximum_c] = maxgrowth(lamguess, lamR);
    end

    function [maxlam, maxc] = maxgrowth(templam, dlam)
        lambda_list = linspace(templam-dlam,templam+dlam,cut);
        c_list = NaN(1,cut);
        % Find eigenvalues of GEP using lambda in lambda_list
        for i = 1:length(lambda_list)
            lam = lambda_list(i);
            fprintf('wavelength = %.5f\n', lam);
            flow1.k = 2*pi/lam;
            flow1.h = 3*lam;
            c = flow1.solvers(sol_vars{:});
            if crange(c) < clim % Howard's semicircle theorm
                c_list(i) = c;
            end
            if ~isnan(flow1.zc)
                sol_vars{end}.zL1 = -flow1.zc;
            end
        end
        % Select the indice of lambda_list with the largest imaginary part of eigenvalue
        [~,ind] = max(imag(c_list)./lambda_list);
        % Decide if lambda converges, if not recursively call this function
        if abs(lambda_list(ind)-templam)<1e-6 && ~isnan(c_list(ind))
            maxlam = lambda_list(ind);
            maxc = c_list(ind);
        else
            if dlam < 1e-10 || ~inrange(lambda_list(ind)) || isnan(c_list(ind))
                maxlam = NaN;
                maxc = NaN;
            else
                new_dlam = lambda_list(2) - lambda_list(1);
                [maxlam, maxc] = maxgrowth(lambda_list(ind), new_dlam);
            end
        end
    end
end