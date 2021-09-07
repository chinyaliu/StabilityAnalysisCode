function [lambda,maximum_c] = findmaxgrowth(init_var, sol_var, cut)
    flow1 = wMorland(init_var{:});
    clim = init_var{3}/2;
    crange = @(c) sqrt((clim-real(c))^2+imag(c)^2);
    [inimin, inimax] = findneutral(init_var{3},init_var{4});
    if isnan(inimin)
        lambda = NaN;
        maximum_c = NaN;
    else
        inrange = @(lam) lam<=inimax & lam>=inimin;
        [lambda, maximum_c] = maxgrowth(0.5*(inimin+inimax), 0.5*(inimax-inimin));
    end

    function [maxlam, maxc] = maxgrowth(templam, dlam)
        lambda_list = linspace(templam-dlam,templam+dlam,cut);
        c_list = NaN(1,cut);
        for i = 1:length(lambda_list)
            lam = lambda_list(i);
            fprintf('wavelength = %.5f\n', lam);
            flow1.k = 2*pi/lam;
%             flow1.h = 2*lam;
%             flow1.h = max(lam,3*init_var{4});
            flow1.h = 3*init_var{4};
            c = flow1.solver(sol_var{:});
            if crange(c) < clim
                c_list(i) = c;
            end
            if ~isnan(flow1.zc)
                sol_var{end}.zL1 = -flow1.zc;
            end
        end
        [~,ind] = max(imag(c_list)./lambda_list);
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
    function [nmin, nmax] = findneutral(ud,dlt)
        polylam = [4*ud^2+pi*dlt 0 -4*pi*dlt*ud^4+4*ud^2+2*pi*dlt 0 pi*dlt];
        lamall = roots(polylam);
        lam = real(lamall(imag(lamall)<eps & real(lamall)>eps));
        if length(lam) == 2
            nmin = min(lam);
            nmax = max(lam);
%         elseif length(lam) == 1 && lam(1) ~= 0
%             nmin = 0;
%             nmax = lam;
        else
            nmin = NaN;
            nmax = NaN;
        end
    end
end