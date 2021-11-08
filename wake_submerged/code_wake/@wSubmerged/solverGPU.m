function [o, an, cA, errGEP, dob] = solverGPU(obj, alg, des, bal, eigspec, funcN, addvar)
% Set initial critical height guess
if isfield(addvar,'zL1')
    obj.zc = max([addvar.zL1 2*obj.H-addvar.zL1]);
else
    obj.zc = 0;
end
% Set subdomain positions and number of colloaction points
obj.setsubd('init', funcN, addvar);
% Choose whether to solve with balancing
it = 11; % maximum iteratiom
c_temp = 0;
for i = 1:it
    % Construct matrix A B
    [a,b] = obj.makeAB();
    A = gpuArray(a);
    B = gpuArray(b);
    % De-singularize
    if strcmpi(des,'y')
        er = -200;
        sinB = all(B<eps,2);
        B(sinB,:) = A(sinB,:)/er;
    end
    % Solve for the eigenvalue(s)
    if nargout > 2
        [o,an,vout{:}] = solveGEP(A,B,eigspec,alg);
    else
        [o,an] = solveGEP(A,B,eigspec,alg);
    end
%    o = gather(o);
    % Determine if the eigenvalue of largest imaginary part converge
    c = o(1)/obj.k;
    c_good = imag(c(1)*obj.k)>1e-8 && real(c(1)) > 0 && real(c(1)) < 1;
    if(abs(c_temp-c(1)) < 1e-8) % converged
        fprintf('converged to zL = %.8f\n', obj.zc);
        break;
    elseif(c_good && i ~= it && obj.criticalH(real(c(1)))<obj.h) % keep iterating
        fprintf('iter %2d, zL = %.8f\n', i, obj.zc);
        c_temp = c(1);
        obj.zc = obj.criticalH(real(c(1)));
    else
        obj.zc = nan;
        fprintf('Didn''t converge.\n');
        break;
    end
    % Update the subdomains with eigenvalues of this iteration
    addvar.c = c;
    obj.setsubd('n', funcN, addvar);
end
o = gather(o);
an = gather(an);
end
