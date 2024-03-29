function [o, an, cA, errGEP, dob] = solver(obj, alg, des, bal, eigspec, funcN, addvar)
if isempty(obj.method)
    error('Numerical method not specified.');
end% Set initial critical height guess
if isfield(addvar,'zL1')
    obj.zc = max([addvar.zL1 2*obj.H-addvar.zL1]);
else
    obj.zc = 0;
end
% Set subdomain positions and number of colloaction points
obj.setsubd('init', funcN, addvar);
% Choose whether to solve with balancing
if strcmpi(bal,'y')
    solfunc = @balanceAB;
    vout = cell(3,1);
else
    solfunc = @solveGEP;
    vout = cell(2,1);
end
it = 11; % maximum iteration
c_temp = 0;
for i = 1:it
    % Construct matrix A B
    [A,B] = obj.makeAB();
    % De-singularize
    er = -200;
    if strcmpi(des,'y')
        sinB = all(B<eps,2);
        B(sinB,:) = A(sinB,:)/er;
    end
    % Solve for the eigenvalue(s)
    if nargout > 2
        [o,an,vout{:}] = solfunc(A,B,eigspec,alg);
    else
        [o,an] = solfunc(A,B,eigspec,alg);
    end
    % Discard artificial eigenvalues(if exists)
    o = o(real(o)>0.8*er); 
    an = an(:,real(o)>0.8*er);
    % Determine if the eigenvalue of largest imaginary part converge
    c = o(1)/obj.k;
    c_good = imag(c(1)*obj.k)>1e-8 && real(c(1)) > 0 && real(c(1)) < 1;
    if(abs(c_temp-c(1)) < 1e-8) % converged
%         fprintf('converged to zL = %.8f\n', obj.zc);
        break;
    elseif(c_good && i ~= it && obj.invbf(real(c(1)))<obj.h) % keep iterating
%         fprintf('iter %2d, zL = %.8f\n', i, obj.zc);
        c_temp = c(1);
        obj.zc = obj.invbf(real(c(1)));
    else
        obj.zc = nan;
%         fprintf('Didn''t converge.\n');
        break;
    end
    % Update the subdomains with eigenvalues of this iteration
    addvar.c = c;
    obj.setsubd('n', funcN, addvar);
end
% Return variables
if nargout > 2
    if strcmpi(bal,'y')
        [dob,errGEP,cA] = deal(vout{:});
    else
        [errGEP,cA] = deal(vout{:});
        dob = 'n';
    end
end
end