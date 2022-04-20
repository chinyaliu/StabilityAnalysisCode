function [c, an, cA, errGEP, dob] = solver(obj, alg, des, bal, eigspec, funcN, addvar)
% Set initial critical height guess
if isfield(addvar,'zL1')
    obj.zc = addvar.zL1;
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
it = 11; % maximum iteratiom
for i = 1:it
    % Construct matrix A B
    [A,B] = obj.makeAB();
    % De-singularize
    er = -200i;
    if strcmpi(des,'y')
        sinB = all(abs(B)<eps,2);
        B(sinB,:) = A(sinB,:)/er;
    end
    sinA = all(abs(A)<eps,2);
    if sum(sinA) > 0
        A(sinA,:) = B(sinA,:)*er;
    end
    % Solve for the eigenvalue(s)
    if nargout > 2
        [c,an,vout{:}] = solfunc(A,B,eigspec,alg);
    else
        [c,an] = solfunc(A,B,eigspec,alg);
    end
    % Discard artificial eigenvalues(if exists)
%         c = c(real(c)>0.8*er); 
%         an = an(:,real(c)>0.8*er);
    c = c(imag(c)>0.8*er); 
    an = an(:,imag(c)>0.8*er);
    if isempty(c)
        c = nan;
    end
    % Determine if the eigenvalue of largest imaginary part converge
    ztemp = obj.invbf(real(c(1)));
    if i == 1 % given initial value on first iteration
        c_converging = 1;
        c_diff = 100;
    else
        c_converging = abs(c_temp-c(1)) < 1.5*c_diff;
        c_diff = abs(c_temp-c(1));
    end
    c_temp = c(1);
    c_good = imag(c_temp*obj.k)>1e-8 && real(c_temp) > 0 && real(c_temp) < max(obj.ud);
    if(c_diff < 1e-8) % converged
        fprintf('converged to zL = %.8f\n', obj.zc);
        break;
    elseif(c_good && i ~= it && c_converging && -ztemp<obj.h) % keep iterating
        fprintf('iter %2d, zL = %.8f\n', i, obj.zc);
        obj.zc = ztemp;
    else
        obj.zc = nan;
        fprintf('Didn''t converge.\n');
        break;
    end
    % Update the subdomains with eigenvalues of this iteration
    addvar.c = c;
    obj.setsubd('n', funcN, addvar);
end
% Return variables
if strcmpi(bal,'y')
    [dob,errGEP,cA] = deal(vout{:});
else
    [errGEP,cA] = deal(vout{:});
    dob = 'n';
end
end