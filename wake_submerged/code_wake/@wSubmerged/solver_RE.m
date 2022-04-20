function [o, an, oa, ana] = solver_RE(obj, alg, des, bal, eigspec, funcN, addvar)
% A solver with eigenvalue choosing criteria (Rayleigh equation only)
eigspec = 'all'; % Whole eigenvalue spectrum is required
% Set initial critical height guess
if isfield(addvar,'zL1')
    obj.zc = addvar.zL1;
else
    obj.zc = 0;
end
ztemp = obj.zc; % given initial value on first iteration
% Set subdomain positions and number of colloaction points
obj.setsubd('init', funcN, addvar);
if obj.subD(end).U(end,1)<eps
    [zz,uu] = obj.getprop('u');
    inds = find(uu(:,1)>eps);
    obj.h = -zz(inds(end));
    obj.setsubd('init', funcN, addvar);
end
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
    if strcmpi(des,'y')
        er = -200;
        sinB = all(B<eps,2);
        B(sinB,:) = A(sinB,:)/er; % Give artificial eigenvalues
    end
    [oa,ana] = solfunc(A,B,eigspec,alg);
    if strcmpi(des,'y')
        % Discard artificial eigenvalues
        oa = oa(real(oa)>0.8*er); 
        ana = ana(:,real(oa)>0.8*er);
    end
    % Choose eigenvalue
    a = 1:length(oa); 
    ca = oa/obj.k;
    crange = ((real(ca)-0.0012>-1e-5) & (real(ca)-1<=1e-7));
    aa = a(crange);
    abch = isoutlier(imag(ca(aa)),'movmedian',5);
    o = oa(aa(abch)); an = ana(:,aa(abch));
    c = ca(aa(abch));
    % Determine if the eigenvalue of largest growth rate converge
    if isempty(o)
        obj.zc = nan;
        fprintf('Didn''t converge.\n');
        break;
    else      
        [~,ind] = sort(imag(o),'descend');   
        c = c(ind); an = an(:,ind); o = o(ind);
        z_diff = abs(obj.invbf(real(c(1)))-ztemp);
        ztemp = obj.invbf(real(c(1)));
        if(z_diff < 1e-10) % converged
            fprintf('converged to zL = %.8f\n', obj.zc);
            break;
        elseif(i ~= it && -ztemp<obj.h && isreal(ztemp)) % keep iterating
            fprintf('iter %2d, zL = %.8f\n', i, obj.zc);
            obj.zc = ztemp;
        else
            obj.zc = nan;
            fprintf('Didn''t converge.\n');
            break;
        end
    end
    % Update the subdomains with eigenvalues of this iteration
    addvar.c = c;
    obj.setsubd('n', funcN, addvar);
end
o = [o; oa(a(~crange))];
an = [an ana(:,a(~crange))];
end