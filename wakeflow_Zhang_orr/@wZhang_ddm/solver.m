function [o, an, cA, errGEP, dob] = solver(obj, alg, bal, eigspec, funcN, addvar)
if isempty(obj.method)
    error('Numerical method not specified.');
end
if isfield(addvar,'zL1')
    obj.zc = addvar.zL1;
else
    obj.zc = 0;
end
%% Create subdomains according to funcN
[N, arr] = funcN(obj,'init',addvar,0);
if isinf(obj.Re)
    subBLK = @subRay;
else
    subBLK = @subOrr;
end
obj.subD = repmat(subBLK(),length(N),1);
for j = 1:length(N)
    obj.subD(j) = subBLK(N(j),arr(j),arr(j+1),obj.dm,obj.k,obj.Re);
    obj.subD(j).makeAB();
end
%% Choose whether to solve with balancing
dob = 0; errGEP = []; cA = [];
if strcmpi(bal,'y')
    solfunc = @balanceAB;
    vout = {dob,errGEP,cA};
else
    solfunc = @solveGEP;
    vout = {errGEP,cA};
end
it = 21; % maximum iteratiom
for i = 1:it
    %% Construct matrix A B
    % Governing equation
    [Age, Bge] = obj.subD(1).match(obj.subD);
    % BC (free surface)
    [Abc1, Bbc1] = obj.subD(1).BC0(obj.Fr2,size(Age,2)-1);
    % BC (truncated)
    [Abc2, Bbc2] = obj.subD(end).BCh(size(Age,2)-1);
    A = [Age; Abc1; Abc2];
    B = [Bge; Bbc1; Bbc2];
    %% Find the eigenvalue(s)
    if nargout > 2
        [o,an,vout{:}] = solfunc(A,B,eigspec,alg);
    else
        [o,an] = solfunc(A,B,eigspec,alg);
    end
    %% Determine if the eigenvalues converge
    ztemp = obj.criticalH(real(o(1))/obj.k);
    if(abs(ztemp+obj.zc) < 1e-8) % converged
        fprintf('converged to zL = %.8f\n', obj.zc);
        break;
    elseif(imag(o(1)) > 0 && i ~= it && real(o(1))>0) % keep iterating
%     elseif(imag(o(1)) > 0 && i ~= it && ~strcmp(eigspec,'all') && real(o(1))>0) % keep iterating
        fprintf('iter %2d, zL = %.8f\n', i, obj.zc);
        obj.zc = -ztemp;
    else
        obj.zc = nan;
        fprintf('Didn''t converge.\n');
        break;
    end
    %% Modify the subdomains for the next iteration
    [N, arr] = funcN(obj,'n',addvar,o);
    for j = 1:length(N)
        if j > length(obj.subD)
            obj.subD(j) = subBLK(N(j),arr(j),arr(j+1),obj.dm,obj.k,obj.Re);
            obj.subD(j).makeAB();
        else
            obj.subD(j).chgL(N(j),arr(j),arr(j+1),obj.dm);
        end
    end
    if length(N) < length(obj.subD)
        obj.subD = obj.subD(1:length(N));
    end
end
obj.zc = -obj.zc;
end      