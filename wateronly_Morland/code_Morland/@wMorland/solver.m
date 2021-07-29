function [c, an, cA, errGEP, dob] = solver(obj, alg, bal, eigspec, funcN, addvar)
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

obj.subD = repmat(obj.subDclass(),length(N),1);
for j = 1:length(N)
    obj.subD(j) = obj.subDclass(N(j),arr(j),arr(j+1),obj.dm,obj.k,obj.ud,obj.delta);
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
it = 11; % maximum iteratiom
for i = 1:it
    %% Construct matrix A B
    % Governing equation
    [Age, Bge] = obj.subD(1).match(obj.subD);
    % BC (free surface)
    [Abc1, Bbc1] = obj.subD(1).BC0(size(Age,2)-1);
    % BC (truncated)
    [Abc2, Bbc2] = obj.subD(end).BCh(size(Age,2)-1);
    A = [Age; Abc1; Abc2];
    B = [Bge; Bbc1; Bbc2];
    %% Find the eigenvalue(s)
    if nargout > 2
        [c,an,vout{:}] = solfunc(A,B,eigspec,alg);
    else
        [c,an] = solfunc(A,B,eigspec,alg);
    end
    %% Determine if the eigenvalues converge
    ztemp = obj.criticalH(real(c(1)));
    ztemp_good = isreal(ztemp) && ztemp < 0 && ztemp > -obj.h && imag(c(1)*obj.k)>1e-5;
    if(abs(ztemp+obj.zc) < 1e-6) % converged
        fprintf('converged to zL = %.8f\n', obj.zc);
        break;
    elseif(imag(c(1)) > 0 && i ~= it && ztemp_good)
%     elseif(imag(c(1)) > 0 && i ~= it && ztemp_good && ~strcmp(eigspec,'all'))
        fprintf('iter %2d, zL = %.8f\n', i, obj.zc);
        obj.zc = -ztemp;
    else
        obj.zc = nan;
        fprintf('Didn''t converge.\n');
        break;
    end
    %% Modify the subdomains for the next iteration
    [N, arr] = funcN(obj,'n',addvar,c);
    for j = 1:length(N)
        if j > length(obj.subD)
            obj.subD(j) = obj.subDclass(N(j),arr(j),arr(j+1),obj.dm,obj.k,obj.ud,obj.delta);
            obj.subD(j).makeAB();
        else
            obj.subD(j).chgL(N(j),arr(j),arr(j+1),obj.dm,obj.ud,obj.delta);
        end
    end
    if length(N) < length(obj.subD)
        obj.subD = obj.subD(1:length(N));
    end
end
obj.zc = -obj.zc;
end      