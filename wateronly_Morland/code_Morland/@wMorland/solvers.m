function [c, an, cA, errGEP, dob] = solvers(obj, alg, des, bal, eigspec, funcN, addvar)
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
    vout = cell(3,1);
else
    solfunc = @solveGEP;
    vout = cell(2,1);
end
it = 11; % maximum iteratiom
c_diff = 100;
c_temp = 0;
for i = 1:it
    %% Construct matrix A B
    % Governing equation
    [Age, Bge] = obj.subD(1).match(obj.subD);
    % BC (free surface)
    [Abc1, Bbc1] = obj.subD(1).BC0(size(Age,2)-1);
   % BC (free slip)
   [Abc2, Bbc2] = obj.subD(end).BCh(size(Age,2)-1);
   A = [Age; Abc1; Abc2];
   B = [Bge; Bbc1; Bbc2];
%     % BC (exponential decay)
%     [Abc2, Bbc2] = obj.subD(end).BCh3(size(Age,2)-1);
%     A = [Age; Abc1; Abc2];
%     B = [Bge; Bbc1; Bbc2];    
      %% De-singularize
    if strcmpi(des,'y')
        er = -200i;
        sinB = all(B<eps,2);
        B(sinB,:) = A(sinB,:)/er;
    end
    %% Find the eigenvalue(s)
    if nargout > 2
        [c,an,vout{:}] = solfunc(A,B,eigspec,alg);
    else
        [c,an] = solfunc(A,B,eigspec,alg);
    end
    %% Determine if the eigenvalues converge
    ztemp = obj.criticalH(real(c(1)));
    c_converging = i==1 || (abs(c_temp-c(1)) < 1.5*c_diff);
    c_diff = abs(c_temp-c(1));
    c_temp = c(1);
    c_good = imag(c_temp*obj.k)>1e-5 && real(c_temp) > 0 && real(c_temp) < obj.ud;
    if(c_diff < 1e-8) % converged
        fprintf('converged to zL = %.8f\n', obj.zc);
        break;
    elseif(c_good && i ~= it && c_converging && ztemp<obj.h)
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
if strcmpi(bal,'y')
    [dob,errGEP,cA] = deal(vout{:});
else
    [errGEP,cA] = deal(vout{:});
end
obj.zc = -obj.zc;
end      