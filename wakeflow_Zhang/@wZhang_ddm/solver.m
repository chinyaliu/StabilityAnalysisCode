function [o, an, cA, errGEP, dob] = solver(obj, alg, bal, funcN, addvar)
%% Iterate for domain height zc
if isfield(addvar,'zL1')
    obj.zc = addvar.zL1;
else
    obj.zc = 0;
end
[N, arr] = funcN(obj,'init',addvar,0);
obj.subD = repmat(subdomain(),length(N),1);
for j = 1:length(N)
    obj.subD(j) = subdomain(N(j),arr(j),arr(j+1),obj.dm,obj.k);
    obj.subD(j).makeAB();
end
it = 21;
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
      %% Find eigenvalue(s)
    switch bal
        case 'y'
            [o,an,dob,errGEP,cA] = balancing(A,B,1,'y',alg);
        case 'all'
            [o,an] = solveGEP(A,B,'all',alg);
        otherwise
            if nargout > 2
                [o,an,errGEP,cA] = solveGEP(A,B,'max',alg);
                dob = 0;
            else
                [o,an] = solveGEP(A,B,'max',alg);
            end
    end
    ztemp = -obj.g(real(o(1))/obj.k);
    if(abs(ztemp+obj.zc) < 1e-8) % converged
        fprintf('converged to zL = %.8f\n', obj.zc);
        break;
    elseif(imag(o(1)) > 0 && i ~= it && ~strcmp(bal,'all')) % keep iterating
        fprintf('iter %2d, zL = %.8f\n', i, obj.zc);
        obj.zc = -ztemp;
    else
        fprintf('Didn''t converge.\n');
        break;
    end
    [N, arr] = funcN(obj,'n',addvar,o);
    for j = 1:length(N)
        if j > length(obj.subD)
            obj.subD(j) = subdomain(N(j),arr(j),arr(j+1),obj.dm,obj.k);
            obj.subD(j).makeAB();
        else
            obj.subD(j).chgL(N(j),arr(j),arr(j+1));
        end
    end
    if length(N) < length(obj.subD)
        obj.subD = obj.subD(1:length(N));
    end
end
obj.zc = -obj.zc;
end      