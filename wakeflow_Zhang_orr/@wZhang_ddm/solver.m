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
end
it = 21;
for i = 1:it
      %% Construct matrix A B
    % Governing equation
    [A, B] = obj.subD(1).match(obj.subD);
    % BC (free surface)
    [A(end-2:end-1,:), B(end-2:end-1,:)] = obj.subD(1).BC0(obj.Fr2,size(A,2)-1);
    % BC (truncated)
    A(end, obj.N-obj.subD(end).N+length(obj.subD):end-1) = obj.subD(end).BCh();    
      %% Find eigenvalue(s)
    if strcmp(bal,'y')
        [o,an,dob,errGEP,cA] = balancing(A,B,1,'y',alg);
    else
        if nargout > 2
            [o,an,errGEP,cA] = solveGEP(A,B,'max',alg);
            dob = 0;
        else
            [o,an] = solveGEP(A,B,'max',alg);
        end
    end
    ztemp = -obj.g(real(o(1))/obj.k);
    if(abs(ztemp+obj.zc) < 1e-9) % converged
        fprintf('converged to zL = %.8f\n', obj.zc);
        break;
    elseif(imag(o(1)) > 0 && i ~= it) % keep iterating
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
            A = blkdiag(A,0); B = blkdiag(B,0);
        else
            obj.subD(j).chgL(N(j),arr(j),arr(j+1));
        end
    end
    if length(N) < length(obj.subD)
        obj.subD = obj.subD(1:length(N));
        A = A(1:end-1,1:end-1); B = B(1:end-1,1:end-1);
    end
end
obj.zc = -obj.zc;
end      