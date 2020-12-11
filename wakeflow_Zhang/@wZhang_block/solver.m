function [o, an, cA, errGEP, dob] = solver(obj, zL1, iter, alg, bal, eps)
%% Iterate for domain height z_L
obj.cL = eps; obj.zL = zL1;
flag = 0;
[N, arr] = obj.setN4sub();
% [N, arr] = obj.setN2sub();
for j = 1:length(N)
    subD(j) = subdomain(N(j),arr(j),arr(j+1),obj.dm,obj.k);
end
for i = 1:21
      %% Construct matrix A B
    [A, B] = obj.matAB(subD);
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
    if (strcmp(iter,'n')) % known critical height
        break;
    elseif(abs(obj.zL+ztemp) < 1e-9) % converged
        fprintf('converged to zL = %.8f\n', obj.zL);
        break;
    elseif(imag(o(1)) > 0 && i ~= 10) % keep iterating
        fprintf('iter %2d, zL = %.8f\n', i, obj.zL);
        obj.zL = -ztemp;
    elseif (flag == 1 || zL1 == 0.74708299) % didn't converge or growth rate !> 0
        obj.zL = zL1;
        break;
    else % try inflection point z = -0.74708299
        obj.zL = 0.74708299;
        flag = 1;
        fprintf('iter %2d, try inflection pt\n', i);
    end
%     if i == 1
%         obj.N = obj.N/2;
%     end
    newcl(obj,imag(o(1)/obj.k),eps);
    [N, arr] = obj.setN4sub();
    for j = 1:length(N)
        if j > length(subD)
            subD(j) = subdomain(N(j),arr(j),arr(j+1),obj.dm,obj.k);
        else
            subD(j).chgL(N(j),arr(j),arr(j+1));
        end
    end
    subD = subD(1:length(N));
end
obj.zc = -obj.g(real(o)/obj.k);
if nargout > 1 % calculate modeshape only when eigenvector is needed
    obj.z = []; obj.phi = [];
    N = [0 N+1];
    for i = 1:length(N)-1
        obj.z = [obj.z;subD(i).z];
        obj.phi = [obj.phi;subD(i).modeshape(an(sum(N(1:i))+1:sum(N(1:i+1))))];
    end
end
end      