function [k, an, cA, errGEP, dob] = solver(obj, alg, bal)
%% Differential matrix & Base flow velocity
Ubase = obj.baseflow();
%% Construct matrix A
A = obj.matAB(obj.Din, Ubase);
%% Find eigenvalue(s)
if strcmp(bal,'y')
    %%[o,an,dob,errGEP,cA] = balancing(A,B,1,'y',alg);
else
    dob = 0;
    if nargout > 2
        [k,an,errGEP,cA] = solvePoly(A,'all',alg);
    else
        [k,an] = solvePoly(A,'all',alg);
    end
end
% obj.zc = -obj.g(obj.o/real(k));
end      