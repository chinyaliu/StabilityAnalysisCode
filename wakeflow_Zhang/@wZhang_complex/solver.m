function [o, an, cA, errGEP, dob] = solver(obj, alg, bal)
%% Differential matrix & Base flow velocity
w = 1./((2*obj.h./((1+obj.zeta).^2)-2*1i*obj.del*obj.zeta).^(0:1:obj.ord));
% w = 1./((obj.h./(obj.zeta+1)-2*1i*obj.del*obj.zeta).^(0:1:obj.ord));
% w = 1./((0.5*obj.h*(1-2*1i*obj.del*obj.zeta)).^(0:1:obj.ord));
% w = 1./((0.5*obj.h-2*1i*obj.del*obj.zeta).^(0:1:obj.ord));
% w = 1./((0.5*obj.h*(1-0.5*pi*1i*obj.del*sinpi(obj.zeta/2))).^(0:1:obj.ord));
obj.D = permute(permute(obj.Din,[1 3 2]).*w, [1 3 2]);
Ubase = obj.baseflow();
%% Construct matrix A B
[A, B] = obj.matAB(Ubase);
%% Find eigenvalue(s)
if strcmp(bal,'y')
    [o,an,dob,errGEP,cA] = balancing(A,B,1,'y',alg);
else
    if nargout > 2
        [o,an,errGEP,cA] = solveGEP(A,B,'all',alg);
        dob = 0;
    else
        [o,an] = solveGEP(A,B,'all',alg);
    end
end
obj.zc = obj.g(real(o(1))/obj.k);
end