function [N,arr] = setN4sub(obj,init,addvar,o,varargin)
if strcmpi(init,'init')
    cL = obj.zc;
else
    c1 = 0.9988; c2 = 0.8814;
    uzz = 2*c1*c2^2*(sech(c2*obj.zc)).^2.*((sech(c2*obj.zc)).^2-2*(tanh(c2*obj.zc)).^2);
    cL = addvar.eps*sqrt(2*abs(imag(o(1)/obj.k)/uzz));
    if isnan(cL)
        cL = obj.zc;
    end
end
if cL > 0.5*obj.h
    arr = [0 obj.zc obj.h];
else
    if cL > 0.5*obj.zc
        arr = [0 obj.zc 2*obj.zc obj.h];
    else
        arr = [0 obj.zc-cL obj.zc obj.zc+cL obj.h];
    end
end
N = zeros(1,length(arr)-1);
a = diff(arr)/obj.h < 0.2;
N(a) = ceil(obj.N/5);
N = N + ceil((obj.N-sum(N))*diff(arr)/obj.h).*(~a);
N(end) = obj.N-sum(N(1:end-1));
end