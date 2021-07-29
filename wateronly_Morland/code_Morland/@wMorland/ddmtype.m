function out = ddmtype(name)
    switch(name)
        case 1
            out = @setN1sub;
        case 2
            out = @setN2sub;
        case 4
            out = @setN4sub;
        otherwise
            error('Function for domain number not specified.');
    end
    % For one domain
    function [N, arr] = setN1sub(obj,varargin)
        arr = [0 obj.h];
        N = obj.N;
    end
    % For two subdomains
    function [N, arr] = setN2sub(obj,varargin)
        arr = [0 obj.zc obj.h];
        N = round(diff(arr)*obj.N/obj.h);
        if sum(N<30)>0
            N(N<30) = 30;
        end
    end
    % For four subdomains
    function [N,arr] = setN4sub(obj,init,addvar,c,varargin)
    if strcmpi(init,'init')
        cL = obj.zc;
    else
        uzz = obj.ud*exp(2*obj.zc/obj.delta)*(2/obj.delta).^2;
        cL = addvar.eps*sqrt(2*abs(imag(c(1))/uzz));
        if isnan(cL)
            cL = obj.zc;
        end
    end
    if cL >= obj.zc
        arr = [0 obj.zc min([2*obj.zc,0.5*(obj.h+obj.zc)]) obj.h];
    else
        if cL >= obj.h-obj.zc
            arr = [0 obj.zc-cL obj.zc obj.h];
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
end