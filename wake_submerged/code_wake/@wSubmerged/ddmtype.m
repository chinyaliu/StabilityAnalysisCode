function out = ddmtype(name)
    switch(name)
        case 1
            out = @setN1sub;
        case 2
            out = @setN2sub;
        case 3
            out = @setN3sub;
        case 4
            out = @setN4sub;
        case 44
            out = @setN4test4;
        case 45
            out = @setN4test5;
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
        arr = [0 obj.H obj.h];
        N = round(diff(arr)*obj.N/obj.h);
        if sum(N<30)>0
            N(N<30) = 30;
        end
    end
    % For three subdomains
    function [N,arr] = setN3sub(obj,varargin)
        if obj.zc > obj.h
            arr = [0 obj.H obj.h];
        elseif obj.zc-obj.H < 0 || 2*obj.H-obj.zc < 0
            arr = [0 obj.zc obj.h];
        else
            arr = [0 2*obj.H-obj.zc obj.zc obj.h];
        end
        N = zeros(1,length(arr)-1);
        b = ceil(obj.N/5);
        Nmin = b*ones(1,length(N));
        N = Nmin + (obj.N-b*length(N))*round(diff(arr)/obj.h);
        try
            mustBePositive(N);
            mustBePositive(diff(arr));
        catch
            [N, arr] = setN2sub(obj);
        end
    end

    % For four subdomains
    function [N,arr] = setN4sub(obj,varargin)
        Nmin = min(100,ceil(obj.N/5));
        if obj.zc > obj.h
            arr = [0 obj.H obj.h];
            N = [Nmin Nmin] + (obj.N-2*Nmin)*round(diff(arr)/obj.h);
        else
            N = [ceil((obj.N-2*Nmin)*obj.H/obj.h) Nmin Nmin];
            N(4) = obj.N-sum(N);
            if obj.zc-obj.H < 0
                arr = [0 obj.zc obj.H 2*obj.H-obj.zc obj.h];
            else
                arr = [0 2*obj.H-obj.zc obj.H obj.zc obj.h];
            end
        end
        mustBePositive(N);
        mustBePositive(diff(arr));
    end

    function [N,arr] = setN4test4(obj,init,addvar,varargin)
        if obj.zc < obj.h
            if strcmpi(init,'init')
                cL = min(obj.zc/2,abs(obj.h-obj.zc)/2);
            else
                uzz = obj.subDclass().baseflow(obj.zc,obj.H);
                cL = addvar.eps*sqrt(2*abs(imag(addvar.c(1))/uzz(3)));
                cL = min([(obj.h-obj.zc)/2, obj.zc/2, cL]);
            end
            arr = [0 obj.zc-cL obj.zc obj.zc+cL obj.h];
            Nmin = 50;
            N = [2*Nmin Nmin Nmin obj.N-4*Nmin];
        else
            [N, arr] = setN1sub(obj); % Or throw error
        end
        mustBePositive(N);
        mustBePositive(diff(arr));
    end
    function [N,arr] = setN4test5(obj,init,addvar,varargin)
        if obj.zc < obj.h
            if strcmpi(init,'init')
                cL = min(obj.zc/2,abs(obj.h-obj.zc)/2);
            else
                uzz = obj.subDclass().baseflow(obj.zc,obj.H);
                cL = addvar.eps*sqrt(2*abs(imag(addvar.c(1))/uzz(3)));
                cL = min([(obj.h-obj.zc)/2, obj.zc/2, cL]);
            end
            arr = [0 obj.zc-cL obj.zc obj.zc+cL obj.h];
            Nmin = 50;
            if obj.N < 300
                Nmin = round(obj.N/5);
            end
            N = [2*Nmin Nmin Nmin obj.N-4*Nmin];
        else
            [N, arr] = setN1sub(obj); % Or throw error
        end
        mustBePositive(N);
        mustBePositive(diff(arr));
    end
end