function out = ddmtype(name)
    switch(name)
        case 1
            out = @setN1sub;
        case 2
            out = @setN2sub;
        case 41
            out = @setN4test1;
        case 44
            out = @setN4test4;
        case 5
            out = @setN5sub;
        case 9
            out = @setNnsub;
        case 92
            out = @setNn2sub;
        otherwise
            error('DDM not specified.');
    end
    % For one domain
    function [N, arr] = setN1sub(obj,varargin)
        arr = [0 obj.h];
        N = obj.N;
    end
    % For two subdomains
    function [N, arr] = setN2sub(obj,varargin)
        arr = [0 -obj.zc obj.h];
        N = round(diff(arr)*obj.N/obj.h);
        if sum(N<30)>0
            N(N<30) = 30;
        end
    end
    % Test first domain
    function [N,arr] = setN4test1(obj,init,addvar,varargin)
        if -obj.zc < obj.h
            if strcmpi(init,'init')
                cL = min(-obj.zc/2,abs(obj.h+obj.zc)/2);
            else
                uzz = obj.baseflow(obj.zc);
                cL = addvar.eps*sqrt(2*abs(imag(addvar.c(1))/uzz(3)));
                cL = min([(obj.h+obj.zc)/2, obj.zc/2, cL]);
            end
            arr = [0 -obj.zc-cL -obj.zc -obj.zc+cL obj.h];
            Nmin = 50;
            N = [obj.N-4*Nmin Nmin Nmin 2*Nmin];
        else
            [N, arr] = setN1sub(obj); % Or throw error
        end
        mustBePositive(N);
        mustBePositive(diff(arr));
    end
    % Test fourth domain
    function [N,arr] = setN4test4(obj,init,addvar,varargin)
        if -obj.zc < obj.h
            if strcmpi(init,'init')
                cL = min(-obj.zc/2,abs(obj.h+obj.zc)/2);
            else
                uzz = obj.baseflow(obj.zc);
                cL = addvar.eps*sqrt(2*abs(imag(addvar.c(1))/uzz(3)));
                cL = min([(obj.h+obj.zc)/2, -obj.zc/2, cL]);
            end
            arr = [0 -obj.zc-cL -obj.zc -obj.zc+cL obj.h];
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
    % For five subdomains
    function [N,arr] = setN5sub(obj,init,addvar,varargin)
        if -obj.zc < obj.h
            if strcmpi(init,'init')
                cL = min(-obj.zc/2,abs(obj.h+obj.zc)/2);
            else
                uzz = obj.baseflow(obj.zc);
                cL = addvar.eps*sqrt(2*abs(imag(addvar.c(1))/uzz(3)));
                cL = min([(obj.h+obj.zc)/2, -obj.zc/2, cL]);
            end
            arr = [0 -obj.zc-cL -obj.zc -obj.zc+cL 0.5*(obj.h-obj.zc+cL) obj.h];
            Nmin = 50;
            N = [2*Nmin Nmin Nmin round((obj.N-4*Nmin)/2) round((obj.N-4*Nmin)/2)];
        else
            [N, arr] = setN1sub(obj); % Or throw error
        end
        mustBePositive(N);
        mustBePositive(diff(arr));
    end
    % automatically cut domain (for exponential profile)
    function [N,arr] = setNnsub(obj,init,addvar,varargin)
        if strcmpi(init,'init')
            cL = min(obj.zc/2,abs(obj.h+obj.zc)/2);
        else
            uzz = obj.baseflow(obj.zc);
            cL = addvar.eps*sqrt(2*abs(imag(addvar.c(1))/uzz(3)));
            cL = min([(obj.h+obj.zc)/2, -obj.zc/2, cL]);
        end
        Nmin = 50;
        zdecayu = obj.delta/2*log(10);
        zin = flip(obj.h:-zdecayu:-obj.zc+cL);
        arr = [0 -obj.zc-cL -obj.zc -obj.zc+cL zin(2:end-1) obj.h];
        cutarr = diff(arr)/(obj.h+obj.zc-cL);
        N = [2*Nmin Nmin Nmin round((obj.N-4*Nmin)*cutarr(4:end))];
        N(end) = N(end) + obj.N - sum(N);
        
        mustBePositive(N);
    end
    % automatically cut domain (for exponential profile)
    function [N,arr] = setNn2sub(obj,init,addvar,varargin)
        if strcmpi(init,'init')
            cL = min(-obj.zc/2,abs(obj.h+obj.zc)/2);
        else
            uzz = obj.baseflow(obj.zc);
            cL = addvar.eps*sqrt(2*abs(imag(addvar.c(1))/uzz(3)));
            cL = min([(obj.h+obj.zc)/2, -obj.zc/2, cL]);
        end
        Nmin = 50;
        zdecayu = obj.delta*log(10);
        zin = flip(obj.h:-zdecayu:-obj.zc+cL);
        arr = [0 -obj.zc-cL -obj.zc -obj.zc+cL zin(2:end-1) obj.h];
        cutarr = diff(arr)/(obj.h+obj.zc-cL);
        N = [2*Nmin Nmin Nmin round((obj.N-4*Nmin)*cutarr(4:end))];
        N(end) = N(end) + obj.N - sum(N);
        
        mustBePositive(N);
    end
end