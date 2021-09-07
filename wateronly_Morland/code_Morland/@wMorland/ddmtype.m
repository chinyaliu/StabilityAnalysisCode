function out = ddmtype(name)
    switch(name)
        case 1
            out = @setN1sub;
        case 2
            out = @setN2sub;
        case 4
            out = @setN4sub;
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
            uzz = obj.subDclass().baseflow(-obj.zc,obj.ud,obj.delta);
            cL = addvar.eps*sqrt(2*abs(imag(c(1))/uzz(3)));
            if isnan(cL) || cL > (obj.h-obj.zc)/2 ||  cL > obj.zc/2
                cL = min((obj.h-obj.zc)/2,obj.zc/2);
            end
        end
        if cL > obj.h
            arr = [0 obj.h];
        elseif cL >= obj.zc
            arr = [0 obj.zc min([2*obj.zc,0.5*(obj.h+obj.zc)]) obj.h];
        else
            if cL >= obj.h-obj.zc
                arr = [0 obj.zc-cL obj.zc obj.h];
            else
                arr = [0 obj.zc-cL obj.zc obj.zc+cL obj.h];
            end
        end
        N = zeros(1,length(arr)-1);
        b = ceil(obj.N/5);
        Nmin = [2*b b*ones(1,length(N)-1)];
        N = Nmin + (obj.N-b*(length(N)+1))*round(diff(arr)/obj.h);

        mustBePositive(N);
    end
    % Test first domain
    function [N,arr] = setN4test1(obj,init,addvar,c,varargin)
        if strcmpi(init,'init')
            cL = obj.zc/2;
        else
            uzz = obj.subDclass().baseflow(-obj.zc,obj.ud,obj.delta);
            cL = addvar.eps*sqrt(2*abs(imag(c(1))/uzz(3)));
            if isnan(cL) || cL > (obj.h-obj.zc)/2 ||  cL > obj.zc/2
                cL = min((obj.h-obj.zc)/2,obj.zc/2);
            end
        end
        arr = [0 obj.zc-cL obj.zc obj.zc+cL obj.h];
    %     Nmin = min(100,round(obj.N/4));
        Nmin = 50;
        N = [obj.N-4*Nmin Nmin Nmin 2*Nmin];
        mustBePositive(N);
    end
    % Test fourth domain
    function [N,arr] = setN4test4(obj,init,addvar,c,varargin)
        if obj.zc < obj.h
            if strcmpi(init,'init')
                cL = min(obj.zc/2,abs(obj.h-obj.zc)/2);
            else
                uzz = obj.subDclass().baseflow(-obj.zc,obj.ud,obj.delta);
                cL = addvar.eps*sqrt(2*abs(imag(c(1))/uzz(3)));
                if isnan(cL) || cL > (obj.h-obj.zc)/2 ||  cL > obj.zc/2
                    cL = min((obj.h-obj.zc)/2,obj.zc/2);
                end
            end
            arr = [0 obj.zc-cL obj.zc obj.zc+cL obj.h];
            Nmin = 50;
            N = [2*Nmin Nmin Nmin obj.N-4*Nmin];
        else
            [N, arr] = setN1sub(obj);
        end
        mustBePositive(N);
        mustBePositive(diff(arr));
    end
    % For five subdomains
    function [N,arr] = setN5sub(obj,init,addvar,c,varargin)
        if strcmpi(init,'init')
            cL = obj.zc/2;
        else
            uzz = obj.subDclass().baseflow(-obj.zc,obj.ud,obj.delta);
            cL = addvar.eps*sqrt(2*abs(imag(c(1))/uzz(3)));
            if isnan(cL) || cL > (obj.h-obj.zc)/2 ||  cL > obj.zc/2
                cL = min((obj.h-obj.zc)/2,obj.zc/2);
            end
        end
        arr = [0 obj.zc-cL obj.zc obj.zc+cL 0.5*(obj.h+obj.zc+cL) obj.h];
        Nmin = 50;
        N = [2*Nmin Nmin Nmin round((obj.N-4*Nmin)/2) round((obj.N-4*Nmin)/2)];
        mustBePositive(N);
    end
    % automatically cut domain (for exponential profile)
    function [N,arr] = setNnsub(obj,init,addvar,c,varargin)
        if strcmpi(init,'init')
            cL = obj.zc/2;
        else
            uzz = obj.subDclass().baseflow(-obj.zc,obj.ud,obj.delta);
            cL = addvar.eps*sqrt(2*abs(imag(c(1))/uzz(3)));
            if isnan(cL) || cL > (obj.h-obj.zc)/2 ||  cL > obj.zc/2
                cL = min((obj.h-obj.zc)/2,obj.zc/2);
            end
        end
        Nmin = 50;
        zdecayu = obj.delta/2*log(10);
        zin = flip(obj.h:-zdecayu:obj.zc+cL);
        arr = [0 obj.zc-cL obj.zc obj.zc+cL zin(2:end-1) obj.h];
        cutarr = diff(arr)/(obj.h-obj.zc-cL);
        N = [2*Nmin Nmin Nmin round((obj.N-4*Nmin)*cutarr(4:end))];
        N(end) = N(end) + obj.N - sum(N);
        
        mustBePositive(N);
    end
    % automatically cut domain (for exponential profile)
    function [N,arr] = setNn2sub(obj,init,addvar,c,varargin)
        if strcmpi(init,'init')
            cL = obj.zc/2;
        else
            uzz = obj.subDclass().baseflow(-obj.zc,obj.ud,obj.delta);
            cL = addvar.eps*sqrt(2*abs(imag(c(1))/uzz(3)));
            if isnan(cL) || cL > (obj.h-obj.zc)/2 ||  cL > obj.zc/2
                cL = min((obj.h-obj.zc)/2,obj.zc/2);
            end
        end
        Nmin = 50;
        zdecayu = obj.delta*log(10);
        zin = flip(obj.h:-zdecayu:obj.zc+cL);
        arr = [0 obj.zc-cL obj.zc obj.zc+cL zin(2:end-1) obj.h];
        cutarr = diff(arr)/(obj.h-obj.zc-cL);
        N = [2*Nmin Nmin Nmin round((obj.N-4*Nmin)*cutarr(4:end))];
        N(end) = N(end) + obj.N - sum(N);
        
        mustBePositive(N);
    end
end