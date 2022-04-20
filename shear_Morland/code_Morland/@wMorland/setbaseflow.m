function setbaseflow(obj,bf)
    switch(lower(bf))
    case 'exponential'
        % exponential profile
        [obj.baseflow,obj.invbf] = exponential();
    case 'error function'
        % error functiona profile
        [obj.baseflow,obj.invbf] = errfunction();
    otherwise
        error('Invalid velocity profile name.');
    end

    function [bff,ibf] = exponential()
        bff = @baseflow;
        ibf = @invbf;
        function Ur = baseflow(z)
            Ur(:,1) = obj.ud*exp(2*z/obj.delta);
            Ur(:,2) = (2/obj.delta)*Ur(:,1);
            Ur(:,3) = (2/obj.delta)*Ur(:,2);
        end
        function out = invbf(x)
            out = 0.5*obj.delta*log(x/obj.ud);
        end
    end
    function [bff,ibf] = errfunction()
        bff = @baseflow;
        ibf = @invbf;
        function Ur = baseflow(z)
            Ur(:,1) = obj.ud*erfc(-2*z/obj.delta/sqrt(pi));
            Ur(:,2) = 4*obj.ud*exp(-4*z.^2/pi/obj.delta^2)/pi/obj.delta;
            Ur(:,3) = Ur(:,2).*(-8*z/pi/obj.delta^2);
        end
        function out = invbf(x)
            out = -erfcinv(x/obj.ud)*sqrt(pi)*obj.delta/2;
        end
    end
end     
