function setbaseflow(obj,bf)
    switch(lower(bf))
    case 'cosh'
        [obj.baseflow,obj.invbf] = hydrofoil();
    case 'pwlinear'
        [obj.baseflow,obj.invbf] = pwlinear();
    case 'tanh'
        [obj.baseflow,obj.invbf] = cylintan();
    otherwise
        error('Invalid velocity profile name.');
    end

    function [bff,ibf] = hydrofoil()
        bff = @baseflow;
        ibf = @invbf;
        function Ur = baseflow(z)
            c1 = 0.9988; c2 = 0.8814;
            Ur(:,1) = 1-c1*cosh(c2*(z+obj.H)).^(-2);
            Ur(:,2) = 2*c1*c2*tanh(c2*(z+obj.H)).*(sech(c2*(z+obj.H))).^2;
            Ur(:,3) = 2*c1*c2^2*(sech(c2*(z+obj.H))).^2.*((sech(c2*(z+obj.H))).^2-2*(tanh(c2*(z+obj.H))).^2);
        end
        function out = invbf(x)
            out = obj.H+abs((5000*acosh((-2497./(625*(x-1))).^(1/2)/2))./4407);
            out(x <= 0.0012) = nan;
            out(x > 1) = nan;
        end
    end
    function [bff,ibf] = pwlinear()
        bff = @baseflow;
        ibf = @invbf;
        function Ur = baseflow(z)
            if obj.H > 0
                error('Profile for half-submerged only');
            end
            u0 = 0.0012; h0 = 0.25; H0 = 1.75;
            Ur = [ones(length(z),1) zeros(length(z),2)];
            Ur(z<H0,1) = u0 + (1-u0)*(z-h0)/(H0-h0);
            Ur(z<h0,1) = u0;
            Ur((z<H0)&(z<h0),2) = (1-u0)/(H0-h0);
        end
        function out = invbf(x)
            u0 = 0.0012; h0 = 0.25; H0 = 1.75;
            if x <= u0
                out = 0;
            elseif x >= 1
                out = H0;
            else
                out = h0+(x-u0)*(H0-h0)/(1-u0);
            end
        end
    end
    function [bf,ibf] = cylintan()
        bf = @baseflow;
        ibf = @invbf;
        function Ur = baseflow(z)
            if obj.H > 0
                error('Profile for half-submerged only');
            end
            c = 0.75; alpha = 4; beta = 0.32;
            Ur(:,1) = 1-c*(1-tanh(alpha*z.^2-beta));
            Ur(:,2) = 2*c*alpha*z.*(sech(beta-alpha*z.^2)).^2;
            Ur(:,3) = 2*c*alpha*(4*alpha*z.^2.*tanh(beta-alpha*z.^2)+1).*(sech(beta-alpha*z.^2)).^2;
        end
        function out = invbf(x)
            c = 0.75; alpha = 4; beta = 0.32;
            out = (beta+atanh((x-1)/c+1)).^0.5/alpha.^0.5;
        end
    end
end     
