function [x_max] = bisearchnan(xL,xU,f,varargin)
%
% Bi-section search method for local mimimum of f(x)
% init     : set to TRUE for first iteration
% [xL, xU] : interval to define the local maximum
% f        : function handle f(x,varargin)
% varargin : additional parameters required to define f(x)
%
d = @(xl,xu) 0.5*(xu+xl); 
err = 1e-5;
x1 = xL;
x2 = xU;
xm = d(xL,xU);
f1 = f(x1,varargin{:});
[f2,zc] = f(x2,varargin{:});
if(~isnan(f1)||isnan(f2))
    error('Wrong upper/lower bound for bi-section search.');
end
[fm,zcm] = f(xm,varargin{:},zc);
for i = 1:100
   if (x2-x1) < err
       x_max = xm;
       break;
   end
   if isnan(fm)
        xm2 = d(d(x1,xm),xm);
        [fm2,zcm2] = f(xm2,varargin{:},zc);
        if isnan(fm2)
            x1 = xm;
        else
            x2 = xm2;
            zc = zcm2;
        end
   else
       x2 = xm;
       zc = zcm;
   end
   xm = d(x1,x2);
   [fm,zcm] = f(xm,varargin{:},zc);
end
end