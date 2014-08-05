function cost1 = funcostrarc(p,fitfun,DST,PST,TST)
%function cost1 = funcostrarc(p,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1,ifast,partials)
%function cost1 = funcostrarc(p,fitfun,varargin)
%function cost1 = funcostrarc(p,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1)
% cost  function for phase model
%   p   == parameter
%   xyzm == distance in m  
%   xd == observed (wrapped) phase in radians [-pi, +pi]
%   yd == dummy, for historical reasons
% 2010-JAN-11
%
% for use with ANNEAL

nargchk(5, 5, nargin);

if numel(p) ~= numel(PST.p0)
    error(sprintf('Dimension mismatch %d %d\n',numel(p),numel(PST.p0)));
end
PST.p1 = p;
% field of costs in DN [0,127]
%costs = funcostsiarc1(p,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1);
%costs = funcostsiarc1(p,fitfun,varargin{:});
%costs = funcostsiarc1(p,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1,ifast,partials);
% field of costs on DN [0, 
%costs = funcostsrarc(p,fitfun,xyzm,tepochs,bpest,dops,DD,unitv,xd,yd,ippix1,ifast,partials);
costs = funcostsrarc(fitfun,DST,PST,TST);
%size(costs)

% number of elements
n=numel(costs);

% average cost is L1 norm in radians
cost1=sum(colvec(costs))/n; 
% convert to cycles
cost1 = cost1/2.0/pi;

return;

