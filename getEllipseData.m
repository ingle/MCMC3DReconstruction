% [XDAT, YDAT, FDAT] = GETELLIPSEDATA(A,B,limX, limY,N)
% xdat: xcoor of data points
% ydat: ycoor of data points
% fdat: ellipsoid data f(xdat, ydat)
%
% limX, limY are box limits along x,y
% a,b are the three semiprincipal axis lengths along x,y
% N number of data points to produce
%
% fdat = 1 if xdat.^2/a^2 + ydat.^2/b^2 <= 1
% else fdat = 0.
%
% v1 Jan 10, 2015
%
function [xdat, ydat, fdat] = getEllipseData(a, b, limX, limY, N)
assert( a<limX && b<limY );

xdat = unifrand(-limX, limX, N);
ydat = unifrand(-limY, limY, N);

fdat = zeros(size(xdat));

fdat( xdat.^2/a^2 + ydat.^2/b^2 <= 1 ) = 1;
end

% generate n unifomly distributed samples in [a,b]
function r = unifrand(a,b,n)
rng('default');
assert( a<b );
r = a + (b-a)*rand(n,1);
end
