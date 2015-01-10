% [XDAT, YDAT, ZDAT, FDAT] = GETELLIPSOIDDATA(A,B,C,limX, limY, limZ,N)
% xdat: xcoor of data points
% ydat: ycoor of data points
% zdat: zcoor of data points
% fdat: ellipsoid data f(xdat, ydat, zdat)
%
% limX, limY, limZ are box limits along x,y,z
% a,b,c are the three semiprincipal axis lengths along x,y,z
% N number of data points to produce
%
% fdat = 1 if xdat.^2/a^2 + ydat.^2/b^2 + zdat.^/c^2 <=1
% else fdat = 0.
%
% v1 Jan 10, 2015
%
function [xdat, ydat, zdat, fdat] = getEllipsoidData(a, b, c, limX, limY, limZ, N)
assert( a<limX && b<limY && c<limZ );

xdat = unifrand(-limX, limX, N);
ydat = unifrand(-limY, limY, N);
zdat = unifrand(-limZ, limZ, N);

fdat = zeros(size(xdat));

fdat( xdat.^2/a^2 + ydat.^2/b^2 + zdat.^2/c^2 <= 1 ) = 1;
end

% generate n unifomly distributed samples in [a,b]
function r = unifrand(a,b,n)
rng('default');
assert( a<b );
r = a + (b-a)*rand(n,1);
end
