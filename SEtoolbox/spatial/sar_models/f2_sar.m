function llike = f2_sar(parm,y,x,W,detval)
% PURPOSE: evaluates log-likelihood -- given ML estimates
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:llike = f2_sar(parm,y,X,W,ldet)
%  where: parm = vector of maximum likelihood parameters
%                parm(1:k-2,1) = b, parm(k-1,1) = rho, parm(k,1) = sige
%         y    = dependent variable vector (n x 1)
%         X    = explanatory variables matrix (n x k)
%         W    = spatial weight matrix
%         ldet = matrix with [rho log determinant] values
%                computed in sar.m using one of Kelley Pace's routines  
%                or a size [1 1] fitted spline model /* RSB */ 
% ---------------------------------------------------
%  RETURNS: a  scalar equal to minus the log-likelihood
%           function value at the ML parameters
%  --------------------------------------------------
%  NOTE: this is really two functions depending
%        on nargin = 4 or nargin = 5 (see the function)
% ---------------------------------------------------
%  SEE ALSO: sar, f2_far, f2_sac, f2_sem
% ---------------------------------------------------

% written by: James P. LeSage 1/2000
% University of Toledo
% Department of Economics
% Toledo, OH 43606
% jlesage@spatial.econometrics.com

n = length(y); 
k = length(parm);
b = parm(1:k-2,1);
rho = parm(k-1,1);
sige = parm(k,1);

if  all(size(detval) == [1 1]) % /* RSB */
  detm = ppval(detval, rho); % /* RSB */
elseif all(size(detval) == [n 1]) % /* RSB */
  detm = sum(log(ones(n,1) - rho .* detval)); % /* RSB */
else % /* RSB */
  gsize = detval(2,1) - detval(1,1);
  i1 = find(detval(:,1) <= rho + gsize);
  i2 = find(detval(:,1) <= rho - gsize);
  i1 = max(i1);
  i2 = max(i2);
  index = round((i1+i2)/2);
  if isempty(index)
    index = 1;
  end;
  detm = detval(index,2);
end; % /* RSB */
%fprintf(1, 'in f2_sar.m/f2_sar\n') % /* RSB */
e = y-x*b-rho*sparse(W)*y;
epe = e'*e;
tmp2 = 1/(2*sige);
llike = -(n/2)*log(2*pi) - (n/2)*log(sige) + detm - tmp2*epe; % /* RSB */

