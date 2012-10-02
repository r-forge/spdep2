function llike = f_sar(rho,detval,epe0,eped,epe0d,n)
% PURPOSE: evaluates concentrated log-likelihood for the 
%  spatial autoregressive model using sparse matrix algorithms
% ---------------------------------------------------
%  USAGE:llike = f_sar(rho,detval,epe0,eped,epe0d,n)
%  where: rho  = spatial autoregressive parameter
%         detval = a matrix with vectorized log-determinant information
%                or a size [1 1] fitted spline model /* RSB */ 
%         epe0   = see below
%         eped   = see below
%         eoe0d  = see below
%         n      = # of obs
%          b0 = AI*xs'*ys;
%          bd = AI*xs'*Wys;
%          e0 = ys - xs*b0;
%          ed = Wys - xs*bd;
%          epe0 = e0'*e0;
%          eped = ed'*ed;
%          epe0d = ed'*e0;
% ---------------------------------------------------
%  RETURNS: a  scalar equal to minus the log-likelihood
%           function value at the parameter rho
% ---------------------------------------------------                         
%  NOTE: this is really two functions depending
%        on nargin = 3 or nargin = 4 (see the function)
%  --------------------------------------------------
%  SEE ALSO: sar, f_far, f_sac, f_sem
% ---------------------------------------------------

% written by: James P. LeSage 1/2000
% University of Toledo
% Department of Economics
% Toledo, OH 43606
% jlesage@spatial-econometrics.com
global funvals;
if nargin == 6
  if  all(size(detval) == [1 1]) % /* RSB */
    detm = ppval(detval, rho); % /* RSB */
  elseif all(size(detval) == [n 1]) % /* RSB */
    detm = sum(log(ones(n,1) - rho .* detval)); % /* RSB */
  else % /* RSB */
    gsize = detval(2,1) - detval(1,1);
% Note these are actually log detvalues
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
%  fprintf(1, 'in f_sar.m/f_sar\n') % /* RSB */

  z = epe0 - 2*rho*epe0d + rho*rho*eped;

  llike = (n/2)*log(z) - detm;
  thisi = [rho detm llike]; % /* RSB */
  funvals = [funvals; thisi]; % /* RSB */
else

  error('f_sar: Wrong # of input arguments');

end;
