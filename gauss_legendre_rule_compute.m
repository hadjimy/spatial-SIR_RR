function [ x, w ] = gauss_legendre_rule_compute (n, a, b )

%*****************************************************************************80
%
%% GAUSS_LEGENDRE_RULE_COMPUTE computes a Gauss-Legendre rule.
%
%  Discussion:
%
%    This routine should be very efficient in MATLAB.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    03 March 2007
%
%  Author:
%
%    Lloyd Trefethen
%
%  Reference:
%
%    Lloyd Trefethen,
%    Is Gauss Quadrature Better than Clenshaw-Curtis?,
%    SIAM Review,
%    Volume 50, Number 1, 2008, pages 67-87.
%
%  Parameters:
%
%    Input, integer N, the order of the rule.
%
%    Output, real X(N), the abscissas of the rule.
%
%    Output, real W(N), the weights of the rule.
%
  beta = 0.5 ./ sqrt ( 1 - ( 2 * ( 1:n-1 ) ) .^ (-2) );

  tri = diag ( beta, 1 ) + diag ( beta, -1 );

  [ eigval, eigvec ] = eig ( tri );

  x = diag ( eigvec );
  [ x, i ] = sort ( x );
  x = x';

  w = 2 * eigval(1,i) .^ 2;

  % Shift abscissas and weights from [-1,1] interval to [a,b]
  x = ( b - a ) / 2 .* x + ( a + b ) / 2;
  w = ( b - a ) / 2 .* w;
  
  return
end
