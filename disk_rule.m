function [ w, x, y, r, theta ] = disk_rule ( nr, nt, xc, yc, rc )

%*****************************************************************************80
%
%% DISK_RULE computes a quadrature rule for a general disk.
%
%  Discussion:
%
%    The general disk is the region:
%
%      ( x - xc ) ^ 2 + ( y - yc ) ^ 2 <= rc ^ 2.
%
%    The integral I(f) is then approximated by
%
%      S(f) = sum ( 1 <= i <= NT * NR ) W(i) * F ( X(i), Y(i) ).
%
%      Area = pi * RC ^ 2
%
%      Q(f) = Area * S(f)
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    18 April 2016
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer NR, the number of points in the radial rule.
%
%    Input, integer NT, the number of angles to use.
%
%    Input, real XC, YC, the coordinates of the disk center.
%
%    Input, real RC, the radius of the disk.
%
%    Output, real W(NR*NT), the weights for the rule.
%
%    Output, real X(NR*NT), Y(NR*NT), the points for the rule.
%
  [ w01, r01, t01 ] = disk01_rule ( nr, nt );
%
%  Recompute the rule for the general circle in terms of X, Y.
%
  w = repmat ( w01, 1, nt );
  x = xc + rc * r01(1:nr) * cos ( t01(1:nt) )';
  y = yc + rc * r01(1:nr) * sin ( t01(1:nt) )';
%
%  Return column vectors.
%
  w = reshape ( w, nr * nt, 1 );
  x = reshape ( x, nr * nt, 1 );
  y = reshape ( y, nr * nt, 1 );
  

  r = meshgrid(rc * r01(1:nr),1:nt)';
  r = reshape ( r, nr * nt, 1 );
  theta = meshgrid(t01(1:nt),1:nr);
  theta = reshape ( theta, nr * nt, 1 );

  return
end
