function [ x, y, istop, itn, rnorm, Anorm, Acond, xnorm ]...
  = craigSOL( m, n, A, b, atol, btol, conlim, itnlim, show )
%
%        [ x, y, istop, itn, rnorm, Anorm, Acond, xnorm ]...
% = craigSOL( m, n, A, b, atol, btol, conlim, itnlim, show );
%
% CRAIG finds a solution x to the linear equation Ax = b, where
% A is a real matrix with m rows and n columns, and b is a real
% m-vector.  If A is square and nonsingular, CRAIG finds the
% unique solution x = A(inv)b.
% If the system Ax = b is under-determined (i.e. there are many
% solutions), CRAIG finds the solution of minimum Euclidean length,
% namely    x = A inv(A A') b.  Thus, CRAIG solves the problem
%
%          min  x'x  subject to  Ax = b.
%
% y returns a vector satisfying A'y = x.  Hence AA'y = b.
%
% A is an m by n matrix (ideally sparse),
% or a function handle such that
%    y = A(x,1) returns y = A*x   (where x will be an n-vector);
%    y = A(x,2) returns y = A'*x  (where x will be an m-vector).

%-----------------------------------------------------------------------
% CRAIG uses an iterative (conjugate-gradient-like) method.
% For further information, see 
% 1. C. C. Paige and M. A. Saunders (1982a).
%    LSQR: An algorithm for sparse linear equations and sparse least squares,
%    ACM TOMS 8(1), 43-71.

% 08 Apr 2003: First craig.m derived from Fortran 77 version of craig1.for.
%              Michael Saunders, Systems Optimization Laboratory,
%              Dept of MS&E, Stanford University.
% 09 Apr 2003: Experimenting on singular systems (for inverse iteration).
%              Separately, on full-rank Ax = b, "Acond" seems to
%              over-estimate cond(A) drastically.
% 02 Oct 2006: Output y such that x = A'y (already in the f77 version).
% 15 Aug 2014: A can now be a matrix or a function handle, as in lsqrSOL.
% 28 Aug 2014: Fixed glitches found by Dominique Orban.
%-----------------------------------------------------------------------

% Initialize.

if isa(A,'numeric')
  explicitA = true;
elseif isa(A,'function_handle')
  explicitA = false;
else
  error('SOL:craigSOL:Atype','%s','A must be numeric or a function handle');
end

msg=['The exact solution is  x = 0                              '
     'Ax - b is small enough, given atol, btol                  '
     'The system Ax = b seems to be incompatible                '
     'The estimate of cond(A) has exceeded conlim               '
     'Ax - b is small enough for this machine                   '
     'The system Ax = b seems to be incompatible                '
     'Cond(A) seems to be too large for this machine            '
     'The iteration limit has been reached                      '];

if show
   disp(' ')
   disp('CRAIG          minimum-length solution of  Ax = b')
   str1 = sprintf('The matrix A has %8g rows  and %8g cols', m, n);
   str3 = sprintf('atol = %8.2e                 conlim = %8.2e', atol, conlim);
   str4 = sprintf('btol = %8.2e                 itnlim = %8g'  , btol, itnlim);
   disp(str1);   disp(str3);   disp(str4);
end

itn    = 0;
istop  = 0;
ctol   = 0;
if conlim > 0, ctol = 1/conlim; end
Anorm  = 0;
Acond  = 0;
xnorm  = 0;

%     Set beta(1) and u(1) for the bidiagonalization.
%         beta*u = b.

v      = zeros(n,1);
x      = zeros(n,1);
w      = zeros(m,1);
y      = zeros(m,1);

beta   = norm(b);
bnorm  = beta;
rnorm  = beta;
if beta==0, disp(msg(istop+1,:)); return, end
u      = (1/beta)*b;

% More initialization.
% aanorm  is norm(L_k)**2, an estimate of norm(A)**2.  It is
%                alfa1**2  +  (alfa2**2 + beta2**2)  +  ...
% ddnorm  is norm(D_k)**2, an estimate of norm( (A'A)inverse ).
% xxnorm  is norm(x_k)**2  =  norm(z_k)**2.

aanorm = 0;
ddnorm = 0;
xxnorm = 0;
alfa   = 1;
z      = -1;

if show
   disp(' ')
   str1 = '   Itn      x(1)       rnorm      xnorm';
   str2 = '     Norm A   Cond A   alpha    beta';
   disp([str1 str2])
   str1   = sprintf( '%6g %12.5e',        itn, x(1)  );
   str2   = sprintf( ' %10.3e',  rnorm );
   disp([str1 str2])
end

%------------------------------------------------------------------
%     Main iteration loop.
%------------------------------------------------------------------
while itn < itnlim
  itn    = itn + 1;
  oldalf = alfa;

  % Perform the next step of the bidiagonalization to obtain the
  % next alfa, v, beta, u.  These satisfy the relations
  %      alfa*v  =  A'*u  -  beta*v.
  %      beta*u  =  A*v   -  alfa*u,

  if explicitA
      v = A'*u   - beta*v;
  else
      v = A(u,2) - beta*v;
  end
  alfa   = norm(v);
  if alfa==0, istop = 2; disp(msg(istop+1,:)); return, end
  v      = (1/alfa)*v;

  aanorm =   aanorm  +  alfa^2;
  z      = - (beta/alfa)*z;
  x      =   x + z*v;
  t1     = - beta/oldalf;
  t2     =   z   /alfa;
  t3     =   1   /alfa;

  w      =   u + t1*w;
  y      =   y + t2*w;
  ddnorm =   ddnorm  + norm(t3*w)^2;

  if explicitA
      u = A*v    - alfa*u;
  else 
      u = A(v,1) - alfa*u; 
  end
  beta = norm(u);
  if beta > 0
      u = (1/beta)*u;
  end

  %===============================================================
  % Test for convergence.
  % We estimate various norms and then see if
  % the quantities test1 or test3 are suitably small.
  %===============================================================
  aanorm =   aanorm + beta^2;
  Anorm  =   sqrt( aanorm );
  Acond  =   sqrt( ddnorm )*Anorm;
  xxnorm =   xxnorm + z^2;
  rnorm  =   abs( beta*z );
  xnorm  =   sqrt( xxnorm );

  test1  =   rnorm/bnorm;
  test3  =   1    /Acond;
  t1     =   test1/(1 + Anorm*xnorm/bnorm);
  rtol   =   btol + atol*Anorm*xnorm/bnorm;

  % The following tests guard against extremely small values of
  % atol, btol  or  ctol.  (The user may have set any or all of
  % the parameters  atol, btol, conlim  to zero.)
  % The effect is equivalent to the normal tests using
  % atol = eps,  btol = eps,  conlim = 1/eps.

  if itn >= itnlim, istop = 7; end
  if test3 <= eps , istop = 6; end
  if t1    <= eps , istop = 4; end
  % Allow for tolerances set by the user.
  if test3 <= ctol, istop = 3; end
  if test1 <= rtol, istop = 1; end

  % See if it is time to print something.

  prnt = 0;
  if n     <= 40       , prnt = 1; end
  if itn   <= 10       , prnt = 1; end
  if itn   >= itnlim-10, prnt = 1; end
  if rem(itn,10) == 0  , prnt = 1; end
  if test3 <=  2*ctol  , prnt = 1; end
  if test1 <= 10*rtol  , prnt = 1; end
  if istop ~=  0       , prnt = 1; end

  if prnt == 1
      if show
          str1 = sprintf( '%6g %12.5e',       itn, x(1 ) );
          str2 = sprintf( ' %10.3e %10.3e', rnorm, xnorm );
          str4 = sprintf( ' %8.1e %8.1e',   Anorm, Acond );
          str5 = sprintf( ' %8.1e %8.1e',   alfa , beta  );
          disp([str1 str2 str4 str5])
      end
  end
  if istop > 0, break, end
end

if show
   disp(' ')
   disp('CRAIG finished')
   disp(msg(istop+1,:))
   disp(' ')
   str1 = sprintf( 'istop =%8g    rnorm =%8.1e',   istop, rnorm );
   str2 = sprintf( 'Anorm =%8.1e',  Anorm );
   str3 = sprintf( 'itn   =%8g    xnorm =%8.1e',     itn, xnorm );
   str4 = sprintf( 'Acond =%8.1e',  Acond );
   disp([str1 '   ' str2])
   disp([str3 '   ' str4])
   disp(' ')
end % function craigSOL
