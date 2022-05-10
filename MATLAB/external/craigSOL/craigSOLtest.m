function x = craigSOLtest( m, n )

%        x = craigSOLtest(  m, n );
%        x = craigSOLtest( 10,10 );
%        x = craigSOLtest( 10,20 );
%
% 15 Aug 2014: craigSOLtest.m derived from lsqrSOLtest.m.

  A      = @(v,mode) Aprodxxx( v,mode,m,n );  % Nested function

  xtrue  = (n : -1: 1)';
  b      = A(xtrue,1);

  atol   = 1.0e-6;
  btol   = 1.0e-6;
  conlim = 1.0e+10;
  itnlim = 10*n;
  show   = 1;

  [ x, y, istop, itn, rnorm, anorm, acond, xnorm]...
      = craigSOL( m, n, A, b, atol, btol, conlim, itnlim, show );  

  disp(' ');   j1 = min(n,5);   j2 = max(n-4,1);
  disp('First elements of x:');  disp(x(1:j1)');
  disp('Last  elements of x:');  disp(x(j2:n)');

  r    = b - A(x,1);
  r1   = norm(r);
  disp(' ')
  str1 = sprintf( 'rnorm (est.)  %10.3e %10.3e', rnorm );
  str2 = sprintf( 'rnorm (true)  %10.3e %10.3e', r1    );
  disp(str1)
  disp(str2)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Nested functions (only 1 here).
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function y = Aprodxxx( x, mode, m, n )

    % Private function.
    % if mode = 1, computes y = A*x
    % if mode = 2, computes y = A'*x
    % for some matrix  A.
    %
    % This is a simple example for testing  LSQR.
    % It uses the leading m*n submatrix from
    % A = [ 1
    %       1 2
    %         2 3
    %           3 4
    %             ...
    %               n ]
    % suitably padded by zeros.
    %
    % 11 Apr 1996: First version for distribution with lsqr.m.
    %              Michael Saunders, Dept of EESOR, Stanford University.
    % 15 Aug 2014: Borrowed for testing craigSOL.m.

    if mode == 1
      d  = (1:n)';  % Column vector
      y1 = [d.*x; 0] + [0;d.*x];
      if m <= n+1
	y = y1(1:m);
      else         
	y = [     y1; 
	  zeros(m-n-1,1)];
      end
    else
      d  = (1:m)';  % Column vector
      y1 = [d.*x] + [d(1:m-1).*x(2:m); 0];
      if m >= n
	y = y1(1:n);
      else
	y = [y1;
	  zeros(n-m,1)];
      end
    end

  end % nested function Aprodxxx

end % function craigSOLtest

