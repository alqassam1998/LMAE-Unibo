function x = solveLinearSystem( A , t, x0)

% Simple function solving equation dxdt = Ax, with x(0) = x0 in the
% time-span t.
% there's a banch of smarter ways to solve these equations, this is the
% simplest (and probably less efficient) way, however this is the one you
% studied.
% a good exercise could be to create a function solving these equations
% using the Jordan form or the diagonal form

    x = zeros(length(x0),length(t));
    for i=1:length(t)
        x(:,i) = expm(A*t(i))*x0;
    end

end

