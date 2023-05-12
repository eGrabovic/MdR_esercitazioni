function [x_next, S, x_colloc] = Direct_collocation_num(x, u, xdot, dt, options)


arguments
   x
   u
   xdot
   dt
   options.degree = 2;
   options.Solver = [];
   options.x_colloc_guess = [];
end

h = dt;
out = x;
S = options.Solver;
d = options.degree;
pp = full(xdot(x,u));
gain = linspace(1,d,d);


%% Matrix for direct collocation
    import casadi.*

    % Degree of interpolating polynomial
    d = options.degree;

    % Get collocation points
    tau_root = [0 collocation_points(d, 'legendre')];

    % Coefficients of the collocation equation
    C = zeros(d+1,d+1);

    % Coefficients of the continuity equation
    D = zeros(d+1, 1);

    % Coefficients of the quadrature function
    B = zeros(d+1, 1);

    % Construct polynomial basis
    for j=1:d+1
      % Construct Lagrange polynomials to get the polynomial basis at the collocation point
      coeff = 1;
      for r=1:d+1
        if r ~= j
          coeff = conv(coeff, [1, -tau_root(r)]);
          coeff = coeff / (tau_root(j)-tau_root(r));
        end
      end
      % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      D(j) = polyval(coeff, 1.0);

      % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      pder = polyder(coeff);
      for r=1:d+1
        C(j,r) = polyval(pder, tau_root(r));
      end

      % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
      pint = polyint(coeff);
      B(j) = polyval(pint, 1.0);
    end



    function F = eq(X,out2,d2,C2,xdot2,control2,h2)
        
        F = casadi.MX(length(X), d2);
        % Loop over collocation points
        for jj=1:d2
           % Expression for the state derivative at the collocation point
           xp = C2(1,jj+1)*out2;
           for rr=1:d2
               xp = xp + C2(rr+1,jj+1)*X(:,rr);
           end

           % Append collocation equations
           fj = xdot2(X(:,jj),control2);
           F(:,jj) = h2*fj - xp;

           % Add contribution to the end state
        end
    end


    
    if isempty(S) % build a casadi findroot solver once
        X = MX.sym('x2', length(x), d); % x_k+1
        x_1 = MX.sym('x2', length(x), 1);
        fun = eq(X,x_1,d,C,xdot, u,h);
        problem = struct('x', X(:), 'p', x_1, 'g', fun(:));
        S = rootfinder('S', 'newton', problem);
        
    end

%     options_fun = optimoptions('fsolve','Display','none');
%     [x_colloc, ~, exF] = fsolve(fun, x+pp*gain*(h/d), options_fun);% [x+pp*(h/3),x+pp*(2*h/3),x+(pp*h)]);
    if isempty(options.x_colloc_guess)
        x0 = out+pp*gain*(h/d);
    else
        x0 = options.x_colloc_guess;
    end
    tic
    sol = S('x0', x0(:), 'p', out);
    toc
    x_colloc = reshape(full(sol.x), length(x), d);
                        
    x_next = D(1)*x;
    x_next = x_next + sum(D(2:end)'.*x_colloc,2);


end