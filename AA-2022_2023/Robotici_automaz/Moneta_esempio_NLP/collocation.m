function [D, C, B] = collocation(degree)
    import casadi.*

    % Degree of interpolating polynomial
    d = degree;

    % Get collocation points
    tau_root = [0 collocation_points(d, 'radau')];

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
      for r=1:d+1 %costruisco le basi considerando come punti quelli di collocazione e il punto iniziale
        if r ~= j %calcolo il prodotto tra le basi di lagrange (dove ho normalizzato per t)
          coeff = conv(coeff, [1, -tau_root(r)]); %prodotto tra i coefficienti dei polinomi
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
end

%N.B. il polinomio viene costruito facendolo passare per il punto di inizio
%intervallo, questa cosa va bene per calcolare D e C ma non per B
%per sfruttare la precisione data dal metodo di gauss-legendre i polinomi
%devono essere fatti passare dai punti di collocazione (radici) e basta. Su
%tali polinomio posso poi calcolare gli integrali e valutarli nei punti di
%collazione. Tali integrali saranno i pesi delle formule di quadratura