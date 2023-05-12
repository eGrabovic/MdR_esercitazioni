omega_n = linspace(2, 500, 2000);
e_dot_0 = max(full(q_d_des_fun(0)));
e_0 = 0;
e_max = @(x) exp(-e_dot_0./(e_0.*x + e_dot_0)).*(e_0 + e_dot_0./x);

figure; hold on; grid on;
plot(omega_n, e_max(omega_n).*180/pi)
xlabel('$\omega_n$')
ylabel('$e_\textrm{max}$ (deg)')
set(gca, 'fontsize', 26)