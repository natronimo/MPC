function [u, du] = MPC(x, u, r, A, B, C, D, Q, R, S, Ts, N_p, N_c, y_max, y_min, u_max, u_min, du_max, du_min)

n = length(x);
m = length(u);
p = length(r);

x = [x; u];
r = repmat([r; u], N_p, 1);

A = [A, B; zeros(m, n), eye(m)];
B = [B; eye(m)];
C = [C, D; zeros(m, n), eye(m)];
D = [D; eye(m)];

n = n + m;
p = p + m;

du_max = du_max*Ts;
du_min = du_min*Ts;

Q_bar = diag(repmat([Q, R], 1, N_p));

R_bar = diag(repmat(S, 1, N_c));

S_bar = zeros(p*N_p, m*N_c);
for i = 1:N_p
	for j = 1:min(i, N_c)
		S_bar(i*p-p+1:i*p, j*m-m+1:j*m) = C*A^(i-j)*B;
	end
	if i < N_c
		S_bar(i*p-p+1:i*p, i*m+1:i*m+m) = D;
	end
end

T_bar = zeros(p*N_p, n);
for i = 1:N_p
	T_bar(i*p-p+1:i*p, :) = C*A^i;
end

H = 2*(R_bar + S_bar'*Q_bar*S_bar);
f = 2*(x'*T_bar' - r')*Q_bar*S_bar;

Aineq_x = [S_bar; -S_bar];
bineq_x = [repmat([y_max; u_max], N_p, 1); repmat([-y_min; -u_min], N_p, 1)] - [T_bar; -T_bar]*x;
Aineq_du = [eye(m*N_c); -eye(m*N_c)];
bineq_du = [repmat(du_max, N_c, 1); repmat(-du_min, N_c, 1)];
Aineq = [Aineq_x; Aineq_du];
bineq = [bineq_x; bineq_du];

Aeq = zeros(0, m*N_c);
beq = zeros(0, 1);

iA0 = false(2*(p*N_p + m*N_c), 1);

opt = mpcActiveSetOptions;

du = mpcActiveSetSolver(H, f', Aineq, bineq, Aeq, beq, iA0, opt);
du = du(1:m);
u = u + du;
