%H = S_z + \theta(t)*\lambda*S_x 

%parameters
lambda = 0.1;
t_end = 1e1;
line = '--';

%in the bases of eigenstates of [1, lambda; lambda, -1]
H = [1, lambda; lambda, -1];
[eigen_vecs, eigen_values] = eig(H);

offdiag_matrix = [-lambda, lambda; lambda, -lambda];
A = zeros(4,4);
A(1:2, 3:4) = offdiag_matrix;
A(3:4, 1:2) = offdiag_matrix;
A(3:4, 3:4) = [2, 0; 0, -2];
disp(A);

%initial rho
%rho0 = [0.5, 0.5, 0.5, 0.5]; %rho_11 rho_22 rho_12 rho_21
rho0_matrix_eigen = [0.3, 0.1; 0.1, 0.7];
%rho0_matrix_eigen = [0, 0; 0, 1];
rho0_matrix = eigen_vecs * rho0_matrix_eigen * eigen_vecs';
rho0 = [rho0_matrix(1,1), rho0_matrix(2,2), rho0_matrix(1,2), rho0_matrix(2,1)]; 

tspan = [0, t_end];
[t, rho] = ode45(@(t, y) drhodt(t, y, A), tspan, rho0);
figure(1)
plot(t, rho(:,1), ['b',line])
hold on
plot(t, rho(:,2), ['r', line])
hold on
plot(t, real(rho(:,3)), ['g', line])
hold on
plot(t, imag(rho(:,3)), ['g', line])
ylim([-0.1,1.1])
legend('\rho_{11}', '\rho_{22}', 'Re(\rho_{12})', 'Im(\rho_{12})')
title("Matrix elements vs time")

entropy_arr = zeros(1, size(t, 1));
rho_eigen = zeros(size(t,1), 4);
for i=1:size(t,1)
    rho_matrix = [rho(i,1), rho(i,3); rho(i,4), rho(i,2)];
    rho_matrix_eigen = eigen_vecs' * rho_matrix * eigen_vecs;
    rho_eigen(i,:) = [rho_matrix_eigen(1,1), rho_matrix_eigen(2,2), rho_matrix_eigen(1,2), rho_matrix_eigen(2,1)];
    entropy_arr(i) = entropy(rho_matrix);
end
figure(2)
plot(t, entropy_arr)
title("entropy vs time")

figure(3)
plot(t, rho_eigen(:,1), ['b', line])
hold on
plot(t, rho_eigen(:,2), ['r', line])
hold on 
plot(t, real(rho_eigen(:,3)), ['g', line])
hold on
plot(t, imag(rho_eigen(:,3)), ['g', line])
ylim([-0.1, 1.1])
legend('\rho_{11}', '\rho_{22}', 'Re(\rho_{12})', 'Im(\rho_{12})')
title("Matrix elements(eigen) vs time")



%thermal density matrix

