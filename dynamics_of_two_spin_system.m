% H = S_1^z + S_2^z + 2\lambda \vec{S_1} \cdot \vec{S_2}

%parameters
lambda = 0.5;
t_end = 1e2;
line = '-';
T = 1;
beta = 1/T;
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

%H in bases of |S_1^z, S_2^z >
H = [lambda/2-1 , 0, 0, 0; 0, -lambda/2, lambda, 0; 0, lambda, -lambda/2, 0; 0, 0, 0, lambda/2+1];
H
[eigen_vecs, eigen_values] = eig(H); 
disp("eigen vectors") 
eigen_vecs %bases |S, S^z> where S = S_1 + S_2 : |1,-1>, |0,0>, |1,-1> |1,1>
disp("eigen values")
eigen_values

%initial density matrix
rho0_matrix_eigen_1 = [ 0.4,   0,   0, 0.1];
rho0_matrix_eigen_2 = [   0, 0.3,   0,   0];
rho0_matrix_eigen_3 = [   0,   0, 0.2,   0];
rho0_matrix_eigen_4 = [ 0.1,   0,   0, 0.1];
rho0_matrix_eigen = [rho0_matrix_eigen_1; rho0_matrix_eigen_2; rho0_matrix_eigen_3; rho0_matrix_eigen_4];
rho0_matrix_eigen
rho0_matrix = eigen_vecs * rho0_matrix_eigen * eigen_vecs'
rho0_vec = reshape(rho0_matrix, 1, size(rho0_matrix,1)*size(rho0_matrix,2));

tspan = [0, t_end];
ode_func = @(t, rho_vec) drhodt_twospins(t, rho_vec, H);
%[t, rho_vecs] = ode45(ode_func, tspan, rho0_vec, options);
[t, rho_vecs] = ode45(ode_func, tspan, rho0_vec);
rho_matrices = reshape(rho_vecs.', size(H,1), size(H,2), []);
total_iter = size(rho_matrices, 3);
rho_matrices_eigen = zeros(4,4,total_iter);

entropy_arr = zeros(1, total_iter);
rho_trace_arr = zeros(1, total_iter);
for i=1:total_iter
   entropy_arr(i) = entropy(rho_matrices(:,:,i));
   rho_trace_arr(i) = trace(rho_matrices(:,:,i));
   rho_matrices_eigen(:,:,i) = eigen_vecs' * rho_matrices(:,:,i) * eigen_vecs;
end

figure(1)
plot(t, entropy_arr, ['r',line])
title("entropy vs time")



figure(2)
plot(t, reshape(rho_matrices(1,1,:), 1, total_iter), ['r',line])
hold on
plot(t, reshape(rho_matrices(2,2,:), 1, total_iter), ['b',line])
hold on
plot(t, reshape(rho_matrices(3,3,:), 1, total_iter), ['b',line])
hold on
plot(t, reshape(rho_matrices(4,4,:), 1, total_iter), ['r',line])
hold on
plot(t, real(reshape(rho_matrices(2,3,:), 1, total_iter)), ['g',line])
hold on
plot(t, imag(reshape(rho_matrices(2,3,:), 1, total_iter)), ['g',line])
legend('\rho_{11}', '\rho_{22}', '\rho_{33}', '\rho_{44}', 'Re(\rho_{23})', 'Im(\rho_{23})')
title("matrix elements vs time(|S_1^z, S_2^z>)")

figure(3)
plot(t, reshape(rho_matrices_eigen(1,1,:), 1, total_iter), ['r',line])
hold on
plot(t, reshape(rho_matrices_eigen(2,2,:), 1, total_iter), ['b',line])
hold on
plot(t, reshape(rho_matrices_eigen(3,3,:), 1, total_iter), ['b',line])
hold on
plot(t, reshape(rho_matrices_eigen(4,4,:), 1, total_iter), ['r',line])
hold on
plot(t, real(reshape(rho_matrices_eigen(2,3,:), 1, total_iter)), ['g',line])
hold on
plot(t, imag(reshape(rho_matrices_eigen(2,3,:), 1, total_iter)), ['g',line])
legend('\rho_{11}', '\rho_{22}', '\rho_{33}', '\rho_{44}', 'Re(\rho_{23})', 'Im(\rho_{23})')
title("matrix elements vs time(|S, S^z>)")

figure(4)
plot(t, rho_trace_arr)

