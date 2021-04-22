function PXPmodel2

%%parameters
global imagunit tolerance 
imagunit = sqrt(-1);
tolerance = 1.0e-12;
total_iter = 5000;
dt = 0.1;
tem = 1;
beta = 1/tem;
N_sites = 10; %num of sites in subsystem A
N_sites_B = 6; %num of sites in subsystem B

%%two seperate problems
%\rho_{11}, first spin of B is down, so the last spin of A can be either up
%or down
N_sites_1 = N_sites;
bases_1 = Bases(N_sites_1);
bases_constrained_1 = Bases_constrained(bases_1, N_sites_1)
%disp("Size of bases_1:")
%size(bases_constrained_1, 1)
%\rho_{22}, first spin of B is up, so the last spin of A must be down
N_sites_2 = N_sites - 1;
bases_2 = Bases(N_sites_2);
bases_constrained_2 = Bases_constrained(bases_2, N_sites_2);
bases_constrained_2_full = [bases_constrained_2, -ones(size(bases_constrained_2, 1), 1)]
%disp("Size of bases_2:")
%size(bases_constrained_2, 1)

%%Hamiltonian A of two seperate problems
%H_A for \rho_{11}
H_a_1 = Hamiltonian_a(bases_constrained_1);
%H_A for \rho_{22}
H_a_2 = Hamiltonian_a(bases_constrained_2);

%%eigenstates and eigenvalues of H_a_1 and H_a_2
[eigenstates_1,eigen_energies_mat_1] = eig(H_a_1);
%eigenstates_1
%eigen_energies_1 = diag(eigen_energies_mat_1)
[eigenstates_2,eigen_energies_mat_2] = eig(H_a_2);
%eigenstates_2
%eigen_energies_2 = diag(eigen_energies_mat_2)

%%subsystem B
bases_B = Bases(N_sites_B);
bases_B_constrained = Bases_constrained(bases_B, N_sites_B);
%disp("Size of bases_B:")
%size(bases_B_constrained, 1)
[bases_B_1, bases_B_2]=split_bases_B(bases_B_constrained);
bases_B_1
%disp("Size of bases_B_1:")
%size(bases_B_1, 1)
bases_B_2
%disp("Size of bases_B_2:")
%size(bases_B_2, 1)XLp1PLp2(bases_B_2, bases_B_1)

%%thermal density matrix of subsystem B
H_b_1 = Hamiltonian_b(bases_B_1);
H_b_2 = Hamiltonian_b(bases_B_2);
[eigenstates_b_1, eigen_energies_mat_b_1] = eig(H_b_1);
[eigenstates_b_2, eigen_energies_mat_b_2] = eig(H_b_2);

frequency_b_1 = diag(eigen_energies_mat_b_1)./(2*pi)
frequency_b_2 = diag(eigen_energies_mat_b_2)./(2*pi)

rho11_b_therm = zeros(size(eigen_energies_mat_b_1));
for i=1:size(rho11_b_therm,1)
    rho11_b_therm(i,i) = exp(-beta*eigen_energies_mat_b_1(i,i));
end
rho11_b_therm = rho11_b_therm./trace(rho11_b_therm);
rho11_b_therm_real = eigenstates_b_1 * rho11_b_therm *(eigenstates_b_1');

rho22_b_therm = zeros(size(eigen_energies_mat_b_2));
for i=1:size(rho22_b_therm,1)
    rho22_b_therm(i,i) = exp(-beta*eigen_energies_mat_b_2(i,i));
end
rho22_b_therm = rho22_b_therm./trace(rho22_b_therm)
rho22_b_therm_real = eigenstates_b_2 * rho22_b_therm * (eigenstates_b_2')

%%bases of full density matrices: rho11, rho22, rho12, rho21
n_bases_full_1 = size(bases_constrained_1,1).*size(bases_B_1,1);
bases_full_1 = cell(1,n_bases_full_1);
k=0;
for i=1:size(bases_constrained_1,1)
    for j=1:size(bases_B_1,1)
    k = k+1;
    bases_full_1{k} = [i,j];
    end
end
n_bases_full_2 = size(bases_constrained_2,1).*size(bases_B_2,1);
bases_full_2 = cell(1,n_bases_full_2);
k=0;
for i=1:size(bases_constrained_2,1)
    for j=1:size(bases_B_2,1)
    k = k+1;
    bases_full_2{k} = [i,j];
    end
end

dim_a_1 = size(bases_constrained_1, 1)
dim_a_2 = size(bases_constrained_2, 1)
dim_b_1 = size(bases_B_1, 1)
dim_b_2 = size(bases_B_2, 1)

H_a_1_full = kron(H_a_1, eye(dim_b_1));
H_a_2_full = kron(H_a_2, eye(dim_b_2));
H_b_1_full = kron(eye(dim_a_1), H_b_1);
H_b_2_full = kron(eye(dim_a_2), H_b_2);

%%density matrix evolution with lambda=0

%initial full density matrix of rho11

rho11_a = zeros(dim_a_1);
rho11_a(1,1) = 1;
%rho11_a = eye(dim_a_1);
%rho11_a = rho11_a./trace(rho11_a);
rho11_a_eigen = eigenstates_1' * rho11_a * eigenstates_1;
rho11_full = kron(rho11_a, rho11_b_therm_real);

%evolution of rho11_full

entropy_full_arr = [];
entropy_A_arr = [];
entropy_A_arr1 = [];
entropy_full_arr = [entropy_full_arr, Entropy(rho11_full)];
entropy_A_arr = [entropy_A_arr, Entropy(traceOverB(rho11_full, dim_a_1, dim_b_1))];
entropy_A_arr1 = [entropy_A_arr1, Entropy(rho11_a)];
params_11 = {H_a_1_full, H_b_1_full};
params_11_a = {H_a_1};
for n_iter = 2: total_iter
    rho11_full = RungeKutta_oneMatrix(rho11_full, @rho11fullEvolution0, params_11, dt);
    rho11_a = RungeKutta_oneMatrix(rho11_a, @rho11aEvolution0, params_11_a, dt);
    entropy_full_arr = [entropy_full_arr, Entropy(rho11_full)];
    entropy_A_arr = [entropy_A_arr, Entropy(traceOverB(rho11_full, dim_a_1, dim_b_1))];
    entropy_A_arr1 = [entropy_A_arr1, Entropy(rho11_a)];
end



%%plots
tvector=0:dt:(total_iter-1)*dt;
fig_num = 1;
%rho11_full: entropy vs time
figure(fig_num)
fig_num = fig_num + 1;
plot(tvector,entropy_full_arr,'r-.')
title('entropy of \rho_{11} vs time')

figure(fig_num)
fig_num = fig_num + 1;
plot(tvector, entropy_A_arr, 'b')
hold on
plot(tvector, entropy_A_arr1, 'r-.')
title('entropy of \rho_{11}^{A} vs time')

traceOverA(rho11_full, dim_a_1, dim_b_1)
(traceOverA(rho11_full, dim_a_1, dim_b_1) - rho11_b_therm_real) < tolerance
%eigenstates_1' * traceOverB(rho11_full, dim_a_1, dim_b_1) * eigenstates_1
rho11_a_eigen = eigenstates_1' * rho11_a * eigenstates_1



%%do perturbation to generate rho12 and rho21
%%operator  P_L
PL_1mat2 = PL(bases_constrained_2_full, bases_constrained_1);
PL_2mat1 = PL(bases_constrained_1, bases_constrained_2_full);
%%operator  X_{L+1} P_{L+2} 
XLp1PLp2_1mat2 = XLp1PLp2(bases_B_2, bases_B_1);
XLp1PLp2_2mat1 = XLp1PLp2(bases_B_1, bases_B_2);
%%operator P_L X_{L+1} P_{L+2}
PL_XLp1_PLp2_2mat1 = kron(PL_2mat1, XLp1PLp2_2mat1);
PL_XLp1_PLp2_1mat2 = kron(PL_1mat2, XLp1PLp2_1mat2);
%%two ways to generate rho21_full
rho11_a = zeros(dim_a_1);
rho11_a(1,1) = 1;
rho11_full = kron(rho11_a, rho11_b_therm_real);
rho21_full = PL_XLp1_PLp2_2mat1 * rho11_full;
rho21_full_1 = kron(PL_2mat1*rho11_a,  XLp1PLp2_2mat1*rho11_b_therm_real);
%(rho21_full - rho21_full_1) < tolerance
rho12_full = rho11_full * PL_XLp1_PLp2_1mat2;
disp('dim of rho21_full:')
size(rho21_full)
disp('dim of rho12_full:')
size(rho12_full)

% rho21_b_therm_real_proj = XLp1PLp2_2mat1*rho11_b_therm_real;
% rho21_b_therm_proj = eigenstates_b_2' * rho21_b_therm_real_proj * eigenstates_b_1;
% rho21_b_therm_proj' * rho21_b_therm_proj
% size(rho21_b_therm_proj' * rho21_b_therm_proj)


lambda = 1;
rho21_full = -imagunit* lambda *dt * rho21_full;
rho12_full = imagunit* lambda *dt *rho12_full;
rho11_flip21 = PL_XLp1_PLp2_1mat2 * rho21_full;
rho11_flip12 = rho12_full * PL_XLp1_PLp2_2mat1;
trace_HABrho21_rho12HAB_arr = [];
trace_HABrho21_rho12HAB_arr = [trace_HABrho21_rho12HAB_arr, trace(traceOverA(-imagunit*(rho11_flip21-rho11_flip12), dim_a_1, dim_b_1))];
params_21 = {H_a_1_full, H_a_2_full, H_b_1_full, H_b_2_full};
params_12 = {H_a_1_full, H_a_2_full, H_b_1_full, H_b_2_full};
for n_iter = 2: total_iter
    rho21_full = RungeKutta_oneMatrix(rho21_full, @rho21fullEvolution0, params_21, dt);
    rho11_flip21 =  PL_XLp1_PLp2_1mat2 * rho21_full;
    rho12_full = RungeKutta_oneMatrix(rho12_full, @rho12fullEvolution0, params_12, dt);
    rho11_flip12 = rho12_full * PL_XLp1_PLp2_2mat1;
    trace_HABrho21_rho12HAB_arr = [trace_HABrho21_rho12HAB_arr, trace(traceOverA(-imagunit*(rho11_flip21-rho11_flip12), dim_a_1, dim_b_1))];
end

eigenstates_b_1'*traceOverA(-imagunit*(rho11_flip21 - rho11_flip12), dim_a_1, dim_b_1)*eigenstates_b_1



%%plot
figure(fig_num)
fig_num = fig_num + 1;
plot(tvector, trace_HABrho21_rho12HAB_arr, 'r-.')
title('trace(X_{L+1} P_{L+2} \rho_{21}^{B} - \rho_{12}^{B} X_{L+1} P_{L+2} )')

%%fast fourier transform
X = trace_HABrho21_rho12HAB_arr(10:end);

Fs = 1/dt;                    % Sampling frequency
T = dt;                     % Sampling period
L = size(X,2);                     % Length of signal
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(fig_num)
fig_num = fig_num + 1;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum')
xlabel('f (Hz)')

end



%%bases
function [bases_B_1, bases_B_2]=split_bases_B(bases_B)
%split the bases of B into two subspaces with the first spin being up or
%down
bases_B_1 = [];
bases_B_2 = [];
n_states = size(bases_B, 1);
for i=1:n_states
    basis_i = bases_B(i,:);
    if basis_i(1) == -1
        bases_B_1 = [bases_B_1;basis_i];
    else
        bases_B_2 = [bases_B_2;basis_i];
    end
end
end

function bases_a = Bases(N_sites)
%to find the bases in subsystem A or subsystem B
%spin down = -1
%spin up = 1
sites_vec = 1:N_sites;
total = 2^(N_sites);
bases_a = -ones(total,N_sites);
start = 1;
for i=1:N_sites
    spinup_ind = nchoosek(sites_vec,i);
    size_spinup_ind = size(spinup_ind,1);
    for j=1:size_spinup_ind
        bases_a(start+j,spinup_ind(j,:))=1;
    end
    start = start + size_spinup_ind;
end
end

function bases_a_constrained = Bases_constrained(bases_a, N_sites)
%filter bases_a to the constrained space, which doesn't allow two adjacent
%up spins.
dim = size(bases_a,1);
bases_a_constrained = [];
for i=1:dim
    basis = bases_a(i,:);
    flag = true;
    for j=1:N_sites-1
        if basis(j)>0 && basis(j+1)>0
            flag = false;
            break;
        end
    end
    if flag 
        bases_a_constrained = [bases_a_constrained;basis];
    end
end
end

%%Hamiltonian
function H_b = Hamiltonian_b(bases)
dimH = size(bases,1);
H_b = zeros(dimH,dimH);
for i=1:dimH
    for j=i+1:dimH
        basis_i = bases(i,:);
        basis_j = bases(j,:);
        H_b(j,i) = PXP_b(basis_i, basis_j);
    end
end
H_b = H_b + H_b';
end

function Hji = PXP_b(basis_i, basis_j)
%Hamiltonian of subsystem B
%H = P_1 X_2 P_3 + ... + P_{L-1} X_L
N_sites = size(basis_i, 2);
diff = basis_i - basis_j;
diff_elements = find(diff);

Hji = 0;
if size(diff_elements,2)==1 
    if diff_elements(1) == N_sites
        if basis_i(N_sites-1) == -1
            Hji = 1;
        end
    else       
        xi = diff_elements(1);
        if basis_i(xi-1)==-1 && basis_i(xi+1)==-1
            Hji = 1;
        end
    end
end
end

function H_a = Hamiltonian_a(bases)
% to find the matrix form of Hamiltonian_a in the bases
% functions used: PXP
dimH = size(bases,1);
H_a = zeros(dimH,dimH);
for i=1:dimH
    for j=i+1:dimH
        basis_i = bases(i,:);
        basis_j = bases(j,:);
        H_a(j,i) = PXP(basis_i,basis_j);
    end
end
H_a = H_a + H_a';
end

function Hji = PXP(basis_i,basis_j)
% the matrix element between basis_i and basis_j
% H = X_1 P_2 + P_1 X_2 P_3 + ... + P_{L-1} X_L
% where L=N_sites
% used by Hamiltonian_a
N_sites = size(basis_i, 2);
diff = basis_i - basis_j;
diff_elements = find(diff);

Hji = 0;
if size(diff_elements,2)==1 
    if diff_elements(1)==1
        if basis_i(2)==-1
            Hji = 1;
        end
    elseif diff_elements(1) == N_sites
        if basis_i(N_sites-1) == -1
            Hji = 1;
        end
    else       
        xi = diff_elements(1);
        if basis_i(xi-1)==-1 && basis_i(xi+1)==-1
            Hji = 1;
        end
    end
end
end


%%write operators in the form of matrices
function PL_matrix = PL(from_subspace,to_subspace)
%matrix form of P_L, going from one subspace to the other subspace of
%subsystem A
%either from_subpace or to_subspace must be in complate shape, i.e. the
%states in subspace 2 of subsystem A cannot be truncated (the last spin)
n_sites = size(from_subspace, 2);
dim_from = size(from_subspace, 1);
dim_to = size(to_subspace, 1);
PL_matrix = zeros(dim_to, dim_from);
for fromi = 1:dim_from
    state_from = from_subspace(fromi,:);
    for toj = 1:dim_to
        state_to = to_subspace(toj,:);
        if state_from(n_sites)==-1 && state_to(n_sites)==-1
            state_abs_diff = abs(state_from - state_to);
            if sum(state_abs_diff) == 0
                PL_matrix(toj,fromi) = 1;
            end
        end
    end
end
end

function XLp1PLp2_matrix = XLp1PLp2(from_subspace, to_subspace)
%matrix form of X_{L+1}P_{L+2}, going from one subspace to the other subspace of
%subsystem A
dim_from = size(from_subspace, 1);
dim_to = size(to_subspace, 1);
XLp1PLp2_matrix = zeros(dim_to, dim_from);
for fromi = 1:dim_from
    state_from = from_subspace(fromi,:);
    for toj = 1:dim_to
        state_to = to_subspace(toj,:);
        if state_from(2) == -1 && state_to(2) == -1
            state_abs_diff = abs(state_from - state_to);
            if state_abs_diff(1) == 2 && sum(state_abs_diff) == 2
                XLp1PLp2_matrix(toj,fromi) = 1;
            end
        end
    end
end
end

%%Runge-Kutta, evolution of density matrices
function rho_next = RungeKutta_oneMatrix(rho, rho_der_fun, params, dt)
%4th order Runge-Kutta
%one matrix 
k1 = rho_der_fun(rho, params);
k2 = rho_der_fun(rho + 0.5*dt*k1, params);
k3 = rho_der_fun(rho + 0.5*dt*k2, params);
k4 = rho_der_fun(rho + dt*k3, params);
rho_next = rho + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

function rho11a_der = rho11aEvolution0(rho11a, params)
%evolution of rho11_a
%params  is a cell, {H_a_1}
global imagunit
H_a_1 = params{1};
rho11a_der = -imagunit * commutator(H_a_1, rho11a);
end

function rho11full_der = rho11fullEvolution0(rho11full, params)
%evolution of full rho11 with lambda=0
%params is a cell, {H_a_1_full, H_b_1_full}
global imagunit
H_a_1_full = params{1};
H_b_1_full = params{2};
rho11full_der = -imagunit * commutator(H_a_1_full + H_b_1_full, rho11full);
end

function rho22full_der = rho22fullEvolution0(rho22full, params)
%evolution of full rho22 with lambda=0
%params is cell, {H_a_2_full, H_b_2_full}
global imagunit
H_a_2_full = params{1};
H_b_2_full = params{2};
rho22full_der = -imagunit * commutator(H_a_2_full + H_b_2_full, rho22full);
end

function rho12full_der = rho12fullEvolution0(rho12full, params)
%evolutioin of rho12 with lambda=0
%params is a cell, {H_a_1_full, H_a_2_full, H_b_1_full, H_b_2_full}
global imagunit
H_a_1_full = params{1};
H_a_2_full = params{2};
H_b_1_full = params{3};
H_b_2_full = params{4};
rho12full_der = -imagunit * (H_a_1_full * rho12full - rho12full * H_a_2_full + H_b_1_full * rho12full - rho12full * H_b_2_full);
end

function rho21full_der = rho21fullEvolution0(rho21full, params)
%evolutioin of rho12 with lambda=0
%params is a cell, {H_a_1_full, H_a_2_full, H_b_1_full, H_b_2_full}
global imagunit
H_a_1_full = params{1};
H_a_2_full = params{2};
H_b_1_full = params{3};
H_b_2_full = params{4};
rho21full_der = -imagunit * (H_a_2_full * rho21full - rho21full * H_a_1_full + H_b_2_full * rho21full - rho21full * H_b_1_full);
end


%%operation on matrices
function com = commutator(A, B)
com = A * B - B * A;
end

function anticom = anticommutator(A, B)
anticom = A * B + B * A;
end

function matrix_B = traceOverA(matrix_full, dim_A, dim_B)
%trace over B of the tensor product of two square matrices
%trace over A: sum of all diagonal dim_B-by-dim_B blocks in the full
%matrix
matrix_B = zeros(dim_B);
for n=0:dim_A-1
    range = n*dim_B+1 : (n+1)*dim_B;
    matrix_B  = matrix_B + matrix_full(range,range);
end
end

function matrix_A = traceOverB(matrix_full, dim_A, dim_B)
%trace over A of the tensor product of two square matrices
%trace over B: (i,j) element of matrix_A is equal to the trace of the
%(i,j) block of dimension dim_B-by-dim_B
matrix_A = zeros(dim_A);
for i=1:dim_A
    ri = (i-1)*dim_B+1 : i*dim_B;
    for j=1:dim_A
        rj = (j-1)*dim_B+1 : j*dim_B;
        matrix_A(i,j) = trace(matrix_full(ri,rj));
    end
end
end

function entropy = Entropy(densityMatrix)
% Boltzmann constant k_B = 1
[eigenstates,diagD] = eig(densityMatrix);
diagd = diag(diagD); 
dim = size(diagD,1);
entropy = 0;
for i = 1:dim
    entropy = entropy - diagd(i).*log(diagd(i));
end
end
