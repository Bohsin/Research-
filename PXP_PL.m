function PXP_PL
%%% P_L in the basis of energy eigenvectors

%%parameters
N_sites = 5; %num of sites in subsystem A

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
eigenstates_1
eigen_energies_1 = diag(eigen_energies_mat_1)
[eigenstates_2,eigen_energies_mat_2] = eig(H_a_2);
%eigenstates_2
%eigen_energies_2 = diag(eigen_energies_mat_2)

PL_spinconfig = PL_insubspace1(bases_constrained_1)

PL_eigenvector = eigenstates_1' * PL_spinconfig * eigenstates_1
diag(PL_eigenvector)

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

function PL_matrix = PL_insubspace1(subspace1)
%matrix form of P_L in subspace 1
n_sites = size(subspace1, 2);
dim = size(subspace1, 1);
PL_matrix = zeros(dim, dim);
for i=1:dim
    if subspace1(i,n_sites) == -1
        PL_matrix(i,i) = 1;
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
