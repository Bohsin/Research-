function s = entropy(rho)
eigen_values = eig(rho);
s=0;
for i=1:size(rho,1)
    s = s + eigen_values(i)*log(eigen_values(i));
end
s = -s;
end