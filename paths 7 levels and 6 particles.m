gamma=1;
semionpath=gamma*[0., 1., 1., 1.33333, 1.5, 1.8, 1.71429, 2.14286, 2.125, 2.25, 2.25, 2.42857, 2.14286, 2.4, 2.25, 2., 2., 2., 1.]
fermionpath=gamma*[0., 1., 1.2, 1.6, 2., 2.33333, 2.54545, 2.85, 3.06667, 3.16667, 3.37778, 3.45, 3.45455, 3.5, 3.5, 3.2, 3.2, 3., 2.]


a=linspace(0,18,19);

figure
plot(a, semionpath, "bo")
hold on
plot(a,fermionpath, "r*")
hold off
legend('Paths of semion','Paths of spinful fermion')
