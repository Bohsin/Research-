function semion

mcmc = [2.,2.,2.,2.,2.,1.99998,1.99994,1.99969,1.99842,1.98789,1.90839,1.40279,0.468753,0.152231,0.051922,0.019353,0.007236,0.002294,0.000818,0.000211,0.000071,0.000017];
Wudist = [2.,2.,2.,2.,2.,1.99999,1.99993,1.99951,1.99636,1.97358,1.8264,1.27221,0.580452,0.221765,0.0820198,0.0301954,0.0111093,0.00408695,0.00150351,0.00055311,0.000203478,0.0000748553];
spinfulFermidist = [1.99998,1.99994,1.99985,1.99959,1.99889,1.997,1.99186,1.97803,1.94138,1.84828,1.63515,1.24492,0.755079,0.364849,0.151716,0.0586241,0.0219738,0.00814023,0.00300235,0.00110555,0.000406852,0.000149692];
len = size(mcmc,2)
levels = 1:1:len;
figure();
plot(levels,mcmc,'r*')
hold on
plot(levels,Wudist,'b--')
legend("Semion(MCMC)","Semion(Yongshi Wu's result)")
title("Particle occupation numbers vs energy levels")
hold on
plot(levels,spinfulFermidist,'black-')
legend("Semion(MCMC)","Semion(Yongshi Wu's result)","Spinful fermion")

figure()
plot(levels, abs(mcmc-Wudist),'--')
title("Difference between MCMC and Wu's result")
legend("Original data")

mcmc_shift_back_1 = [mcmc(2:len), 0];
disp(mcmc)
disp(mcmc_shift_back_1)
mcmc_shift_forward_1 = [2, mcmc(1:len-1)];
disp(mcmc_shift_forward_1)

mcmc_back1_avg = (mcmc+mcmc_shift_forward_1)./2;
figure();
plot(levels,mcmc_back1_avg,'r*')
hold on
plot(levels,Wudist,'b--')
legend("Semion(Backward averaged MCMC)","Semion(Yongshi Wu's result)")
title("Particle occupation numbers vs energy levels")
%hold on
%plot(levels,spinfulFermidist,'black-')
%legend("Semion(MCMC)","Semion(Yongshi Wu's result)","Spinful fermion")

figure()
plot(levels, abs(mcmc_back1_avg-Wudist),'--')
title("Difference between MCMC and Wu's result")
legend("backward-1 average")

mcmc_forward1_avg = (mcmc+mcmc_shift_back_1)./2;
figure();
plot(levels,mcmc_forward1_avg,'r*')
hold on
plot(levels,Wudist,'b--')
legend("Semion(forward averaged MCMC)","Semion(Yongshi Wu's result)")
title("Particle occupation numbers vs energy levels")
%hold on
%plot(levels,spinfulFermidist,'black-')
%legend("Semion(MCMC)","Semion(Yongshi Wu's result)","Spinful fermion")

figure()
plot(levels, abs(mcmc_forward1_avg-Wudist),'--')
title("Difference between MCMC and Wu's result")
legend("forward-1 average")

mcmc_adjacent3_avg = (mcmc+mcmc_shift_back_1+mcmc_shift_forward_1)./3;
figure();
plot(levels,mcmc_adjacent3_avg,'r*')
hold on
plot(levels,Wudist,'b--')
legend("Semion(adjacent-3 averaged MCMC)","Semion(Yongshi Wu's result)")
title("Particle occupation numbers vs energy levels")
%hold on
%plot(levels,spinfulFermidist,'black-')
%legend("Semion(MCMC)","Semion(Yongshi Wu's result)","Spinful fermion")

figure()
plot(levels, abs(mcmc_adjacent3_avg-Wudist),'--')
title("Difference between MCMC and Wu's result")
legend("Adjacent-3 average")
end



