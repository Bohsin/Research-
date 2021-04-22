function g_2_3
mcmc = [2.,1.,2.,1.,1.99999,0.999986,1.99955,0.99915,1.98685,0.985644,1.75053,0.749684,0.323919,0.128411,0.047957,0.018139,0.006616,0.002332,0.000794,0.000337,0.000088,0.000025];
Wudist = [1.5,1.5,1.5,1.49999,1.49997,1.49987,1.49941,1.49737,1.48832,1.44945,1.30354,0.938098,0.49004,0.206672,0.0798417,0.0298935,0.0110679,0.0040812,0.00150269,0.000552982,0.000203455,0.00007485];
len = size(mcmc,2);
levels = 1:1:len;

figure()
plot(levels,mcmc,"r-*")
hold on
plot(levels,Wudist,"b--")
legend("g=2/3(MCMC)","g=2/3(Yongshi Wu's result)")
title("Particle occupation numbers vs energy levels")


figure()
plot(levels, abs(mcmc-Wudist),'--')
title("Difference between MCMC and Wu's result")
legend("Original data")

mcmc_shift_back_1 = [mcmc(2:len), 0];
disp(mcmc)
disp(mcmc_shift_back_1)
mcmc_shift_forward_1 = [1, mcmc(1:len-1)];
disp(mcmc_shift_forward_1)

mcmc_back1_avg = (mcmc+mcmc_shift_forward_1)./2;
figure();
plot(levels,mcmc_back1_avg,'r*-')
hold on
plot(levels,Wudist,'b--')
legend("g=2/3(Backward averaged MCMC)","g=2/3(Yongshi Wu's result)")
title("Particle occupation numbers vs energy levels")

figure()
plot(levels, abs(mcmc_back1_avg-Wudist),'--')
title("Difference between MCMC and Wu's result")
legend("backward-1 average")

mcmc_forward1_avg = (mcmc+mcmc_shift_back_1)./2;
figure();
plot(levels,mcmc_forward1_avg,'r*-')
hold on
plot(levels,Wudist,'b--')
legend("g=2/3(forward averaged MCMC)","g=2/3(Yongshi Wu's result)")
title("Particle occupation numbers vs energy levels")

figure()
plot(levels, abs(mcmc_forward1_avg-Wudist),'--')
title("Difference between MCMC and Wu's result")
legend("forward-1 average")

mcmc_adjacent3_avg = (mcmc+mcmc_shift_back_1+mcmc_shift_forward_1)./3;
figure();
plot(levels,mcmc_adjacent3_avg,'r*-')
hold on
plot(levels,Wudist,'b--')
legend("g=2/3(adjacent-3 averaged MCMC)","g=2/3(Yongshi Wu's result)")
title("Particle occupation numbers vs energy levels")

figure()
plot(levels, abs(mcmc_adjacent3_avg-Wudist),'--')
title("Difference between MCMC and Wu's result")
legend("Adjacent-3 average")


end