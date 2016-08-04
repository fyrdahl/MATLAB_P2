

Mz(1,:) = 1-2*exp(-[1:4500]/282);

 TI = [100:50:600,700:100:2500];
 
 opts.fiteff = 1;
 [pd r1 eff res] = qmap_t1_fit_ir(abs(Mz(1,:)), TI, opts);
 1/r1
 
 figure, plot(TI,Mz(1,:))