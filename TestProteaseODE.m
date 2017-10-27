tic
clear all
tps =  20;
P1 = struct('s0', 1000, 'mRNAPool',[1 1],...
    'x_cs',[4 5], 'k_cs',[5 5], 'k_cf',[25 25],...
    'RBS',[2 2],'L',[5 10], 'epsi', 0.01,'Acat',20, 'km',1);

%Initial values
y0 = [10 5 0 0];
%ODE solver
[T1, x] = ode23(@(t,y) ProteaseODE(t, y, P1), [0 tps], y0);
%Parameter Extract
for i = 1: length(T1)
    [x2, paraout] = ProteaseODE(T1(i),x(i,:),P1);
    Time(i,1)=T1(i);
    para_init(i,1) = paraout(1);
    para_release(i,1) = paraout(2);
    para_RBSr(i,1) = paraout(3);
end
xplot =T1;
y1plot = x(:,3); %[TC]
y2plot = x(:,1); %[AA]
y3plot =  para_release;
y4plot = para_init;

figure
subplot(3,1,1);
plot(xplot,y1plot, ' r'),hold on
xlabel('Time');
ylabel('[TC]');
title({'k_{cs}=5 , k_{cf}=25, mRNA length=10, Slow Codon Location=5, Km=1, Acat=49';'Initial values:AA Pool=10, Ribo Pool=5, TC=0'});
subplot(3,1,2);
plot(xplot,y2plot, ' b'),hold on
xlabel('Time');
ylabel('[Amino Acids]');
subplot(3,1,3);
plot(xplot,y3plot),hold on
plot(xplot,y4plot)
xlabel('Time');
ylabel('Rates');
legend('Release', 'Initiation');
toc