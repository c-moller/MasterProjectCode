tic
% Parameter Sweep
clear
tps =  50;
y = 0:1:50;
%y = [0.5, 1, 1.5, 2, 2.5, 3];
endSweep = length(y);
j=1;

%Parameter Extract
for var = 1:endSweep
    
    P1 = struct('s0', 1000, 'mRNAPool', 1,...
        'x_cs',5, 'k_cs', 5, 'k_cf', 25,...
        'RBS',2,'L', 10, 'epsi', 1e-2,'Acat',y(j), 'km',1);
    
    %Initial values
    y0 = [0 1 0 0];
    %ODE solver
    [T1, x] = ode23(@(t,y) AminoODE(t, y, P1), [0 tps], y0);
    
    %Parameter Extract
    for i = 1: length(T1)
        [x2, paraout] = AminoODE(T1(i),x(i,:),P1);
        Time(i,j)=T1(i);
        para_init(i,j) = paraout(1);
        para_release(i,j) = paraout(2);
        para_RBSr(i,j) = paraout(3);
    end
    last_TC(j,:)= x(i,3);
    last_Amino(j,:)= x(i,1);
    last_init(j,:)= para_init(i,j);
    last_release(j,:)= para_release(i,j);
    last_RBSr(j,:)= para_RBSr(i,j);
    
    j=j+1;
end
xplot =y;
y1plot = last_TC;
y2plot = last_Amino;
y3plot =  last_release;

figure
subplot(3,1,1);
plot(xplot,y1plot, '-- +  r'),hold on
xlabel('Rate AAcat');
ylabel('TC');
title({'k_{cs}=5 , k_{cf}=25, mRNA length=10, Slow Codon Location=5, Km=1';'Initial values:AA Pool=0, Ribo Pool=1, TC=0'});
subplot(3,1,2);
plot(xplot,y2plot, '-- + k'),hold on
xlabel('Rate AAcat');
ylabel('Amino Acids');
subplot(3,1,3);
plot(xplot,y3plot, '-- +  g'),hold on
xlabel('Rate AAcat');
ylabel('Release Rate');

%% Time Plot
clear all
tps =  20;
P1 = struct('s0', 1000, 'mRNAPool', 1,...
    'x_cs',5, 'k_cs', 5, 'k_cf', 25,...
    'RBS',2,'L', 10, 'epsi', 0.01,'Acat',49, 'km',1);

%Initial values
y0 = [10 5 0 0];
%ODE solver
[T1, x] = ode23(@(t,y) AminoODE(t, y, P1), [0 tps], y0);
%Parameter Extract
for i = 1: length(T1)
    [x2, paraout] = AminoODE(T1(i),x(i,:),P1);
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