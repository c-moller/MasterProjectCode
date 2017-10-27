function Transient_Sweep(P2)

tps=P2.tps;
y = [0.1, 0.5, 1, 1.5, 2];
endSweep = length(y);
j=1;
figure
for var = 1:endSweep
    
P1 = struct('s0', 1000, 'mRNAPool', 1,...
	'x_cs',25, 'k_cs', 5, 'k_cf', 25,...
    'RBS',y(j)*5,'L', 50, 'epsi', 10.^-2, 'km',5);
    
    %Initial values
    y0 = [0 1 0 0];
    
    %ODE solver
    [T1, x] = ode23(@(t,y) ODEdistRibo(t, y, P1), [0 tps], y0);
    %[T1, x] = ode23(@(t,y) ProteaseODE(t, y, P1), [0 tps], y0);
    
    %Parameter Extract
    
    for i = 1: length(T1)
        [x2, paraout] = ODEdistRibo(T1(i),x(i,:),P1);
        %[x2, paraout] = ProteaseODE(T1(i),x(i,:),P1);
        Time(i,j)=T1(i);
        para_init(i,j) = paraout(1);
        para_release(i,j) = paraout(2);
        para_RBSr(i,j) = paraout(3);
    end
    last_init(j,:)= para_init(i,j);
    last_release(j,:)= para_release(i,j);
    last_RBSr(j,:)= para_RBSr(i,j);
    subplot(3,1,1); 
    x1plot = x(:, 3);
    plot(T1,x1plot, '-- ','LineWidth', 2),hold on
    xlabel('Time');
    ylabel('TR');
       
    j=j+1;
end
title({'k_{cs} =5 , k_{cf} =25, mRNA length =50, Slow Codon Location =25';'Initial values:AA Pool=0, Ribo Pool=1, TC=5'});
legend('RBS=0.1*k_{cs}', 'RBS=0.5*k_{cs}','RBS=1*k_{cs}','RBS=1.5*k_{cs}','RBS=2*k_{cs}'); 


%  Transient SWEEP GRAPH
x2plot = para_init;
last_x2plot = last_init;
x3plot = para_release;
last_x3plot = last_release;

for l=1:(j-1)
    for k=10:length(x2plot)
        check = x2plot(k,l)==0;
        check3 = x3plot(k,l)==0;
        check2 = Time(k,l)==0;
        x2plot(k,l) = x2plot(k,l)+check*last_x2plot(l);
        x3plot(k,l) = x3plot(k,l)+check3*last_x3plot(l);
        Time(k,l) = Time(k,l)+check2*tps;
    end
end

subplot(3,1,2);
plot(Time,x2plot, '-- ','LineWidth', 2),hold on
xlabel('Time');
ylabel('Initiation Rate');
legend('RBS=0.1*k_{cs}', 'RBS=0.5*k_{cs}','RBS=1*k_{cs}','RBS=1.5*k_{cs}','RBS=2*k_{cs}');
subplot(3,1,3); 
plot(Time,x3plot, '-- ','LineWidth', 2),hold on
xlabel('Time');
ylabel('Release Rate');
legend('RBS=0.1*k_{cs}', 'RBS=0.5*k_{cs}','RBS=1*k_{cs}','RBS=1.5*k_{cs}','RBS=2*k_{cs}');
end




