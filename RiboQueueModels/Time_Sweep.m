function Time_Sweep(P2)

tps=P2.tps;
y = [0.1, 0.5, 1, 1.5, 2];
endSweep = length(y);
j=1;
figure
for var = 1:endSweep
    
    P1 = struct('s0', 1000, 'mRNAPool', 1,...
        'x_cs',3, 'k_cs', 5, 'k_cf', 25,...
        'ribo_density_max',4,...
        'RBS',y(j)*5,'L', 10, 'n', 10, 'a', 20, ...
        'CellMass', 100);
    
    %Initial values
    %y0 = ones(1,(1+2*P1.ngenes));
    y0 = [0 5 0 0];
    
    %ODE solver
    [T1, x] = ode23(@(t,y) ODErho(t, y, P1), [0 tps], y0);
    
    %Parameter Extract
    
    for i = 1: length(T1)
        [x2, paraout] = ODErho(T1(i),x(i,:),P1);
        Time(i,var)=T1(i);
        para_init(i,var) = paraout(1);
        para_prod(i,var) = paraout(2);
        para_Queue(i,var) = paraout(3);
        para_density(i,var) = paraout(4);
        para_RBSr(i,var) = paraout(5);
        para_MaxDens(i,var) = paraout(6);
    end
    last_TC(var,:)= x(i,3);
    last_init(var,:)= para_init(i,var);
    last_prod(var,:)= para_prod(i,var);
    last_density(var,:)= para_density(i,var);
    last_RBSr(var,:)= para_RBSr(i,var);
    last_MaxDens(var,:) = para_MaxDens(i,var);
    
    x1plot = x(:, 3);
    subplot(3,1,1);
    plot(T1,x1plot, '-- ','LineWidth', 2),hold on
    xlabel('Time');
    ylabel('TR');
    
    j=j+1;
end
title({'k_{cs}=5, k_{cf}=25, mRNA length=10, Slow Codon Loc=3, Hill Coeff=10 and 20, Ribo Density Max=4';'Initial values: Ribo Pool=5, TR=0'});
legend('RBS=0.1*k_{cs}', 'RBS=0.5*k_{cs}','RBS=1*k_{cs}','RBS=1.5*k_{cs}','RBS=2*k_{cs}');


%  Transient SWEEP GRAPH
x2plot = para_init;
last_x2plot = last_init;
x3plot = para_prod;
last_x3plot = last_prod;

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
title({'k_{cs}=5, k_{cf}=25, mRNA length=10, Slow Codon Loc=7, Hill Coeff=10 and 20, Ribo Density Max=4';'Initial values: Ribo Pool=5, TC=0'});
subplot(3,1,3);
plot(Time,x3plot, '-- ','LineWidth', 2),hold on
xlabel('Time');
ylabel('Release Rate');
legend('RBS=0.1*k_{cs}', 'RBS=0.5*k_{cs}','RBS=1*k_{cs}','RBS=1.5*k_{cs}','RBS=2*k_{cs}');
end




