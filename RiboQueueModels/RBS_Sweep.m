function RBS_Sweep(P2)

tps=1000;
y = 0:0.1:10;
endSweep = length(y);
j=1;
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
    j = j+1;
end
figure
xplot =last_RBSr;
y1plot = last_TC;
y2plot = last_init;
y3plot =  last_prod;
subplot(3,1,1);
plot(xplot, P1.x_cs*ones(size(xplot)), 'LineStyle','--','LineWidth', 2),hold on
plot(xplot,y1plot, '+ r')
xlabel('RBS Binding Rate (k_{RBS}*Ribo_{free})');
ylabel('TR');
legend('Slow Codon Location ');
title({'k_{cs}=5, k_{cf}=25, mRNA length=10, Slow Codon Loc=3, Hill Coeff=10 and 20, Ribo Density Max=4';'Initial values: Ribo Pool=5, TR=0'});
subplot(3,1,2);
plot(xplot, P1.k_cs*ones(size(xplot)),'m', 'LineStyle','--','LineWidth', 2),hold on
plot(xplot,y2plot, '+ b')
xlabel('RBS Binding Rate (k_{RBS}*Ribo_{free})');
ylabel('Initiation Rate');
legend('Rate Slow Codon ');
subplot(3,1,3);
plot(xplot, P1.k_cs*ones(size(xplot)),'m', 'LineStyle','--','LineWidth', 2),hold on
plot(xplot,y3plot, '+ k ')
axis([0 100 0 8])
xlabel('RBS Binding Rate (k_{RBS}*Ribo_{free})');
ylabel('Release Rate');
legend('Rate Slow Codon');

