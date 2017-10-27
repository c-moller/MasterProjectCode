function [xdot, outpara] = ODErho(t, y, P1)

Amino = y(1);              %Amino acids
Ribo = y(2);  %Ribosomes
TR = y(3); %TC Meta.
Proteins = y(4);  %Metabolic Proteins

mRNA = P1.mRNAPool(1);  %mRNA Metabolic Proteins
TR = y(3);%

%Rates

ribo_density = TR./P1.mRNAPool;
Ribo_RBS = P1.RBS * Ribo;
k_trans_slow =P1.k_cs;% * Amino;
k_trans_fast = P1.k_cf;% * Amino;
%ribo_density_max = P1.x_cs + (P1.L - P1.x_cs)*(k_trans_slow/k_trans_fast);
ribo_density_max = P1.x_cs; %Max queue length, assume length of mRNA in codon increments so L = 100 is mRNA can have 100TC on it

queue_function = ribo_density_max^P1.n ./(ribo_density_max^P1.n + ribo_density^P1.n);    %is 1 when no queue, n steepness => higher value queue close to 0
k_init = P1.RBS .* Ribo .* P1.mRNAPool .* queue_function; %Ribosome initiation rate
max_trans = ribo_density_max/2;
k_release =  k_trans_slow.*(ribo_density^P1.a/(max_trans^P1.a + ribo_density^P1.a));


%-------------------------------------------------------------------------------------------------------------------------------------------------
dAmino = 0;
dRibo =  0;
dTR = k_init*mRNA - k_release;
dProtein = k_release;

%-------------------------------------------------------------------------------------------------------------------------------------------------

xdot = [dAmino; %x(1)
    dRibo; %x(2)
    dTR; %x(3)
    dProtein; %x(4)
    ];
outpara = [k_init;
    k_release;
    queue_function;
    ribo_density;
    Ribo_RBS;
    ribo_density_max;
    ];

end
