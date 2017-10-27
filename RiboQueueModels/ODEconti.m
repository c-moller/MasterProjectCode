function [xdot, outpara] = ODEconti(t, y, P1)

Amino = y(1); %Amino acids
Ribo = y(2);  %Ribosomes
TC_Meta = y(3); %TC Meta.
Meta = y(4);  %Metabolic Proteins

mRNA_Meta = P1.mRNAPool(1);  %mRNA Metabolic Proteins
TC = y(3);%

%Rates

ribo_density = TC./mRNA_Meta;
RiboRBS = P1.RBS * Ribo;
k_codon_slow =P1.k_cs;% * Amino;
k_codon_fast = P1.k_cf;% * Amino;
k_trans_mRNA = (k_codon_slow/k_codon_fast);
ribo_density_max = P1.x_cs; %Max queue length, assume length of mRNA in codon increments so L = 100 is mRNA can have 100TC on it
max_density = ribo_density_max/2;
queue_function = ribo_density_max^P1.n ./(ribo_density_max^P1.n + ribo_density^P1.n);    %is 1 when no queue, n steepness => higher value queue close to 0
k_init = P1.RBS .* Ribo .* P1.mRNAPool .* queue_function; %Ribosome initiation rate
k_Prod_Meta = (ribo_density*k_trans_mRNA - k_codon_slow).*(max_density^P1.a/(max_density^P1.a + ribo_density^P1.a))+  k_codon_slow;


%-------------------------------------------------------------------------------------------------------------------------------------------------
dAmino = 0;% k_Acat - k_Adep - dilution*Amino; %k_ProDeg
dRibo =  0;% k_ProdTotal - Total_TC  - dilution*Ribo; %k_Prod_Ribo
dTC_Meta = k_init*mRNA_Meta - k_Prod_Meta;% - dilution*TC_Meta;
dMeta = k_Prod_Meta;% - dilution*Meta;

%-------------------------------------------------------------------------------------------------------------------------------------------------

xdot = [dAmino; %x(1)
    dRibo; %x(2)
    dTC_Meta; %x(3)
    dMeta; %x(4)
    ];
outpara = [k_init;
    k_Prod_Meta;
    queue_function;
    ribo_density;
    RiboRBS;
    ribo_density_max;
    ];

end
