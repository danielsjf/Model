clc
close all

load Results_stability2days(2).mat

Nr = size(resultsOld,2);

Profit = zeros(Nr,1);
Profit_opt = zeros(Nr,1);
NCHP = zeros(Nr,1);
R_b = zeros(Nr,1);
R_ir = zeros(Nr,1);
S_gen = zeros(Nr,1);
S_red = zeros(Nr,1);
N = zeros(Nr,1);
S = zeros(Nr,3);

for k = 1:Nr
    Profit(k) = resultsOld(k).a.obj;
    R_b(k) = sum(resultsOld(k).a_OPT.R_b);
    R_ir(k) = sum(resultsOld(k).a_OPT.R_ir);
    Profit_opt(k) = resultsOld(k).a_OPT.obj;
    NCHP(k) = resultsOld(k).a_NCHP.obj;
    S_gen(k) = resultsOld(k).Par.S_gen;
    S_red(k) = resultsOld(k).Par.S_red;
    N(k) = resultsOld(k).Par.N;
    S(k,N(k)/6) = S_gen(k) + S_red(k)/100;
end
N_s = unique(N);
N_t = numel(N_s);
Profit_opt_s = sort(unique(roundn(Profit_opt,-1)),'descend');

Profit6 = Profit(find(S(:,1)),1);
S6 = S(find(S(:,1)),1);
S6_s = changem(S6,1:numel(unique(S6)),unique(S6));

Profit12 = Profit(find(S(:,2)),1);
S12 = S(find(S(:,2)),2);
S12_s = changem(S12,1:numel(unique(S12)),unique(S12));

Profit18 = Profit(find(S(:,3)),1);
S18 = S(find(S(:,3)),3);
S18_s = changem(S18,1:numel(unique(S18)),unique(S18));

figure();
axes1 = axes('XTickLabel',unique(S6),...
    'XTick',1:numel(unique(S6)));
box(axes1,'on');
hold(axes1,'all');
plot(S6_s,Profit6, '.', S6_s, Profit_opt_s(1)*ones(numel(S6_s),1),':')
title(['Stability for different numbers of generated scenarios with ',num2str(N_s(1)),' units']);

figure();
axes2 = axes('XTickLabel',unique(S12),...
    'XTick',1:numel(unique(S12)));
box(axes2,'on');
hold(axes2,'all');
plot(S12_s,Profit12, '.', S12_s, Profit_opt_s(2)*ones(numel(S12_s),1),':')
title(['Stability for different numbers of generated scenarios with ',num2str(N_s(2)),' units']);

figure();
axes3 = axes('XTickLabel',unique(S18),...
    'XTick',1:numel(unique(S18)));
box(axes3,'on');
hold(axes3,'all');
plot(S18_s,Profit18, '.', S18_s, Profit_opt_s(3)*ones(numel(S18_s),1),':')
title(['Stability for different numbers of generated scenarios with ',num2str(N_s(3)),' units']);