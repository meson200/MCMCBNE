function [m_632p,CI,m_list] = kyu_632plusbootstrap(m_trn,m_test,gamma)

% calculate 632+ bootstrap estimate of a performance index
% input:
% m_trn: performance metric of a model trained/tested in the entire data (a single value)
% m_test: (Nx1) vector of out-of-bag performance from N bootstrap samples
% gamma: no-information performance index, depends on which metric (eg.AUC:gamma = 0.5)
%
% output:
% m_632p: .632+ averaged metric
% CI: 95% CI of m_632p
% m_list(b): a list of .632+ adjusted metric for each boostrap replicate
% [1-a(b)AUC(x,x)k+a(b)AUC(x^b,x^b(0))]


% calculate the bias of resubtitution (BR)
Nrand = length(m_test);
BR = zeros(Nrand,1);
R = zeros(Nrand,1);
for i = 1:Nrand
    % calculate relative overfitting rate
    if m_test(i)>gamma && m_test(i)<m_trn
        R(i) = (m_trn-m_test(i))/(m_trn-gamma);
    elseif m_test(i) < gamma
        R(i) = 1;
    else
        R(i) = 0;
    end
    alpha = 0.632/(1-0.368*R(i));
    BR(i) = alpha*(m_test(i)-m_trn);
end
m_632p = m_trn + nanmean(BR);
m_list = m_trn + BR;
eq_sample_size = numel(find(~isnan(BR)));
CI = 2*nanstd(BR)/sqrt(eq_sample_size);

