function [m_632p,CI,m_list] = kyu_632plusbootstrap_prob(m_trn,m_test,Pev)

% calculate 632+ bootstrap estimate of a probability estimate
% m_trn is a probability estimate of a model trained/tested in the entire
% data (a single value)
% m_test is a (Nx1) vector of out-of-bag P from N bootstrap samples
% Pev : average event rate in the entire datset

% calculate the bias of resubtitution (BR)
Nrand = length(m_test);
BR = zeros(Nrand,1);
R = zeros(Nrand,1);
alpha = zeros(Nrand,1);
for i = 1:Nrand
    % calculate relative overfitting rate
    % a function R is built upon these assumptions:
    % 1: R is zero when P_trn ~ P_test
    % 2: R increases when P_trn ->0 or P_test -> 1
    % 3: R decreases when P_trn -> Pev (average event rate)
    % a recti-linear function 1/f is made to meet the condition 2 and 3
    f = 1-abs(m_trn-Pev);
    R(i) = abs(m_test(i)-m_trn)/f;
    %R(i) = 1; 
    alpha(i) = 0.632/(1-0.368*R(i));
end
% R should be between 0 and 1. If not, scale it down so that
% max(R) == 1
if max(R) >1
    R = R/max(R);
end
alpha = arrayfun(@(x) 0.632/(1-0.368*x), R);
BR = alpha.*(m_test-m_trn);
m_632p = m_trn + nanmean(BR);
m_list = m_trn + BR;
CI = 2*nanstd(BR)/sqrt(Nrand);

