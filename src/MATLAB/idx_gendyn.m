function [XL, RA, XD, XDP, XDPP, TDP, TDPP, XQ, XQP, XQPP, TQP, TQPP,...
    H, D0, D1, S1, S2, RGOV, T1, T2, T3, T4, T5, F] = idx_gendyn

%IDX_GENDYR   Defines constants for named column indices to gen matrix.
%   Example:
%
%   [XL, RA, XD, XDP, XDPP, TDP, TDPP, XQ, XQP, XQPP, TQP, TQPP,...
%     H, D0, D1, S1, S2, RGOV, T1, T2, T3, T4, T5, F] = idx_gendyn;
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    R  = gen(4, RGOV);   % get the governor droop of generator 4
%    gen(:, xd) = 0;      % set to zero the d axis synchronous resistance of all gens
% 
%   The index, name and meaning of each column of the gen matrix is given
%   below:
%
%      26. leakage reactance x_l(pu),
%      27. resistance r_a(pu),
%      28. d-axis sychronous reactance x_d(pu),
%      29. d-axis transient reactance x'_d(pu),
%      30. d-axis subtransient reactance x"_d(pu),
%      31. d-axis open-circuit time constant T'_do(sec),
%      32. d-axis open-circuit subtransient time constant T"_do(sec),
%      33. q-axis sychronous reactance x_q(pu),
%      34. q-axis transient reactance x'_q(pu),
%      35. q-axis subtransient reactance x"_q(pu),
%      36. q-axis open-circuit time constant T'_qo(sec),
%      37. q-axis open circuit subtransient time constant
%                T"_qo(sec),
%      38. inertia constant H(sec),
%      39. damping coefficient d_o(pu),
%      40. dampling coefficient d_1(pu),
%      41. saturation factor S(1.0)
%      42. saturation factor S(1.2)
%      43. governor droop R(pu),
%      44. time constant T_1(s),
%      45. time constant T_2(s),
%      46. time constant T_3(s),
%      47. time constant T_4(s),
%      48. time constant T_5(s),
%      49. F
%
%   See also DEFINE_CONSTANTS.


%% define the indices
XL     = 26;    %% leakage reactance x_l(pu)
RA     = 27;    %% resistance r_a(pu)
XD     = 28;    %% d-axis sychronous reactance x_d(pu)
XDP    = 29;    %% d-axis transient reactance x'_d(pu)
XDPP   = 30;    %% d-axis subtransient reactance x"_d(pu)
TDP    = 31;    %% d-axis open-circuit time constant T'_do(sec)
TDPP   = 32;    %% d-axis open-circuit subtransient time constant T"_do(sec)
XQ     = 33;    %% q-axis sychronous reactance x_q(pu)
XQP    = 34;   %% q-axis transient reactance x'_q(pu)
XQPP   = 35;   %% q-axis subtransient reactance x"_q(pu)
TQP    = 36;   %% q-axis open-circuit time constant T'_qo(sec)
TQPP   = 37;   %% q-axis open circuit subtransient time constant T"_qo(sec)
H      = 38;   %% inertia constant H(sec)
D0     = 39;   %% damping coefficient d_o(pu)
D1     = 40;   %% damping coefficient d_1(pu)
S1     = 41;   %% saturation factor S(1.0)
S2     = 42;   %% saturation factor S(1.2)
RGOV   = 43;   %% governor droop R(pu)
T1     = 44;   %% time constant T_1(s)
T2     = 45;   %% time constant T_2(s)
T3     = 46;   %% time constant T_3(s)
T4     = 47;   %% time constant T_4(s)
T5     = 48;   %% time constant T_5(s)
F      = 49;   %% F
end

