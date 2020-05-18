function mpc = case39mac(mpc)
%case39mac Loads the machine and governor data for the 39 bus test case

% Machine data format
%       1. machine number,
%       2. bus number,
%       3. base mva,
%       4. leakage reactance x_l(pu),
%       5. resistance r_a(pu),
%       6. d-axis sychronous reactance x_d(pu),
%       7. d-axis transient reactance x'_d(pu),
%       8. d-axis subtransient reactance x"_d(pu),
%       9. d-axis open-circuit time constant T'_do(sec),
%      10. d-axis open-circuit subtransient time constant
%                T"_do(sec),
%      11. q-axis sychronous reactance x_q(pu),
%      12. q-axis transient reactance x'_q(pu),
%      13. q-axis subtransient reactance x"_q(pu),
%      14. q-axis open-circuit time constant T'_qo(sec),
%      15. q-axis open circuit subtransient time constant
%                T"_qo(sec),
%      16. inertia constant H(sec),
%      17. damping coefficient d_o(pu),
%      18. dampling coefficient d_1(pu),
%      19. saturation factor S(1.0)
%      20. saturation factor S(1.2)
% note: all the following machines use transient reactance model
%
% Governor data format
%       1. machine number,
%       2. bus number,
%       3. governor droop R(pu),
%       4. time constant T_1(s),
%       5. time constant T_2(s),
%       6. time constant T_3(s),
%       7. time constant T_4(s),
%       8. time constant T_5(s),
%       9. F


macdat = [
1  30 100 0.0030 0.0  0.0200  0.0060  0  7.00  0  0.019  0.0080 0 0.70 0 500.0  4.00 0 0 0;
2  31 100 0.0350 0.0  0.2950  0.0697  0  6.56  0  0.282  0.1700 0 1.50 0  30.3  9.75 0 0 0;
3  32 100 0.0304 0.0  0.2495  0.0531  0  5.70  0  0.237  0.0876 0 1.50 0  35.8 10.00 0 0 0;
4  33 100 0.0295 0.0  0.2620  0.0436  0  5.69  0  0.258  0.1660 0 1.50 0  28.6 10.00 0 0 0;
5  34 100 0.0540 0.0  0.6700  0.1320  0  5.40  0  0.620  0.1660 0 0.44 0  26.0  3.00 0 0 0;
6  35 100 0.0224 0.0  0.2540  0.0500  0  7.30  0  0.241  0.0814 0 0.40 0  34.8 10.00 0 0 0;
7  36 100 0.0322 0.0  0.2950  0.0490  0  5.66  0  0.292  0.1860 0 1.50 0  26.4  8.00 0 0 0;
8  37 100 0.0280 0.0  0.2900  0.0570  0  6.70  0  0.280  0.0911 0 0.41 0  24.3  9.00 0 0 0;
9  38 100 0.0298 0.0  0.2106  0.0570  0  4.79  0  0.205  0.0587 0 1.96 0  34.5 14.00 0 0 0;
10 39 100 0.0125 0.0  0.1000  0.0310  0 10.20  0  0.069  0.0080 0 0.00 0  42.0  5.56 0 0 0];

govdat = [
1  30  0.00850 0.080 0.000 0.150 0.050 10.000 0.280;
2  31  0.00600 0.180 0.030 0.200 0.000  8.000 0.300;
3  32  0.00600 0.180 0.030 0.200 0.000  8.000 0.300;
4  33  0.00600 0.180 0.030 0.200 0.000  8.000 0.300;
5  34  0.00600 0.180 0.030 0.200 0.000  8.000 0.300;
6  35  0.00600 0.180 0.030 0.200 0.000  8.000 0.300;
7  36  0.00600 0.180 0.030 0.200 0.000  8.000 0.300;
8  37  0.00600 0.180 0.030 0.200 0.000  8.000 0.300;
9  38  0.00548 0.100 0.000 0.200 0.100  8.720 0.300;
10 39  0.00548 0.100 0.000 0.200 0.100  8.720 0.300];

mpc.mac = macdat;
mpc.gov = govdat;
end

