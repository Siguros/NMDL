clc;
clear;
close all;
exp_LTP_raw=[0	1.1356400E-04
1	1.1355160E-04
2	1.1399247E-04
3	1.1419967E-04
4	1.1436800E-04
5	1.1451153E-04
6	1.1466753E-04
7	1.1485113E-04
8	1.1499787E-04
9	1.1507953E-04
10	1.1522900E-04
11	1.1532987E-04
12	1.1541360E-04
13	1.1549493E-04
14	1.1550687E-04
15	1.1561653E-04
16	1.1573253E-04
17	1.1581720E-04
18	1.1584000E-04
19	1.1589533E-04
20	1.1593733E-04
21	1.1595893E-04
22	1.1601947E-04
23	1.1610133E-04
24	1.1613953E-04
25	1.1615900E-04
26	1.1616487E-04
27	1.1626233E-04
28	1.1626607E-04
29	1.1627693E-04
30	1.1631820E-04
31	1.1639687E-04
32	1.1645493E-04
33	1.1647007E-04
34	1.1654680E-04
35	1.1656127E-04
36	1.1657393E-04
37	1.1660653E-04
38	1.1664027E-04
39	1.1661727E-04
40	1.1668833E-04
41	1.1678160E-04
42	1.1677480E-04
43	1.1674993E-04
44	1.1680833E-04
45	1.1677340E-04
46	1.1683687E-04
47	1.1689200E-04
48	1.1698887E-04
49	1.1686767E-04
];

exp_LTD_raw=[ 0	1.135640E-04
1	1.136020E-04
2	1.136245E-04
3	1.136022E-04
4	1.136549E-04
5	1.136868E-04
6	1.137234E-04
7	1.137039E-04
8	1.136935E-04
9	1.137334E-04
10	1.137559E-04
11	1.137273E-04
12	1.137444E-04
13	1.138445E-04
14	1.138029E-04
15	1.138583E-04
16	1.138584E-04
17	1.139341E-04
18	1.139717E-04
19	1.139707E-04
20	1.140633E-04
21	1.140233E-04
22	1.140451E-04
23	1.140951E-04
24	1.141246E-04
25	1.142045E-04
26	1.141276E-04
27	1.142188E-04
28	1.142542E-04
29	1.142954E-04
30	1.143482E-04
31	1.143889E-04
32	1.143985E-04
33	1.144853E-04
34	1.145051E-04
35	1.145181E-04
36	1.145098E-04
37	1.146009E-04
38	1.146737E-04
39	1.146791E-04
40	1.147049E-04
41	1.148253E-04
42	1.148969E-04
43	1.150649E-04
44	1.150780E-04
45	1.152319E-04
46	1.153341E-04
47	1.155708E-04
48	1.159592E-04
49	1.168482E-04
];
xf_ltp = floor(max(exp_LTP_raw(:,1)));
xf_ltd = floor(max(exp_LTD_raw(:,1)));
x_step_ltp = 1/xf_ltp;
x_step_ltd = 1/xf_ltd;
yf_ltp = max(exp_LTP_raw(:,2));
yf_ltd = max(exp_LTD_raw(:,2));
yi_ltp = min(exp_LTP_raw(:,2));
yi_ltd = min(exp_LTD_raw(:,2));

exp_LTP(:,1) = exp_LTP_raw(:,1)/xf_ltp;
exp_LTD(:,1) = exp_LTD_raw(:,1)/xf_ltd;
exp_LTP(:,2) = (exp_LTP_raw(:,2)-yi_ltp)/(yf_ltp-yi_ltp);
exp_LTD(:,2) = (exp_LTD_raw(:,2)-yi_ltd)/(yf_ltd-yi_ltd);

plot(exp_LTP(:,1), exp_LTP(:,2), 'bo', 'LineWidth', 1);
hold on;
plot(exp_LTD(:,1), exp_LTD(:,2), 'ro', 'LineWidth', 1);

xf = 1;
A_LTP = 0.4;
B_LTP = 1./(1-exp(-1./A_LTP));
A_LTD = -0.2;
B_LTD = 1./(1-exp(-1./A_LTD));

% LTP fitting
var_amp = 0.02;    % LTP cycle-to-cycle variation
rng(103);
x_ltp(1) = 0;
y_ltp(1) = 0;
for n=1:1/x_step_ltp+1
    x_ltp(n+1) = x_ltp(n)+x_step_ltp;
    y_ltp(n+1) = B_LTP(1)*(1-exp(-x_ltp(n+1)/A_LTP(1)));
    delta_y = (y_ltp(n+1)-y_ltp(n)) + randn*var_amp;
    y_ltp(n+1) = y_ltp(n) + delta_y;   
    if y_ltp(n+1)>=1
        y_ltp(n+1)=1;
    elseif y_ltp(n+1)<=0
        y_ltp(n+1)=0;
    end
    x_ltp(n+1) = -A_LTP(1)*log(1-(y_ltp(n+1))/B_LTP(1));
end
plot((0:n-1)/(n-1), y_ltp(1:n), 'b', 'linewidth', 2);

% LTD fitting
var_amp = 0.035;    % LTD cycle-to-cycle variation
rng(898);
x_ltd(1) = 1;
y_ltd(1) = 1;
for n=1:1/x_step_ltd+1
    x_ltd(n+1) = x_ltd(n)-x_step_ltd;
    y_ltd(n+1) = B_LTD(1)*(1-exp(-x_ltd(n+1)/A_LTD(1)));
    delta_y = (y_ltd(n+1)-y_ltd(n)) + randn*var_amp;
    y_ltd(n+1) = y_ltd(n) + delta_y;
    if y_ltd(n+1)>=1
        y_ltd(n+1)=1;
    elseif y_ltd(n+1)<=0
        y_ltd(n+1)=0;
    end
    x_ltd(n+1) = -A_LTD(1)*log(1-(y_ltd(n+1))/B_LTD(1));
end
x_start = numel(x_ltd(:));
x_end = numel(x_ltd(:)) - n;
plot((n-1:-1:0)/(n-1), y_ltd(1:n), 'r', 'linewidth', 2);

xlabel('Normalized Pulse #');
ylabel('Normalized Conductance');
legend('Exp data (LTP)','Exp data (LTD)', 'Fit (LTP)', 'Fit (LTD)', 'location', 'southeast');