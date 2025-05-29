% ===================================================== %
%                                                       %
%           Supplementary Material to:                  %
%                                                       %
%   Effective spin Hamiltonians for the quantum-rotor   %
%   tunneling problem in pulse EPR                      %
%                                                       %
%           By: S.Kuzin and G.Jeschke                   %
%                                                       %
%   This MATLAB script generates effective tunneling    %
%   Hamiltonians ht___ in four cases:                   %
%       (1) 'CH3' rotor                                 %
%       (2) 'CD3' rotor                                 %
%       (3) 'CH4' rotor                                 %
%       (4) 'CD4' rotor                                 %
%                                                       %
%   Please pay attention to the tunneling frequency     %
%                                                       %
%   Requires installation of the EasySpin package       %
%   See: https://www.easyspin.org/                      %
%                                                       %
%   2025                                                %
%                                                       %
% ===================================================== %

%% CH3 = protonated C3 rotor

[x1, y1, z1] = sop([1/2, 1/2, 1/2], 'xee', 'yee', 'zee');
[x2, y2, z2] = sop([1/2, 1/2, 1/2], 'exe', 'eye', 'eze');
[x3, y3, z3] = sop([1/2, 1/2, 1/2], 'eex', 'eey', 'eez');

Isym = x1*x2 + y1*y2 + z1*z2 + x2*x3 + y2*y3 + z2*z3 + x3*x1 + y3*y1 + z3*z1 + 1/4*eye(8);
Iasym = 2*( (x1*y2 - y1*x2)*z3 + (x2*y3 - y2*x3)*z1 + (x3*y1 - y3*x1)*z2 );

IA = Isym;
IE1 = cos(2*pi/3)*(Isym) + sin(2*pi/3)*(Iasym);
IE2 = cos(4*pi/3)*(Isym) + sin(4*pi/3)*(Iasym);

% Effective tunneling Hamiltonians

wt = 1; % tunneling frequency

htA =  -2*wt/3*Isym;
htE1 = -2*wt/3*(cos(2*pi/3)*Isym + sin(2*pi/3)*Iasym);
htE2 = -2*wt/3*(cos(4*pi/3)*Isym + sin(4*pi/3)*Iasym);

%% CD3 = deuterated C3 rotor

[x1, y1, z1] = sop([1, 1, 1], 'xee', 'yee', 'zee');
[x2, y2, z2] = sop([1, 1, 1], 'exe', 'eye', 'eze');
[x3, y3, z3] = sop([1, 1, 1], 'eex', 'eey', 'eez');

% Pair-exchange operators
I12 = x1*x2 + y1*y2 + z1*z2;
I23 = x2*x3 + y2*y3 + z2*z3;
I31 = x3*x1 + y3*y1 + z3*z1;

% Creating auxiliary operators L
L0 = eye(27);
L1 = 1/3*(I12 + I23 + I31);
L2a = 1/3*(I12^2 + I23^2 + I31^2);
L2b = 1/3*((I12*I23 + I23*I12) + (I23*I31 + I31*I23) + (I31*I12 + I12*I31));
L3 = 1/3*(I12*I23*I31 + I23*I31*I12 + I31*I12*I23 + I31*I23*I12 + I12*I31*I23 + I23*I12*I31);
Lm2 = 1i/3*((I12*I23 - I23*I12) + (I23*I31 - I31*I23) + (I31*I12 - I12*I31));
Lm3 = 1i/3*(I12*I23*I31 + I23*I31*I12 + I31*I12*I23 - I31*I23*I12 - I12*I31*I23 - I23*I12*I31);

Isym = -3*L0 + L1 + 3*L2a + 0.5*L2b - 0.5*L3;
Iasym = 0.5*Lm2 + 1.5*Lm3;

% Effective tunneling Hamiltonians

wt = 1; % tunneling frequency

htA = -2*wt/3*Isym;
htE1 = -2*wt/3*(cos(2*pi/3)*Isym + sin(2*pi/3)*Iasym);
htE2 = -2*wt/3*(cos(4*pi/3)*Isym + sin(4*pi/3)*Iasym);

%% CH4 = protonated T rotor

[x1, y1, z1] = sop([1/2, 1/2, 1/2, 1/2], 'xeee', 'yeee', 'zeee');
[x2, y2, z2] = sop([1/2, 1/2, 1/2, 1/2], 'exee', 'eyee', 'ezee');
[x3, y3, z3] = sop([1/2, 1/2, 1/2, 1/2], 'eexe', 'eeye', 'eeze');
[x4, y4, z4] = sop([1/2, 1/2, 1/2, 1/2], 'eeex', 'eeey', 'eeez');

% Pair-exchange operators
I12 = x1*x2 + y1*y2 + z1*z2;
I13 = x1*x3 + y1*y3 + z1*z3;
I14 = x1*x4 + y1*y4 + z1*z4;
I23 = x2*x3 + y2*y3 + z2*z3;
I24 = x2*x4 + y2*y4 + z2*z4;
I34 = x3*x4 + y3*y4 + z3*z4;

Isym = 4*(I12+I13+I14+I23+I24+I34) + 2*eye(16);
C12_34 = (eye(16) + 4*I12)*(eye(16) + 4*I34)/4;
C13_24 = (eye(16) + 4*I13)*(eye(16) + 4*I24)/4;
C14_23 = (eye(16) + 4*I14)*(eye(16) + 4*I23)/4;
IV0 =  C12_34 + C13_24 + C14_23;
IV1 = -C12_34 - C13_24 + C14_23;
IV2 = -C12_34 + C13_24 - C14_23;
IV3 =  C12_34 - C13_24 - C14_23;
Iasym = 4*(x1*(y2*z3-z2*y3) - x2*(y1*z3-z1*y3) + x3*(y1*z2-z1*y2) ...
    - x1*(y2*z4-z2*y4) + x2*(y1*z4-z1*y4) - x4*(y1*z2-z1*y2) ...
    + x1*(y3*z4-z3*y4) - x3*(y1*z4-z1*y4) + x4*(y1*z3-z1*y3) ...
    - x2*(y3*z4-z3*y4) + x3*(y2*z4-z2*y4) - x4*(y2*z3-z2*y3));

% Effective tunneling Hamiltonians

a = -1; % see the meaning of a and b in the original paper
b = 1;

htA = a*Isym + b*IV0;
htE1 = a*(cos(2*pi/3)*Isym + sin(2*pi/3)*Iasym) + b*IV0;
htE2 = a*(cos(4*pi/3)*Isym + sin(4*pi/3)*Iasym) + b*IV0;
htT1_A = 2*a*(eye(16) + IV1) + b*IV1;
htT1_E1 = -a*(eye(16) + IV1) + b*IV1;
htT1_E2 = -a*(eye(16) + IV1) + b*IV1;
htT2_A = 2*a*(eye(16) + IV2) + b*IV2;
htT2_E1 = -a*(eye(16) + IV2) + b*IV2;
htT2_E2 = -a*(eye(16) + IV2) + b*IV2;
htT3_A = 2*a*(eye(16) + IV3) + b*IV3;
htT3_E1 = -a*(eye(16) + IV3) + b*IV3;
htT3_E2 = -a*(eye(16) + IV3) + b*IV3;

%% CD4 = deuterated T rotor

[x1, y1, z1] = sop([1, 1, 1, 1], 'xeee', 'yeee', 'zeee');
[x2, y2, z2] = sop([1, 1, 1, 1], 'exee', 'eyee', 'ezee');
[x3, y3, z3] = sop([1, 1, 1, 1], 'eexe', 'eeye', 'eeze');
[x4, y4, z4] = sop([1, 1, 1, 1], 'eeex', 'eeey', 'eeez');

xs = {x1, x2, x3, x4};
ys = {y1, y2, y3, y4};
zs = {z1, z2, z3, z4};

% Pair-exchange operators
I12 = x1*x2 + y1*y2 + z1*z2;
I13 = x1*x3 + y1*y3 + z1*z3;
I14 = x1*x4 + y1*y4 + z1*z4;
I23 = x2*x3 + y2*y3 + z2*z3;
I24 = x2*x4 + y2*y4 + z2*z4;
I34 = x3*x4 + y3*y4 + z3*z4;

C12_34 = (eye(81) + I12 + I12^2)*(eye(81) + I34 + I34^2);
C13_24 = (eye(81) + I13 + I13^2)*(eye(81) + I24 + I24^2);
C14_23 = (eye(81) + I14 + I14^2)*(eye(81) + I23 + I23^2);
IV0 =  C12_34 + C13_24 + C14_23;
IV1 = -C12_34 - C13_24 + C14_23;
IV2 = -C12_34 + C13_24 - C14_23;
IV3 =  C12_34 - C13_24 - C14_23;

[L0_123, L1_123, L2a_123, L2b_123, L3_123, Lm2_123, Lm3_123] = getLOps(xs, ys, zs, 1, 2, 3);
[L0_142, L1_142, L2a_142, L2b_142, L3_142, Lm2_142, Lm3_142] = getLOps(xs, ys, zs, 1, 4, 2);
[L0_134, L1_134, L2a_134, L2b_134, L3_134, Lm2_134, Lm3_134] = getLOps(xs, ys, zs, 1, 3, 4);
[L0_243, L1_243, L2a_243, L2b_243, L3_243, Lm2_243, Lm3_243] = getLOps(xs, ys, zs, 2, 4, 3);

L0 =   L0_123 +  L0_142 +  L0_134 +  L0_243;
L1 =   L1_123 +  L1_142 +  L1_134 +  L1_243;
L2a = L2a_123 + L2a_142 + L2a_134 + L2a_243;
L2b = L2b_123 + L2b_142 + L2b_134 + L2b_243;
L3 =   L3_123 +  L3_142 +  L3_134 +  L3_243;
Lm2 = Lm2_123 + Lm2_142 + Lm2_134 + Lm2_243;
Lm3 = Lm3_123 + Lm3_142 + Lm3_134 + Lm3_243;

Isym = -3*L0 + L1 + 3*L2a + 0.5*L2b - 0.5*L3;
Iasym = 0.5*Lm2 + 1.5*Lm3;

a = -1; % see the meaning of a and b in the original paper
b = 1;

htA = a*Isym + b*IV0;
htE1 = a*(cos(2*pi/3)*Isym + sin(2*pi/3)*Iasym) + b*IV0;
htE2 = a*(cos(4*pi/3)*Isym + sin(4*pi/3)*Iasym) + b*IV0;
htT1_A = 2*a*(eye(81) + IV1) + b*IV1;
htT1_E1 = -a*(eye(81) + IV1) + b*IV1;
htT1_E2 = -a*(eye(81) + IV1) + b*IV1;
htT2_A = 2*a*(eye(81) + IV2) + b*IV2;
htT2_E1 = -a*(eye(81) + IV2) + b*IV2;
htT2_E2 = -a*(eye(81) + IV2) + b*IV2;
htT3_A = 2*a*(eye(81) + IV3) + b*IV3;
htT3_E1 = -a*(eye(81) + IV3) + b*IV3;
htT3_E2 = -a*(eye(81) + IV3) + b*IV3;

function [L0, L1, L2a, L2b, L3, Lm2, Lm3] = getLOps(xs, ys, zs, i1, i2, i3)
    % get all L operators with substitution 1 -> i1, 2 -> i2, 3 -> i3
    % xs, ys and zs contain Cartesian spin operators
    % e.g. xs = {x1, x2, x3, x4}
    I12 = xs{i1}*xs{i2} + ys{i1}*ys{i2} + zs{i1}*zs{i2};
    I31 = xs{i1}*xs{i3} + ys{i1}*ys{i3} + zs{i1}*zs{i3};
    I23 = xs{i2}*xs{i3} + ys{i2}*ys{i3} + zs{i2}*zs{i3};

    L0 = eye(81);
    L1 = 1/3*(I12 + I23 + I31);
    L2a = 1/3*(I12^2 + I23^2 + I31^2);
    L2b = 1/3*((I12*I23 + I23*I12) + (I23*I31 + I31*I23) + (I31*I12 + I12*I31));
    L3 = 1/3*(I12*I23*I31 + I23*I31*I12 + I31*I12*I23 + I31*I23*I12 + I12*I31*I23 + I23*I12*I31);
    Lm2 = 1i/3*((I12*I23 - I23*I12) + (I23*I31 - I31*I23) + (I31*I12 - I12*I31));
    Lm3 = 1i/3*(I12*I23*I31 + I23*I31*I12 + I31*I12*I23 - I31*I23*I12 - I12*I31*I23 - I23*I12*I31);
end