function fem = FEMUpt_4CSM_v5(fem,theta)
% 7 x 6 y 8 j
b_col = 20;  % width (along b_for)
h_col = 0;   % width (along h_for)
t_col = 2;   % thickness

pu = 0.5;
Ix = (h_col^3)*t_col/12+(t_col^3)*b_col/12;
Iy = ((b_col^3)*t_col/12+(t_col^3)*h_col/12)*pu;
Ci = 0.313;
J1  = Ci*(t_col^3)*b_col*2000;

fem.Elem_pro(1:4,7) = (theta(1)+1)*Ix;
fem.Elem_pro(5:8,7) = (theta(2)+1)*Ix;
fem.Elem_pro(9:12,7) = (theta(3)+1)*Ix;
fem.Elem_pro(13:16,7) = (theta(4)+1)*Ix;
fem.Elem_pro(17:20,7) = (theta(5)+1)*Ix;
fem.Elem_pro(21:24,7) = (theta(6)+1)*Ix;
fem.Elem_pro(25:28,7) = (theta(7)+1)*Ix;
fem.Elem_pro(29:32,7) = (theta(8)+1)*Ix;

fem.Elem_pro(1:4,6) = (theta(9)+1)*Iy;
fem.Elem_pro(5:8,6) = (theta(10)+1)*Iy;
fem.Elem_pro(9:12,6) = (theta(11)+1)*Iy;
fem.Elem_pro(13:16,6) = (theta(12)+1)*Iy;
fem.Elem_pro(17:20,6) = (theta(13)+1)*Iy;
fem.Elem_pro(21:24,6) = (theta(14)+1)*Iy;
fem.Elem_pro(25:28,6) = (theta(15)+1)*Iy;
fem.Elem_pro(29:32,6) = (theta(16)+1)*Iy;
% % 
%  % fem.Elem_pro(1:32,8) = (theta(17)+1)*J1;
% 
fem.Elem_pro(1:4,8) = (theta(17)+1)*J1;
fem.Elem_pro(5:8,8) = (theta(18)+1)*J1;
fem.Elem_pro(9:12,8) = (theta(19)+1)*J1;
fem.Elem_pro(13:16,8) = (theta(20)+1)*J1;
fem.Elem_pro(17:20,8) = (theta(21)+1)*J1;
fem.Elem_pro(21:24,8) = (theta(22)+1)*J1;
fem.Elem_pro(25:28,8) = (theta(23)+1)*J1;
fem.Elem_pro(29:32,8) = (theta(24)+1)*J1;
%%
% fem.Elem_pro(1:4,8) = (theta(17)+1)*J1;
% fem.Elem_pro(5:8,8) = (theta(17)+1)*J1;
% fem.Elem_pro(9:12,8) = (theta(18)+1)*J1;
% fem.Elem_pro(13:16,8) = (theta(18)+1)*J1;
% fem.Elem_pro(17:20,8) = (theta(19)+1)*J1;
% fem.Elem_pro(21:24,8) = (theta(19)+1)*J1;
% fem.Elem_pro(25:28,8) = (theta(20)+1)*J1;
% fem.Elem_pro(29:32,8) = (theta(20)+1)*J1;

%%
% fem.Elem_pro(1:4,7) = theta(1);
% fem.Elem_pro(5:8,7) = theta(2);
% fem.Elem_pro(9:12,7) = theta(3);
% fem.Elem_pro(13:16,7) = theta(4);
% fem.Elem_pro(17:20,7) = theta(5);
% fem.Elem_pro(21:24,7) = theta(6);
% fem.Elem_pro(25:28,7) = theta(7);
% fem.Elem_pro(29:32,7) = theta(8);

% fem.Elem_pro(1:32,6) = theta(9);
% fem.Elem_pro(1:32,8) = theta(10);
% 
% fem.Elem_pro(1:4,6) = theta(9);
% fem.Elem_pro(5:8,6) = theta(9);
% fem.Elem_pro(9:12,6) = theta(10);
% fem.Elem_pro(13:16,6) = theta(10);
% fem.Elem_pro(17:20,6) = theta(11);
% fem.Elem_pro(21:24,6) = theta(11);
% fem.Elem_pro(25:28,6) = theta(12);
% fem.Elem_pro(29:32,6) = theta(12);

% fem.Elem_pro(1:32,8) = theta(13);


% fem.Elem_pro(1:4,8) = theta(13);
% fem.Elem_pro(5:8,8) = theta(13);
% fem.Elem_pro(9:12,8) = theta(14);
% fem.Elem_pro(13:16,8) = theta(14);
% fem.Elem_pro(17:20,8) = theta(15);
% fem.Elem_pro(21:24,8) = theta(15);
% fem.Elem_pro(25:28,8) = theta(16);
% fem.Elem_pro(29:32,8) = theta(16);

% fem.Elem_pro(1:4,6) = theta(9);
% fem.Elem_pro(5:8,6) = theta(10);
% fem.Elem_pro(9:12,6) = theta(11);
% fem.Elem_pro(13:16,6) = theta(12);
% fem.Elem_pro(17:20,6) = theta(13);
% fem.Elem_pro(21:24,6) = theta(14);
% fem.Elem_pro(25:28,6) = theta(15);
% fem.Elem_pro(29:32,6) = theta(16);

%fem.Elem_pro(1:32,8) = theta(17);

% fem.Elem_pro(1:4,8) = theta(17);
% fem.Elem_pro(5:8,8) = theta(17);
% fem.Elem_pro(9:12,8) = theta(18);
% fem.Elem_pro(13:16,8) = theta(18);
% fem.Elem_pro(17:20,8) = theta(19);
% fem.Elem_pro(21:24,8) = theta(19);
% fem.Elem_pro(25:28,8) = theta(20);
% fem.Elem_pro(29:32,8) = theta(20);

% fem.Elem_pro(1:4,8) = theta(17);
% fem.Elem_pro(5:8,8) = theta(18);
% fem.Elem_pro(9:12,8) = theta(19);
% fem.Elem_pro(13:16,8) = theta(20);
% fem.Elem_pro(17:20,8) = theta(21);
% fem.Elem_pro(21:24,8) = theta(22);
% fem.Elem_pro(25:28,8) = theta(23);
% fem.Elem_pro(29:32,8) = theta(24);

% elseif a == 2
% fem.Elem_pro(1:4,6) = theta(1);
% fem.Elem_pro(5:8,6) = theta(2);
% fem.Elem_pro(9:12,6) = theta(3);
% fem.Elem_pro(13:16,6) = theta(4);
% fem.Elem_pro(17:20,6) = theta(5);
% fem.Elem_pro(21:24,6) = theta(6);
% fem.Elem_pro(25:28,6) = theta(7);
% fem.Elem_pro(29:32,6) = theta(8);

% elseif a == 3
% fem.Elem_pro(1:2,7) = theta(1);
% fem.Elem_pro(3:4,7) = theta(2);
% fem.Elem_pro(5:6,7) = theta(3);
% fem.Elem_pro(7:8,7) = theta(4);
% fem.Elem_pro(9:10,7) = theta(5);
% fem.Elem_pro(11:12,7) = theta(6);
% fem.Elem_pro(13:14,7) = theta(7);
% fem.Elem_pro(15:16,7) = theta(8);
% fem.Elem_pro(17:18,7) = theta(9);
% fem.Elem_pro(19:20,7) = theta(10);
% fem.Elem_pro(21:22,7) = theta(11);
% fem.Elem_pro(23:24,7) = theta(12);
% fem.Elem_pro(25:26,7) = theta(13);
% fem.Elem_pro(27:28,7) = theta(14);
% fem.Elem_pro(29:30,7) = theta(15);
% fem.Elem_pro(31:32,7) = theta(16);

% % elseif
% fem.Elem_pro([1,3],7) = theta(1);
% fem.Elem_pro([2,4],7) = theta(2);
% fem.Elem_pro([5,7],7) = theta(3);
% fem.Elem_pro([6,8],7) = theta(4);
% fem.Elem_pro([9,11],7) = theta(5);
% fem.Elem_pro([10,12],7) = theta(6);
% fem.Elem_pro([13,15],7) = theta(7);
% fem.Elem_pro([14,16],7) = theta(8);
% fem.Elem_pro([17,19],7) = theta(9);
% fem.Elem_pro([18,20],7) = theta(10);
% fem.Elem_pro([21,23],7) = theta(11);
% fem.Elem_pro([22,24],7) = theta(12);
% fem.Elem_pro([25,27],7) = theta(13);
% fem.Elem_pro([26,28],7) = theta(14);
% fem.Elem_pro([29,31],7) = theta(15);
% fem.Elem_pro([30,32],7) = theta(16);

%%
% fem.Elem_pro(1,7) = theta(1);
% fem.Elem_pro(3,7) = theta(1);
% fem.Elem_pro(2,7) = theta(2);
% fem.Elem_pro(4,7) = theta(2);
% fem.Elem_pro(5,7) = theta(3);
% fem.Elem_pro(7,7) = theta(3);
% fem.Elem_pro(6,7) = theta(4);
% fem.Elem_pro(8,7) = theta(4);
% fem.Elem_pro(9,7) = theta(5);
% fem.Elem_pro(11,7) = theta(5);
% fem.Elem_pro(10,7) = theta(6);
% fem.Elem_pro(12,7) = theta(6);
% fem.Elem_pro(13,7) = theta(7);
% fem.Elem_pro(15,7) = theta(7);
% fem.Elem_pro(14,7) = theta(8);
% fem.Elem_pro(16,7) = theta(8);
% fem.Elem_pro(17,7) = theta(9);
% fem.Elem_pro(19,7) = theta(9);
% fem.Elem_pro(18,7) = theta(10);
% fem.Elem_pro(20,7) = theta(10);
% fem.Elem_pro(21,7) = theta(11);
% fem.Elem_pro(23,7) = theta(11);
% fem.Elem_pro(22,7) = theta(12);
% fem.Elem_pro(24,7) = theta(12);
% fem.Elem_pro(25,7) = theta(13);
% fem.Elem_pro(27,7) = theta(13);
% fem.Elem_pro(26,7) = theta(14);
% fem.Elem_pro(28,7) = theta(14);
% fem.Elem_pro(29,7) = theta(15);
% fem.Elem_pro(31,7) = theta(15);
% fem.Elem_pro(30,7) = theta(16);
% fem.Elem_pro(32,7) = theta(16);

%%
% fem.Elem_pro(1,7) = theta(1);
% fem.Elem_pro(3,7) = theta(2);
% fem.Elem_pro(2,7) = theta(2);
% fem.Elem_pro(4,7) = theta(1);
% fem.Elem_pro(5,7) = theta(3);
% fem.Elem_pro(7,7) = theta(4);
% fem.Elem_pro(6,7) = theta(4);
% fem.Elem_pro(8,7) = theta(3);
% fem.Elem_pro(9,7) = theta(5);
% fem.Elem_pro(11,7) = theta(6);
% fem.Elem_pro(10,7) = theta(6);
% fem.Elem_pro(12,7) = theta(5);
% fem.Elem_pro(13,7) = theta(7);
% fem.Elem_pro(15,7) = theta(8);
% fem.Elem_pro(14,7) = theta(8);
% fem.Elem_pro(16,7) = theta(7);
% fem.Elem_pro(17,7) = theta(9);
% fem.Elem_pro(19,7) = theta(10);
% fem.Elem_pro(18,7) = theta(10);
% fem.Elem_pro(20,7) = theta(9);
% fem.Elem_pro(21,7) = theta(11);
% fem.Elem_pro(23,7) = theta(12);
% fem.Elem_pro(22,7) = theta(12);
% fem.Elem_pro(24,7) = theta(11);
% fem.Elem_pro(25,7) = theta(13);
% fem.Elem_pro(27,7) = theta(14);
% fem.Elem_pro(26,7) = theta(14);
% fem.Elem_pro(28,7) = theta(13);
% fem.Elem_pro(29,7) = theta(15);
% fem.Elem_pro(31,7) = theta(16);
% fem.Elem_pro(30,7) = theta(16);
% fem.Elem_pro(32,7) = theta(15);

%%
% fem.Elem_pro(1:2,7) = theta(1);
% fem.Elem_pro(3:4,7) = theta(2);
% fem.Elem_pro(5:6,7) = theta(3);
% fem.Elem_pro(7:8,7) = theta(4);
% fem.Elem_pro(9:12,7) = theta(5);
% fem.Elem_pro(13:16,7) = theta(6);
% fem.Elem_pro(17:20,7) = theta(7);
% fem.Elem_pro(21:24,7) = theta(8);
% fem.Elem_pro(25:28,7) = theta(9);
% fem.Elem_pro(29:32,7) = theta(10);

%%
% fem.Elem_pro(1,7) = theta(1);
% fem.Elem_pro(3,7) = theta(1);
% fem.Elem_pro(2,7) = theta(2);
% fem.Elem_pro(4,7) = theta(2);
% fem.Elem_pro(5:8,7) = theta(3);
% fem.Elem_pro(9:12,7) = theta(4);
% fem.Elem_pro(13:16,7) = theta(5);
% fem.Elem_pro(17:20,7) = theta(6);
% fem.Elem_pro(21:24,7) = theta(7);
% fem.Elem_pro(25:28,7) = theta(8);
% fem.Elem_pro(29:30,7) = theta(9);
% fem.Elem_pro(31:32,7) = theta(10);

% %%
% fem.Elem_pro(1:4,7) = theta(1);
% fem.Elem_pro(5:6,7) = theta(2);
% fem.Elem_pro(7:8,7) = theta(3);
% fem.Elem_pro(9:12,7) = theta(4);
% fem.Elem_pro(13:16,7) = theta(5);
% fem.Elem_pro(17:20,7) = theta(6);
% fem.Elem_pro(21:22,7) = theta(7);
% fem.Elem_pro(23:24,7) = theta(8);
% fem.Elem_pro(25:28,7) = theta(9);
% fem.Elem_pro(29:32,7) = theta(10);

% end
% N_col = 32;
% fem.Ele_pro(1:N_col,7) = theta(1:N_col);
% fem.Ele_pro(1:N_col,7) = theta(N_col+1:2*N_col);
end