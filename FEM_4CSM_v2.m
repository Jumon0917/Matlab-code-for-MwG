function fem = FEM_4CSM_v2
%% FEM model construction
% -----------------Input:parameter of the model-------------------
% number of floors
Nfloor = 8;
% hight of the model
height = 1200;
% size of the floor
b_flr = 250; % width
h_flr = 350; % length
t_flr = 20;  % thickness
% size of the column
b_col = 20;  % width (along b_for)
h_col = 0;   % width (along h_for)
t_col = 2;   % thickness
% number of elements per column
Nele_col = 1;
% --------------------------------------------------------------
% number of layers (each layer has 4 columns)
Nlayer = Nfloor*Nele_col;
% length of the column elements
Len_colele = height/Nlayer;
% number of nodes
Nnode = Nlayer*4+4;
% number of elements
Nele = Nlayer*4+Nfloor*4;
% coordinates of nodes
Node = zeros(Nnode,6);
% the nodes index in elements
Ele = zeros(Nele,2);
% element property
Ele_pro = zeros(Nele,9);
% -----------------------Node coordination-------------------------
for isle = 1:Nlayer+1
    Node((isle-1)*4+1,:) = [0,0,(isle-1)*Len_colele,0,0,0];
    Node((isle-1)*4+2,:) = [b_flr,0,(isle-1)*Len_colele,0,0,0];
    Node((isle-1)*4+3,:) = [0,h_flr,(isle-1)*Len_colele,0,0,0];
    Node((isle-1)*4+4,:) = [b_flr,h_flr,(isle-1)*Len_colele,0,0,0];
end
% -------------------Elements section settings---------------------
% matrial property
rho = 2.7e-9; % density
E   = 7.0e4;  % young's module
mu  = 0.3;    % poission ratio
% section geometry
% column
A1  = (b_col+h_col)*t_col;
Iy1 = ((b_col^3)*t_col/12+(t_col^3)*h_col/12)*0.5;% strong
Iz1 = (h_col^3)*t_col/12+(t_col^3)*b_col/12;
Ci = 0.313;
J1  = Ci*(t_col^3)*b_col*2000;
T1 = 0;
% floor
A2 = t_flr*b_flr/4; Iy2 = 1e10; Iz2 = 1e10; J2 = 1e10;
A3 = t_flr*h_flr/4; Iy3 = 1e10; Iz3 = 1e10; J3 = 1e10;
T2 = 0; T3 = 0;
% -----------------------Nodes of elements-------------------------
% columns
for isle = 1:Nlayer
    Ele((isle-1)*4+1,:) = [(isle-1)*4+1,(isle-1)*4+5];
    Ele((isle-1)*4+2,:) = [(isle-1)*4+2,(isle-1)*4+6];
    Ele((isle-1)*4+3,:) = [(isle-1)*4+3,(isle-1)*4+7];
    Ele((isle-1)*4+4,:) = [(isle-1)*4+4,(isle-1)*4+8];
end
% rigid plates (floors)
for iflr = 1:Nfloor
    Ele(Nlayer*4+(iflr-1)*4+1,:) = [iflr*Nele_col*4+1,iflr*Nele_col*4+2];
    Ele(Nlayer*4+(iflr-1)*4+2,:) = [iflr*Nele_col*4+2,iflr*Nele_col*4+4];
    Ele(Nlayer*4+(iflr-1)*4+3,:) = [iflr*Nele_col*4+4,iflr*Nele_col*4+3];
    Ele(Nlayer*4+(iflr-1)*4+4,:) = [iflr*Nele_col*4+3,iflr*Nele_col*4+1];
end
% ----------------------Element property--------------------------
% columns
for isle = 1:Nlayer
    Ele_pro((isle-1)*4+1,:) = [6,rho,A1,E,mu,Iy1,Iz1,J1,T1];
    Ele_pro((isle-1)*4+2,:) = [6,rho,A1,E,mu,Iy1,Iz1,J1,T1];
    Ele_pro((isle-1)*4+3,:) = [6,rho,A1,E,mu,Iy1,Iz1,J1,T1];
    Ele_pro((isle-1)*4+4,:) = [6,rho,A1,E,mu,Iy1,Iz1,J1,T1];
end
% rigid plates (floors)
for iflr = 1:Nfloor
    Ele_pro(Nlayer*4+(iflr-1)*4+3,:) = [6,rho,A3,E,mu,Iy3,Iz3,J3,T3];
    Ele_pro(Nlayer*4+(iflr-1)*4+1,:) = [6,rho,A3,E,mu,Iy3,Iz3,J3,T3];
    Ele_pro(Nlayer*4+(iflr-1)*4+4,:) = [6,rho,A2,E,mu,Iy2,Iz2,J2,T2];
    Ele_pro(Nlayer*4+(iflr-1)*4+2,:) = [6,rho,A2,E,mu,Iy2,Iz2,J2,T2];
end
% ---------------------Boundary conditions-------------------------
BC = [ 1   1   0
    1   2   0
    1   3   0
    1   4   0
    1   5   0
    1   6   0
    2   1   0
    2   2   0
    2   3   0
    2   4   0
    2   5   0
    2   6   0
    3   1   0
    3   2   0
    3   3   0
    3   4   0
    3   5   0
    3   6   0
    4   1   0
    4   2   0
    4   3   0
    4   4   0
    4   5   0
    4   6   0];
fem = FEM(Node,Ele,Ele_pro,BC);
return