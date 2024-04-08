function [y_Up, y_Low] = Catia_IGP_Thickness(X_T, T, Rho0, Beta_TE, z_TE, x0, y_Camber)

%********************************************************************
%  c1,c2：coefficients of camber-line-abscissa parameter equation
%  c3,c4：coefficients of camber-line-ordinate parameter equation
%      T：maximum thickness
%    X_T：chordwise location of maximum thickness
%   Rho0：relative quantity of leading edge radius
%Beta_TE：relative quantity of trailing edge boat-tail angle
%     x0：abscissa of airfoil control points
%
%|取值范围|
%     c1：[0.010, 0.960]
%     c2：[0.020, 0.970]
%     c3：[-0.074, 0.247]
%     c4：[-0.102, 0.206]
% 
%    X_T：[0.2002, 0.4813]
%      T：[0.0246, 0.3227]
%   rho0：[0.1750, 1.4944]
%Beta_TE：[0.1452, 4.8724]
%********************************************************************


%% 厚度函数
Rho0 = Rho0*T*T/X_T/X_T;
Beta_TE = Beta_TE*atan(T/(1-X_T));

syms x;
Pt = [T, 0, -tan(Beta_TE/2), sqrt(2*Rho0), z_TE];
g0 = [
     x^0.5
     x
     x^2
     x^3
     x^4
     ];
dg0 = [
      0.5/x^0.5
      1
      2*x
      3*x^2
      4*x^3
      ];
k0 = [
     1
     0
     0
     0
     0
     ];
G = [subs(g0,x,X_T), subs(dg0,x,X_T), subs(dg0,x,1)/2, k0, subs(g0,x,1)];
G = double(G);
g = @(t,K)(t(1)*K.^0.5 + t(2)*K + t(3)*K.^2 + t(4)*K.^3 + t(5)*K.^4);
t = Pt/G;


%% 整合 
y_Up = y_Camber + g(t,x0)/2;
y_Low = y_Camber - g(t,x0)/2;


end










