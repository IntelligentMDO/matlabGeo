function y_Camber = Catia_IGP_Camber(c1, c2, c3, c4, x0)

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


%% 弯度函数
Control_X = [c1 c1+(1-c1)*c2];
Control_Y = [c3 c4];

syms k X;
finv = solve((3*k*(1-k)^2*Control_X(1) + 3*(1-k)*k*k*Control_X(2) + k^3) - X,k);
finv = real(double( subs(finv,X,x0) ));
if finv(1,2)>=0 && finv(1,2)<1
    K = finv(1,:);
elseif finv(2,2)>=0 && finv(2,2)<1
    K = finv(2,:);
elseif finv(3,2)>=0 && finv(3,2)<1
    K = finv(3,:);
end
for n = 3:size(K,2)
    if sign( K(1,n-1)-K(1,n-2) ) ~= sign( K(1,n)-K(1,n-1) ) && sign( K(1,n-1)-K(1,n-2) ) ~= 0
        if sign( finv(1,n-1)-finv(1,n-2) ) == sign( finv(2,n)-finv(2,n-1) )
            K(n:end) = finv(2,n:end);
        else
            K(n:end) = finv(3,n:end);
        end
        break;
    end
end
fy = @(my,K)(3*K.*(1-K).^2*my(1) + 3*(1-K).*K.*K*my(2));
y_Camber = fy(Control_Y,K);


end










