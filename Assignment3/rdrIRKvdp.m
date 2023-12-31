function [r_vdp,dr_vdp] = rdrIRKvdp(deltaT,in2,in3,in4,t)
%rdrIRKvdp
%    [R_VDP,DR_VDP] = rdrIRKvdp(deltaT,IN2,IN3,IN4,T)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    18-Oct-2023 16:52:59

A1_1 = in4(1);
A1_2 = in4(3);
A2_1 = in4(2);
A2_2 = in4(4);
K1_1 = in3(1,:);
K1_2 = in3(3,:);
K2_1 = in3(2,:);
K2_2 = in3(4,:);
xk1 = in2(1,:);
xk2 = in2(2,:);
t2 = A1_1.*deltaT;
t3 = A1_2.*deltaT;
t4 = A2_1.*deltaT;
t5 = A2_2.*deltaT;
t14 = -xk1;
t6 = K1_1.*t2;
t7 = K1_2.*t3;
t8 = K2_1.*t2;
t9 = K1_1.*t4;
t10 = K2_2.*t3;
t11 = K1_2.*t5;
t12 = K2_1.*t4;
t13 = K2_2.*t5;
t15 = t6+t7+xk1;
t16 = t9+t11+xk1;
t17 = t8+t10+xk2;
t18 = t12+t13+xk2;
t19 = t15.^2;
t20 = t16.^2;
t21 = t19.*5.0;
t22 = t20.*5.0;
t23 = t21-5.0;
t24 = t22-5.0;
r_vdp = [-K1_1+t17;-K2_1-t6-t7+t14-t17.*t23;-K1_2+t18;-K2_2-t9-t11+t14-t18.*t24];
if nargout > 1
    dr_vdp = reshape([-1.0,-t2-t2.*t15.*t17.*1.0e+1,0.0,-t4-t4.*t16.*t18.*1.0e+1,t2,-t2.*t23-1.0,t4,-t4.*t24,0.0,-t3-t3.*t15.*t17.*1.0e+1,-1.0,-t5-t5.*t16.*t18.*1.0e+1,t3,-t3.*t23,t5,-t5.*t24-1.0],[4,4]);
end
