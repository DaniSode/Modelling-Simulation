function [tRK2,xRK2] = RK2(tFinal, deltaT, x_0, lambda,BT)

a= BT(2,2);
c1= BT(1,1);
c2 = BT(2,1); 
b1 = BT(3,2);
b2 = BT(3,3);

nRK2 = tFinal / deltaT;

tRK2 = zeros(nRK2,1);
xRK2 = zeros(nRK2,1);
xRK2(1) = x_0;

for j = 1:nRK2
    tRK2(j+1) = tRK2(j) + deltaT;
    K1 = testfunction(xRK2(j), lambda + c1*deltaT);
    K2 = testfunction(xRK2(j)+a*deltaT*K1, lambda + c2*deltaT);
    xRK2(j+1) = xRK2(j) + deltaT * b1 * K1 + deltaT * b2 * K2 ;
end

end