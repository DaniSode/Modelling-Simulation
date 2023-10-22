function [tRK1,xRK1] = RK1(tFinal, deltaT, x_0, lambda,BT)

c1= BT(1,1); 
b1 = BT(2,2);

nRK1 = tFinal / deltaT;
tRK1 = zeros(nRK1,1);
xRK1 = zeros(nRK1,1);
xRK1(1) = x_0;

for j = 1:nRK1
    tRK1(j+1) = tRK1(j) + deltaT;
    K1 = testfunction(xRK1(j), lambda);
    xRK1(j+1) = xRK1(j) + deltaT * b1 * K1 ;
end

end