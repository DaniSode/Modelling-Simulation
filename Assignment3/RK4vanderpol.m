function [tRK4,xRK4] = RK4(tFinal, deltaT, x_0, BT)

a1 = BT(2,2);
a2 = BT(3,3);
a3 = BT(4,4);
c1 = BT(1,1);
c2 = BT(2,1);
c3 = BT(3,1);
c4 = BT(4,1);
b1 = BT(5,2);
b2 = BT(5,3);
b3 = BT(5,4);
b4 = BT(5,5);

nRK4 = tFinal / deltaT;

tRK4 = zeros(nRK4,1);
xRK4 = zeros(nRK4,2);
xRK4(1,:) = x_0;

for j = 1:nRK4
    tRK4(j+1) = tRK4(j) + deltaT;
    K1 = vanderpol(tRK4(j) + c1*deltaT, xRK4(j,:))';
    K2 = vanderpol(tRK4(j) + c2*deltaT, xRK4(j,:) +a1*deltaT*K1)';
    K3 = vanderpol(tRK4(j) + c3*deltaT, xRK4(j,:) +a2*deltaT*K2)';
    K4 = vanderpol(tRK4(j) + c4*deltaT, xRK4(j,:) +a3*deltaT*K3)';
    xRK4(j+1,:) = xRK4(j,:) + deltaT * b1 * K1 + deltaT * b2 * K2 ...
        + deltaT * b3 * K3 + deltaT * b4 * K4;
end

end