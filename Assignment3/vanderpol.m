function dydt = vanderpol(t,y)

u=5;
dydt = [y(2); 
        u*(1-y(1)^2)*y(2)-y(1)];
end