function Y = plot_thrust (V,n)

for i = 1:1:n
    
   R_r(i) = V(i,1);
   
   R_theta(i) = V(i,2);
   
   R_normal(i) = V(i,3);
   
end

figure (10)
plot(R_r)
title('Radial Thrust')

figure (11)
plot(R_theta)
title('Theta Thrust')

figure (12)
plot(R_normal)
title('Normal to plane Thrust')

Y = 0;

end