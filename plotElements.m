function ans = plotElements (a, ecc, RA ,Inclination,ArgPer, TrueA)

figure(3)
plot(ecc)
title('Eccentricity')

figure(4)
plot(Inclination)
title('Inclination')

figure(5)
plot(RA)
title('Right ascension')

figure(6)
plot(TrueA)
title('True Anomaly')

figure(7)
plot(a)
title('Semi-major axis')

figure(8)
plot(ArgPer)
title('Argument of Perigee')

ans = 1; 
end
