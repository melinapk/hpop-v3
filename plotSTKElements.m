function [ DifferenceEcc , DifferenceIncdeg ,DifferenceRAANdeg ,DifferenceTrueAnomaly ,DifferenceSMA ,DifferenceAofPer ] = plotSTKElements (a , ecc , inc , raan , true , sma , argumentofperigee )

figure(9)
plot(a.Eccentricity)
title('Eccentricity STK')

figure(10)
plot(a.Inclinationdeg)
title('Inclination STK')

figure(11)
plot(a.RAANdeg)
title('Right ascension  STK')

figure(12)
plot(a.TrueAnomalydeg)
title('True Anomaly STK')

figure(13)
plot(a.SemimajorAxiskm)
title('Semi-major axis STK')

figure(14)
plot(a.ArgofPerigeedeg)
title('Argument of Perigee STK')

for i = 1:1:1440

DifferenceEcc(i) =  abs( a.Eccentricity(i+2) - ecc(i)) ;
DifferenceIncdeg(i) = abs( a.Inclinationdeg(i+2) - inc(i) ) ;
DifferenceRAANdeg(i) = abs( a.RAANdeg(i+2) - raan(i) ) ;
DifferenceTrueAnomaly(i) = abs( a.TrueAnomalydeg(i+2) - true(i) ) ;
DifferenceSMA(i) = abs( a.SemimajorAxiskm(i+2) - sma(i) );
DifferenceAofPer(i) = abs(a.ArgofPerigeedeg(i+2) - argumentofperigee(i) ); 

end

ans = 1; 

end
