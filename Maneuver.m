%Maneuver : Computes the accelaration because of thrust, in order to commit phasing maneuver 
%             
%
% Inputs:
%  
% Output:
%   a               Acceleration (a=d^2r/dt^2)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------

function a = Maneuver(mu, tau, Y(1:3), Y(4:6))


tolerance = 1.0d-12;

u = 0;

imax = 20;

umax = +realmax;

umin = -realmax;

orbits = 0;

tdesired = tau;

threshold = tolerance * abs (tdesired);

r0 = norm(Y(1:3));

n0 = dot(Y(1:3), Y(4:6));

beta = 2 * (mu / r0) - dot(Y(4:6),Y(4:6));

if (beta ~= 0)
    umax = +1 / sqrt(abs(beta));
    
    umin = -1 / sqrt(abs(beta));
end

if (beta > 0)
    orbits = beta * tau - 2 * n0;
    
    orbits = 1 + (orbits * sqrt(beta)) / (pi * mu);
    
    orbits = floor (orbits / 2);
end

for i = 1:1:imax

    q = beta * u * u;
    
    q = q / (1 + q);

    n =  0;
    
    r = 1;
    
    l = 1;
    
    s = 1;
    
    d = 3;
    
    gcf = 1;
    
    k = -5;
    
    gold = 0;

    while (gcf ~= gold)

        k = -k;
        
        l = l + 2;
        
        d = d + 4 * l;
        
        n = n + (1 + k) * l;
        
        r = d / (d - n * r * q);
        
        s = (r - 1) * s;

        gold = gcf;
        
        gcf  = gold + s;

    end

    h0 = 1 - 2 * q;
    
    h1 = 2 * u * (1 - q);

    u0 = 2 * h0 * h0 - 1;
    
    u1 = 2 * h0 * h1;
    
    u2 = 2 * h1 * h1;
    
    u3 = 2 * h1 * u2 * gcf / 3;

    if (orbits ~= 0)
        u3 = u3 + 2 * pi * orbits / (beta * sqrt(beta));
    end

    r1 = r0 * u0 + n0 * u1 + mu * u2;
    
    dt = r0 * u1 + n0 * u2 + mu * u3;
    
    slope = 4 * r1 / (1 + beta * u * u);
    
    terror = tdesired - dt;

    if (abs (terror) < threshold)
        break;
    end;
    
    if ((i > 1) && (u  ==  uold) )
        break;
    end
    
    if ((i > 1) && (dt == dtold) )
        break;
    end

    uold  = u;
    
    dtold = dt;
    
    ustep = terror / slope;

    if (ustep > 0)
        umin = u;
        
        u = u + ustep;
        
        if (u > umax)
            u = (umin + umax) / 2;
        end
    else
        umax = u;
        
        u = u + ustep;
        
        if (u < umin)
            u = (umin + umax) / 2;
        end
    end

    if (i == imax)
        fprintf('\n\nmax iterations in twobody2 function');
    end

end

usaved = u;

f = 1.0 - (mu / r0) * u2;

gg = 1.0 - (mu / r1) * u2;

g  =  r0 * u1 + n0 * u2;

ff = -mu * u1 / (r0 * r1);

% final position and velocity vectors

for i = 1:1:3
    rf(i) = f  * ri(i) + g  * vi(i);

    vf(i) = ff * ri(i) + gg * vi(i);
end



%Delta - V calculation






    


% Acceleration 
a = 
end