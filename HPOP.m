%--------------------------------------------------------------------------
%
%               High Precision Orbit Propagator
%
% References:
% Montenbruck O., Gill E.; Satellite Orbits: Models, Methods and 
% Applications; Springer Verlag, Heidelberg; Corrected 3rd Printing (2005).
%
% Montenbruck O., Pfleger T.; Astronomy on the Personal Computer; Springer 
% Verlag, Heidelberg; 4th edition (2000).
%
% Seeber G.; Satellite Geodesy; Walter de Gruyter, Berlin, New York; 2nd
% completely revised and extended edition (2003).
%
% Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill
% , New York; 3rd edition(2007).
%
% http://sol.spacenvironment.net/jb2008/
%
% http://ssd.jpl.nasa.gov/?ephemerides
%
% Last modified:   2018/02/11   M. Mahooti
%
%--------------------------------------------------------------------------
clc
clear all                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
format long g
close all
profile on
tic

global const Cnm Snm AuxParam eopdata swdata SOLdata DTCdata APdata PC init_oe master_oe current_orbit K

SAT_Const
constants
load DE436Coeff.mat
PC = DE436Coeff;

%% read Earth gravity field coefficients
Cnm = zeros(181,181);
Snm = zeros(181,181);
fid = fopen('GGM03S.txt','r');
for n=0:180
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end
fclose(fid);

%% read Earth orientation parameters
fid = fopen('eop19990101.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
fclose(fid);

%% read space weather data
fid = fopen('sw19990101.txt','r');
%  ---------------------------------------------------------------------------------------------------------------------------------
% |                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
% | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
%  ---------------------------------------------------------------------------------------------------------------------------------
swdata = fscanf(fid,'%4i %3d %3d %5i %3i %3i %3i %3i %3i %3i %3i %3i %3i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4i %4f %2i %4i %6f %2i %6f %6f %6f %6f %6f',[33 inf]);
fclose(fid);

%% read space weather data
fid = fopen('SOLFSMY.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
%  ------------------------------------------------------------------------
SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);

%% READ GEOMAGNETIC STORM DTC VALUE
fid = fopen('DTCFILE.txt','r');
%  ------------------------------------------------------------------------
% | DTC YYYY DDD   DTC1 to DTC24v
%  ------------------------------------------------------------------------
DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
fclose(fid);

%% read space weather data
fid = fopen('SOLRESAP.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD  F10 F10B Ap1 to Ap8
%  ------------------------------------------------------------------------
APdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[12 inf]);
fclose(fid);

%% model parameters
AuxParam = struct('Mjd_UTC',0,'area_solar',0,'area_drag',0,'mass',0,'Cr',0,'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0,'Thrust',0,'direction',0, 'ThrustMag',1);

%% epoch state (X2)
fid = fopen('InitialState.txt','r');
tline = fgetl(fid);
year = str2num(tline(1:4));
mon = str2num(tline(6:7));
day = str2num(tline(9:10));
hour = str2num(tline(12:13));
min = str2num(tline(15:16));
sec = str2num(tline(18:23));
tline = fgetl(fid);
Y0(1) = str2num(tline);
tline = fgetl(fid);
Y0(2) = str2num(tline);
tline = fgetl(fid);
Y0(3) = str2num(tline);
tline = fgetl(fid);
Y0(4) = str2num(tline);
tline = fgetl(fid);
Y0(5) = str2num(tline);
tline = fgetl(fid);
Y0(6) = str2num(tline);
tline = fgetl(fid);
AuxParam.area_solar = str2num(tline(49:end));
tline = fgetl(fid);
AuxParam.area_drag = str2num(tline(38:end));
tline = fgetl(fid);
AuxParam.mass = str2num(tline(19:end));
tline = fgetl(fid);
AuxParam.Cr = str2num(tline(5:end));
tline = fgetl(fid);
AuxParam.Cd = str2num(tline(5:end));
fclose(fid);

%% epoch
Mjd_UTC = Mjday(year, mon, day, hour, min, sec);
Y0 = ECEF2ECI(Mjd_UTC, Y0);

AuxParam.Mjd_UTC = Mjd_UTC;
AuxParam.n       = 40;
AuxParam.m       = 40;
AuxParam.sun     = 1;
AuxParam.moon    = 1;
AuxParam.planets = 1;
AuxParam.sRad    = 1;
AuxParam.drag    = 1;
AuxParam.SolidEarthTides = 1;
AuxParam.OceanTides = 1;
AuxParam.Relativity = 1;
AuxParam.Thrust = 1;
AuxParam.direction = 1;
AuxParam.ThrustMag = 1;


Mjd0   = Mjd_UTC;
Step   = 60;   % [s]
N_Step = 1440; % 24 hours


%% shorten PC, eopdata, swdata, Cnm, and Snm
num = fix(N_Step*Step/86400)+2;
JD = Mjd_UTC+2400000.5;
i = find(PC(:,1)<=JD & JD<=PC(:,2),1,'first');
PC = PC(i:i+num,:);
mjd = (floor(Mjd_UTC));
i = find(mjd==eopdata(4,:),1,'first');
eopdata = eopdata(:,i:i+num);
i = find((year==swdata(1,:)) & (mon==swdata(2,:)) & (day==swdata(3,:)),1,'first');
swdata = swdata(:,i-3:i+num);
Cnm = Cnm(1:AuxParam.n+1,1:AuxParam.n+1);
Snm = Snm(1:AuxParam.n+1,1:AuxParam.n+1);

y0 = Y0;
sum = -60;
cta = 0;
k =0;

[sm_y0, ecc_y0, RA_y0, Inclination_y0, ArgPer_y0, TrueA_y0, T_y0, b_y0, p_y0, r_y0, h_y0 ] = Orbital_Elem(y0);

% set initial orbit 
init_oe = [sm_y0, ecc_y0, RA_y0, Inclination_y0, ArgPer_y0, TrueA_y0, T_y0,  b_y0, p_y0, r_y0, h_y0];

% set master orbit 
master_oe = [ sm_y0, 0 , RA_y0 ,  Inclination_y0 , ArgPer_y0, TrueA_y0, T_y0,  b_y0, p_y0, r_y0, h_y0];


% % % % % % Master = MasterOrbit ; 
% % % % % % Master.Value = [sm, ecc, RA, Inclination, ArgPer, TrueA, T];

csm = sm_y0;
K = [1, 1, 1, 1, 1];


% desired ascension 
ascension = 100;  


% calculate mean semi major-axis 
mean_sm = (7095+sm_y0)/2;

% propagation
% separation of the main propagation in smaller steps
for i = 0:Step:N_Step*Step
    stp = 1;
    n_stp = 60;
    
    %% OLD
% %     pta = cta;
% %     psm = csm;
    
    %% measure orbital elements
    [sm, ecc, RA, Inclination, ArgPer, TrueA, T,  b, p, r, h] = Orbital_Elem(y0);
    
    
    %% estimate the current orbit // input for input_matrix()
    current_orbit = [sm, ecc, RA, Inclination, ArgPer, TrueA, T,  b, p, r, h];
    
    
%% OLD
% %     cta = TrueA;
% %     csm = sm;
% %     cta
% %     sm
    
    % apply control 
% %     flag = Altitude_Control(sm,cta,pta,k);
% %     if flag == 1
% %         k = k + 1;
% %     end
% %     k;
% %     
    %% calculate Eph for 1 min
    Sub_Eph = Ephemeris(y0,n_stp,stp);
    vector = Sub_Eph(end,:);
    Eph(i/60+1,1) = vector(1) + sum ; 
    Eph(i/60+1,2:7) = vector(2:7);
    
    
    %% OLD
% %     if i >= 60
% %         
% %         delta_sm = csm - psm;
% %         % calculate thrust per step
% %         % check decrease of thrust power
% %         reduce_flag = calculate_thrust(AuxParam.mass, delta_sm , mean_sm , sum+60);
% %         reduce_flag
% %         AuxParam.ThrustMag
% %     
% %     end
       
%% 
    t = Eph(i/60+1,1);
    y0 = vector(2:7)';
    
    % timing sum
    sum = sum + vector(1);
    
    %% calculate gains 
if i > 60 
    K = gains (current_orbit, sum)
end
end

fid = fopen('SatelliteStates.txt','w');
for i=1:N_Step+1
    [year,mon,day,hr,min,sec] = invjday(Mjd0+Eph(i,1)/86400+2400000.5);
    fprintf(fid,'  %4d/%2.2d/%2.2d  %2d:%2d:%6.3f',year,mon,day,hr,min,sec);
    fprintf(fid,'  %14.3f%14.3f%14.3f%12.3f%12.3f%12.3f\n', Eph(i,2),Eph(i,3),Eph(i,4),Eph(i,5),Eph(i,6),Eph(i,7));       
end
fclose(fid);

[n, m] = size(Eph);
Eph_ecef = zeros(n,m);
for i=1:n
    Eph_ecef(i,1) = Eph(i,1);
    Eph_ecef(i,2:7) = ECI2ECEF(Mjd0+Eph_ecef(i,1)/86400, Eph(i,2:7)); 
end

% % % % % % True_EnvisatStates
% % % % % % dd = True_Eph-Eph_ecef(:,2:7);

% Plot orbit in ECI reference
figure(1)
plot3(Eph(:,2),Eph(:,3),Eph(:,4),'o-r')
grid;
title('Orbit ECI (inertial) (m)')

% Plot orbit in ECEF reference
figure(2) 
plot3(Eph_ecef(:,2),Eph_ecef(:,3),Eph_ecef(:,4),'-')
title('Orbit ECEF (m)')
xlabel('X');ylabel('Y');zlabel('Z');
grid

% % Plot Discrepancy of Precise and Propagated orbits
% figure(3)
% subplot(3,1,1);
% plot(dd(:,1));
% title('Discrepancy of Precise and Propagated ICEYE X2 Positions for 26.47 hours');
% axis tight
% xlabel('Time')
% ylabel('dX[m]')
% hold on
% subplot(3,1,2);
% plot(dd(:,2));
% axis tight
% xlabel('Time')
% ylabel('dY[m]')
% subplot(3,1,3);
% plot(dd(:,3));
% axis tight
% xlabel('Time')
% ylabel('dZ[m]')
% toc
% 
% tic
% lamda = zeros(n,1);
% phi = zeros(n,1);
% height = zeros(n,1);
% for i=1:n
%     [lamda(i),phi(i),height(i)] = Geodetic(Eph_ecef(i,2:4));
% end
% 
% figure(4)
% geoshow('landareas.shp','FaceColor',[0.5 1 0.5]);
% title('Satellite''s Ground Track')
% hold on
% plot(lamda*(180/pi),phi*(180/pi),'.r')
% % animation
% an = animatedline('Marker','*');
% for k = 1:n
%     addpoints(an,lamda(k)*(180/pi),phi(k)*(180/pi));
%     drawnow
%                                                                                                                                                                                      pause(0.01);
%     clearpoints(an);
% end
% toc
% 
% profile viewer
% profile off
