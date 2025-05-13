% calculating average temperature, pressure derivatives 
% using abers and hacker 2016 https://doi.org/10.1002/2015GC006171 
% copy this script to the unzipped directory from the supplemental materials and run it


% Load mineral database
[minpropar, compar]=ah16_loaddb('AbersHackerMacroJan2016.txt');
if isempty(minpropar)
    disp('ERROR on AH16_LOADDB')
    return 
end

% Set notional rock composition
mins={'an',};
modes=[100 ];       % These are volume fraction modes

% temperature dependence
tvals = linspace(0.01, 1200, 100);
pp = 0.01;  % GPa
gvals_t = zeros(size(tvals));
kvals_t = zeros(size(tvals));
for it = 1:numel(tvals)
    tt = tvals(it);  %deg C

    %% Calcuate - single P,T  -- also get rhos
    [modu,modhsm,modhsp,modvrm,modvrp,rhos]=ah16_rockvel(tt,pp, minpropar, mins,modes);

    gvals_t(it) = modu.g; 
    kvals_t(it) = modu.k; 
end 

gfit = polyfit(tvals, gvals_t, 1);
kfit = polyfit(tvals, kvals_t, 1);
figure()
plot(tvals, gvals_t,'k')
hold on 
plot(tvals, gfit(1) * tvals + gfit(2), '--k')
plot(tvals, kvals_t,'r')
plot(tvals, kfit(1) * tvals + kfit(2), '--r')

dGdT = gfit(1);
dKdT = kfit(1);


% pressure dependence
pvals = linspace(0, 4, 100);
tt = 0.01;  
gvals_p = zeros(size(tvals));
kvals_p = zeros(size(tvals));
for ip = 1:numel(pvals)
    pp = pvals(ip);  %deg C

    %% Calcuate - single P,T  -- also get rhos
    [modu,modhsm,modhsp,modvrm,modvrp,rhos]=ah16_rockvel(tt,pp, minpropar, mins,modes);

    gvals_p(ip) = modu.g; 
    kvals_p(ip) = modu.k; 
end 

% linear regressions
gfit_p = polyfit(tvals, gvals_p, 1);
kfit_p = polyfit(tvals, kvals_p, 1);

figure()
plot(pvals, gvals_p,'k')
hold on 
plot(pvals, gfit_p(1) * tvals + gfit_p(2), '--k')
plot(pvals, kvals_p,'r')
plot(pvals, kfit_p(1) * tvals + kfit_p(2), '--r')

dGdP = gfit_p(1);
dKdP = kfit_p(1);

disp('pressure and temperature dependence of shear modulus:')
disp('G = Go + dGdT * (T - 273) + dGdP * P')
disp(['   Go = ', num2str(gfit(2))])
disp(['   dGdT = ', num2str(gfit(1))])
disp(['   dGdP = ', num2str(gfit_p(1))])

disp('pressure and temperature dependence of bulk modulus:')
disp('K = Ko + dKdT * (T - 273) + dKdP * P')
disp(['   Ko = ', num2str(kfit(2))])
disp(['   dKdT = ', num2str(kfit(1))])
disp(['   dKdP = ', num2str(kfit_p(1))])
