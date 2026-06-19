%% 
%Initial declarations, don't worry about these. Only used internally.
clear;
syms x gmax
assume(gmax, "positive")
assumeAlso(gmax, "real")
curvesInward = false;
firstNegativeConcavityDistance = 0;
negativeConcavityVibLevel = -0.5;
useSmoothedCurveForMorseFit = false;

%Select input file
inDir = uigetfile("Choose input spreadsheet");
options = readcell(inDir);

%Get location to dump output
outDir = uigetdir(pwd,"Choose output directory");

%Kaiser correction - true if the bottom of the well seems suspect, otherwise
%leave false.
UseKaiserCorrection = logical(options{ind2sub(size(options),find(strcmp(options,'kaiser'))),2});

%Graphing parameters
plotVibLevelEvery = options{ind2sub(size(options),find(strcmp(options,'ladderspace'))),2}; %E.g. if set to 5, the graph will show a horizontal line for v=0, v=5, v=10, and so forth
vmax = options{ind2sub(size(options),find(strcmp(options,'vmax'))),2}; %Maximum vibrational level to try and calculate for. If the RKR results begin to curve inward the
%program will immediately halt, but still extrapolate up to vmax with the
%smoothing procedure.
tol = options{ind2sub(size(options),find(strcmp(options,'errortol'))),2}; %Error tolerance for integrals in case of singularity
space = options{ind2sub(size(options),find(strcmp(options,'space'))),2}; %Spacing between vibrational level increments (treating v as a continuous variable). Smaller = better.
UseSingularityMethod = options{ind2sub(size(options),find(strcmp(options,'intmethod'))),2}; %Watson (use only if you have two or more rovibrational corrections), Fleming (only first derivatives), or anything else for default

%Molecular isotope masses in amu
IsotopeAMass = options{ind2sub(size(options),find(strcmp(options,'mass1'))),2};
IsotopeBMass = options{ind2sub(size(options),find(strcmp(options,'mass2'))),2};

%Molecule name + state for output file naming
name = convertCharsToStrings(options{ind2sub(size(options),find(strcmp(options,'name'))),2});
state = convertCharsToStrings(options{ind2sub(size(options),find(strcmp(options,'state'))),2});

diary(convertCharsToStrings(outDir) + "\OUT_SYSLOG_" + name + "_" + state + "_" + string(datetime(now,'Format','uuuu-MM-dd_HHmmss','ConvertFrom','datenum')) + ".txt")

%Net charge on molecule in units of electron charge 1.6x10^-19 C
NetCharge = options{ind2sub(size(options),find(strcmp(options,'netcharge'))),2};

%Spectroscopic constants. THE MORE THE BETTER.
we = options{ind2sub(size(options),find(strcmp(options,'we'))),2}; %Fundamental vibrational frequency in cm^-1
xwe = options{ind2sub(size(options),find(strcmp(options,'xwe'))),2}; %First-order anharmonicity correction in cm^-1
ywe = options{ind2sub(size(options),find(strcmp(options,'ywe'))),2}; %Second-order anharmonicity correction
zwe = options{ind2sub(size(options),find(strcmp(options,'zwe'))),2}; %Third-order anharmonicity correction
awe = options{ind2sub(size(options),find(strcmp(options,'1we'))),2}; 
bwe = options{ind2sub(size(options),find(strcmp(options,'2we'))),2};
cwe = options{ind2sub(size(options),find(strcmp(options,'3we'))),2};
dwe = options{ind2sub(size(options),find(strcmp(options,'4we'))),2};
ewe = options{ind2sub(size(options),find(strcmp(options,'5we'))),2};
fwe = options{ind2sub(size(options),find(strcmp(options,'6we'))),2}; %Last anharmonicity correction

Be = options{ind2sub(size(options),find(strcmp(options,'Be'))),2}; %Rotational constant in cm^-1
ae = options{ind2sub(size(options),find(strcmp(options,'ae'))),2}; %First-order rotational correction in cm^-1
ye = options{ind2sub(size(options),find(strcmp(options,'ye'))),2}; %Second-order rotational correction
arot = options{ind2sub(size(options),find(strcmp(options,'1e'))),2}; %Third-order
brot = options{ind2sub(size(options),find(strcmp(options,'2e'))),2}; %Fourth-order

Te = options{ind2sub(size(options),find(strcmp(options,'Te'))),2}; %Electronic term energy

%EXPLICIT MORSE CONSTANTS (optional, for comparison graph only)
Die = options{ind2sub(size(options),find(strcmp(options,'De'))),2}; %Dissociation energy in cm^-1
re = options{ind2sub(size(options),find(strcmp(options,'re'))),2}; %Equilibrium internuclear bond length in Angstrom
ke = options{ind2sub(size(options),find(strcmp(options,'ke'))),2}; %Bond stiffness constant in cm^-1 / Angstrom^2

UseQuasiMorse = true;
if Die == 0 || re == 0 || ke == 0
    UseQuasiMorse = false;
end

%Constants
Cu = 16.857629206; % hbar^2/2 in amu Ang^2 cm^-1
me = 0.000548579909; %electron mass in amu

%System calculations
mu = (IsotopeAMass * IsotopeBMass)/(IsotopeAMass + IsotopeBMass - NetCharge*me); %Watson's charge-modified reduced mass in amu
a = sqrt(ke / (2*Die));
Y10 = we; %Just renaming energy derivative coefficients. Ref: Dunham energy levels
Y01 = Be;
Y20 = -1*xwe;
Y11 = -1*ae;
Y00 = 0.25*(Y01+Y20) - ((Y11*Y10)/(12*Y01)) + (1/Y01)*((Y11*Y10)/(12*Y01))^2; %Ref: Leroy, RKR1

disp("Calculated Y00 " + Y00);
disp("Molecular state" + state);
%stop();
if UseKaiserCorrection == true
    vmin = -0.5 - Y00/Y10;
else
    vmin = -0.5;
end

%Vibration-dependent parameters
G = @(vib) we.*(vib+(1/2)) - xwe.*((vib+(1/2)).^2) + ywe.*((vib+(1/2)).^3) - zwe.*((vib+(1/2)).^4) + awe.*((vib+(1/2)).^5) - bwe.*((vib+(1/2)).^6) ...
    + cwe.*((vib+(1/2)).^7) - dwe.*((vib+(1/2)).^8) + ewe.*((vib+(1/2)).^9) - fwe.*((vib+(1/2)).^10) + Y00; %Vibrational energy
dG = matlabFunction(diff(G(x))); %First derivative of vibrational energy
d2G = matlabFunction(diff(dG(x))); %Second derivative of vibrational energy

B = @(vib) Be - ae.*(vib+0.5) + ye.*(vib+0.5).^2 - arot.*(vib+0.5).^3 + brot.*(vib+0.5).^4; %Rotational constant
dB = matlabFunction(diff(B(x))); %First derivative of rotational constant

f = @(vib) sqrt(Cu / mu).*integral(@(x) 1./(sqrt(G(vib)-G(x))),vmin,vib,"RelTol",0,"AbsTol",tol); %First Klein integral
g = @(vib) sqrt(mu / Cu).*integral(@(x) (B(x))./(sqrt(G(vib)-G(x))),vmin,vib,"RelTol",0,"AbsTol",tol);%Second Klein integral

%Useful methods to remove singularity
fWatson = @(vib) sqrt(Cu / mu).*( (2/dG(vmin)) .* sqrt((G(vib)+Y00)) ...
    - 2.*integral(@(vp) sqrt(G(vib) - G(vp)) .* (d2G(vp) ./ (dG(vp).^2)), vmin, vib) );
gWatson = @(vib) sqrt(mu / Cu).*( (2*B(vmin)/dG(vmin)) .* sqrt((G(vib)+Y00)) ...
    - 2.*integral(@(vp) sqrt(G(vib) - G(vp)) .* (B(vp).*d2G(vp) - dB(vp).*dG(vp)) .* (1 ./ (dG(vp).^2)), vmin, vib) );
fFleming = @(vib) sqrt(Cu / mu) .* (1./dG(vib)) .* ( 2.*sqrt(G(vib)+Y00) ...
    + integral(@(vp) (dG(vib)-dG(vp))./(sqrt(G(vib)-G(vp))), vmin, vib));
gFleming = @(vib) sqrt(mu / Cu) .* (1./dG(vib)) .* ( 2.*B(vib).*sqrt(G(vib)+Y00) ...
    + integral(@(vp) (B(vp).*dG(vib)-B(vib).*dG(vp))./(sqrt(G(vib)-G(vp))), vmin, vib));

if UseSingularityMethod == "Watson"
    r1 = @(vib) -1.*fWatson(vib) + sqrt((fWatson(vib).^2) + (fWatson(vib)./gWatson(vib))); %First turning point for a vibrational level vib in Angstrom
    r2 = @(vib) fWatson(vib) + sqrt((fWatson(vib).^2) + (fWatson(vib)./gWatson(vib))); %Second turning point for a vibrational level vib in Angstrom
elseif UseSingularityMethod == "Fleming"
    r1 = @(vib) -1.*fFleming(vib) + sqrt((fFleming(vib).^2) + (fFleming(vib)./gFleming(vib))); %First turning point for a vibrational level vib in Angstrom
    r2 = @(vib) fFleming(vib) + sqrt((fFleming(vib).^2) + (fFleming(vib)./gFleming(vib))); %Second turning point for a vibrational level vib in Angstrom
else
    r1 = @(vib) -1.*f(vib) + sqrt((f(vib).^2) + (f(vib)./g(vib))); %First turning point for a vibrational level vib in Angstrom
    r2 = @(vib) f(vib) + sqrt((f(vib).^2) + (f(vib)./g(vib))); %Second turning point for a vibrational level vib in Angstrom

end
VMorse = @(r) Te+Die.*(1-exp(-a.*(r-re))).^2; %Morse potential energy function

if UseQuasiMorse
    Vsym = VMorse(rsym);
    %Potential asymptote
    asy = limit(Vsym, inf);
end

%Check for maximum of vibrational energy, if vmax exceeds this value then
%truncate vmax accordingly
extrema = double(solve(dG(gmax) == 0));
extrema = extrema(extrema > 0);
Gmax = extrema(1);
if vmax > Gmax
vmax = Gmax;
end

turningPoints = [];
energies = [];
purevibenergies=[];

computed_Bv = [];


computed_dGv = [];
computed_dBv = [];

computed_gFleming = [];
vib_n = [];
r1_s = [];
r2_s = [];
fF_s = [];
Gv_s = [];

tic;
disp("Computing potential energies and turning points...")
outputEvery = 1/space;
loopNo = 0;
lastInnerTurningPoint = inf;
lastVibrationalEnergy = -inf;
InwardCurvatureDetected = false;
disp("The minimum and maximum "+ vmin +" "+ vmax);
%stop();
for vIndex=vmin:space:vmax

    loopNo = loopNo + 1;
    if mod(loopNo, outputEvery) == 0
        disp("Entering loop " + loopNo);% + "vib-->" + vIndex);
    end

    vib = vIndex;

    %disp("vib->"+ vib);
    %disp("f->" +fFleming(vib))


    turningPoint1 = r1(vib);
    turningPoint2 = r2(vib);

    computed_Gv =  G(vib);

    computed_fFleming =  fFleming(vib);


    if ~isnan(turningPoint1) && ~isnan(turningPoint2)

        turningPoints(end + 1) = turningPoint1;
        turningPoints(end + 1) = turningPoint2;
        vib_n(end+1) =  vIndex;
        r1_s(end+1)=turningPoint1;
        r2_s(end+1)=turningPoint2;
        Gv_s(end+1)=computed_Gv;
        fF_s(end+1)=computed_fFleming;

        if G(vib) < lastVibrationalEnergy
            error("Vibrational energy maximum for approximation reached at v=" + vib + ". Try again with a lower vmax.")
        end

        if turningPoint1 > lastInnerTurningPoint
            firstNegativeConcavityDistance = r1(round(0.65*vib));
            negativeConcavityVibLevel = round(0.65*vib);
            disp("Inward curvature detected at vibrational level " + round(vib) +", stopping loops and extrapolating manually from v=" + round(0.65*vib))
            InwardCurvatureDetected = true;
            turningPoints = turningPoints(1:end-2);
            disp("I stop here!");
            break
            
        else
            lastInnerTurningPoint = turningPoint1;
        end

        if UseQuasiMorse
            E1 = VMorse(turningPoint1);
            E2 = VMorse(turningPoint2);
            energies(end+1) = (Te+E1); %- (Cu * J*(J+1) /(mu*turningPoint1^2));
            energies(end+1) = (Te+E2); %- (Cu * J*(J+1) /(mu*turningPoint2^2));
        end
    end

    purevibenergies(end+1) = Te+G(vib);
    purevibenergies(end+1) = Te+G(vib);
    lastVibrationalEnergy = G(vib);
  
    
    %computed_gFleming(end + 1) = gFleming(vib);
    
    %vib_n(end+1) =  vib;


end
data_2_file = [vib_n(:), Gv_s(:), r1_s(:), r2_s(:), fF_s(:)];
filename = sprintf('C:\\msys64\\home\\dadel\\test_data_matlab_%s_%s.dat', name, state);
writematrix(data_2_file, filename, 'Delimiter', 'tab');


%Debug printing.....

%disp("turning points, energies, purevibenergies")
%dlmwrite('..\c++_exec_rkr_method\vib_n.txt', vib_n, 'delimiter', '\n');  % tab-separated

%dlmwrite('..\c++_exec_rkr_method\computed_G.txt', computed_Gv, 'delimiter', '\n');  % tab-separated
%dlmwrite('..\c++_exec_rkr_method\computed_B.txt', computed_Bv, 'delimiter', '\n');  % tab-separated
%dlmwrite('..\c++_exec_rkr_method\computed_dG.txt', computed_dGv, 'delimiter', '\n');  % tab-separated
%dlmwrite('..\c++_exec_rkr_method\computed_dB.txt', computed_dBv, 'delimiter', '\n');  % tab-separated
%dlmwrite('..\c++_exec_rkr_method\computed_fFleming_mlab.txt', real(computed_fFleming), 'delimiter', '\n');  % tab-separated
%dlmwrite('..\c++_exec_rkr_method\computed_gFleming_mlab.txt', real(computed_gFleming), 'delimiter', '\n');  % tab-separated
%dlmwrite('turningPoints.txt', turningPoints, 'delimiter', '\n');  % tab-separated
%dlmwrite('energies.txt', energies, 'delimiter', '\n');  % tab-separated
dlmwrite('purevibenergies.txt', purevibenergies, 'delimiter', '\n');  % tab-separated

disp("done saving calculations..")


disp("Potential curve computed. Time elapsed: " + toc + " seconds.")
disp("Avg loop time: " + (toc/loopNo) + " seconds.")

tic;
disp("Interpolating curve...")
%Organize by ascending internuclear distance
[turningPointsOrdered,I] = sort(turningPoints);

purevibenergiesOrdered = purevibenergies(I);

data_2_file = [turningPointsOrdered(:), purevibenergiesOrdered(:)];
filename = sprintf('C:\\msys64\\home\\dadel\\EvsrOrdered_%s_%s.dat', name, state);
writematrix(data_2_file, filename, 'Delimiter', 'tab');


%stop();

%dlmwrite('..\c++_exec_rkr_method\turningPointsOrdered.txt', turningPointsOrdered, 'delimiter', '\n');  % tab-separated
%dlmwrite('..\c++_exec_rkr_method\purevibenergiesOrdered.txt', purevibenergiesOrdered, 'delimiter', '\n');  % tab-separated

disp("Ordered calculations file write finished")


if UseQuasiMorse
    energiesOrdered = energies(I);
end
%Find equilibrium internuclear distance

%Differentiation matrix for cubic spline
D = diag(3:-1:1,1);

%Interpolate RKR potential
isNan = isnan(turningPointsOrdered);
%VInterp = csape(real(turningPointsOrdered(~isNan)), purevibenergiesOrdered(~isNan), 'variational');
VInterp = spline(real(turningPointsOrdered(~isNan)), purevibenergiesOrdered(~isNan));
dVInterp = VInterp;
d2VInterp = VInterp;

dVInterp.coefs = VInterp.coefs*D;
d2VInterp.coefs = VInterp.coefs*D*D;

distanceDomain = min(turningPointsOrdered):0.0001:max(turningPointsOrdered);
V = ppval(VInterp, distanceDomain); %Potential energy spline array
dV = ppval(dVInterp, distanceDomain); %First derivative of potential energy spline array with respect to distance
d2V = ppval(d2VInterp, distanceDomain); %Second derivative of potential energy spline array with respect to distance


% You already have:
% VInterp = csape(real(turningPointsOrdered(~isNan)), purevibenergiesOrdered(~isNan), 'variational');

% Read C++ spline data
cpp_file = 'C:\msys64\home\dadel\V_dV_d2V_spline_cpp.dat';  % adjust path if needed
cpp_data = readmatrix(cpp_file);

r_cpp   = cpp_data(:,1);
V_cpp   = cpp_data(:,2);
dV_cpp  = cpp_data(:,3);
d2V_cpp = cpp_data(:,4);

% Evaluate MATLAB splines on the same r_cpp grid
V_mat   = ppval(VInterp,  r_cpp);
dV_mat  = ppval(dVInterp, r_cpp);
d2V_mat = ppval(d2VInterp, r_cpp);

% ---- Plot V ----
figure;
plot(r_cpp, V_mat, 'k-', 'LineWidth', 1.5); hold on;
plot(r_cpp, V_cpp, 'r--', 'LineWidth', 1.2);
xlabel('r'); ylabel('V(r)');
legend('MATLAB VInterp', 'C++ V\_spline', 'Location', 'Best');
title('V(r): MATLAB vs C++');
grid on;

% ---- Plot dV ----
figure;
plot(r_cpp, dV_mat, 'k-', 'LineWidth', 1.5); hold on;
plot(r_cpp, dV_cpp, 'r--', 'LineWidth', 1.2);
xlabel('r'); ylabel('dV/dr');
legend('MATLAB dVInterp', 'C++ dV\_spline', 'Location', 'Best');
title('dV/dr: MATLAB vs C++');
grid on;

% ---- Plot d2V ----
figure;
plot(r_cpp, d2V_mat, 'k-', 'LineWidth', 1.5); hold on;
plot(r_cpp, d2V_cpp, 'r--', 'LineWidth', 1.2);
xlabel('r'); ylabel('d^2V/dr^2');
legend('MATLAB d2VInterp', 'C++ d2V\_spline', 'Location', 'Best');
title('d^2V/dr^2: MATLAB vs C++');
grid on;

hold on;




VInterpData = [V(:)];
% Choose a filename (mimicking your other outputs)
outVInterpFile = "C:\\msys64\\home\\dadel\\EInterp_data_matlab" + name + "_" + state + ".txt";
% Write as tab-separated text
writematrix(VInterpData, outVInterpFile, 'Delimiter', 'tab');
disp("Wrote VInterp interpolated data to: " + outVInterpFile);

%disp("First derivative" + dV + "\n Second Derivative \n"+ d2V);

% Build 2-column [r, E(r)] array from spline-evaluated data
VInterpData = [distanceDomain(:), V(:)];
% Choose a filename (mimicking your other outputs)
outVInterpFile = "C:\\msys64\\home\\dadel\\EvsrInterp_data_matlab" + name + "_" + state + ".txt";
% Write as tab-separated text
writematrix(VInterpData, outVInterpFile, 'Delimiter', 'tab');
disp("Wrote VInterp interpolated data to: " + outVInterpFile);


VInterpData = [dV(:)];
filename = sprintf('C:\\msys64\\home\\dadel\\first_deriv_%s_%s.txt', name, state);
writematrix(VInterpData, filename, 'Delimiter', 'tab');

VInterpData = [distanceDomain(:), d2V(:)];
filename = sprintf('C:\\msys64\\home\\dadel\\second_deriv_%s_%s.txt', name, state);
writematrix(VInterpData, filename, 'Delimiter', 'tab');


% for tenPow=0:6
% for expandVal=0:1/(10^tenPow):10000
%     wellBottomIndex = find(dV < expandVal & dV > -expandVal);
%     if length(wellBottomIndex) >= 2
%         break
%     end
% end
% if length(wellBottomIndex) <= 2
%     break
% end
% end

%FIND THE INDEX OF THE WELL BOTTOM WITH FIRST DERIVATIVE
for p=1:1:length(distanceDomain)-1
        if dV(p) < 0 && dV(p+1) > 0
            wellBottomIndex = [p p+1];
            break
        end
end
[wellMinimum, indexOfWellMinimum] = mink(abs(V),2);

%FIND FIRST OCCURRENCE OF NEGATIVE CONCAVITY
preBottomNegativeConcavityDistances = [];
if InwardCurvatureDetected == false
    for l=1:1:length(distanceDomain)-1
        if d2V(l) < 0 && d2V(l+1) > 0 && distanceDomain(l) < distanceDomain(wellBottomIndex(2))
            preBottomNegativeConcavityDistances(end+1) = distanceDomain(l);
        end
    end
end

if ~isempty(preBottomNegativeConcavityDistances)
firstNegativeConcavityDistance = preBottomNegativeConcavityDistances(end);
end

if firstNegativeConcavityDistance ~= 0
    [~, indexInTurningPointList] = min(abs(turningPointsOrdered - firstNegativeConcavityDistance));
    if negativeConcavityVibLevel == -0.5
        negativeConcavityVibLevel = floor(vmax - indexInTurningPointList*space);
        disp("Negative concavity located near v=" + negativeConcavityVibLevel + ", proceeding to smooth. This may take a long time.")
    end

    %Finds the best integer inverse-power fit using the ratio of each
    %derivative to the previous if concavity ever went below zero

    scale = 0;
    numDerivLoop = 0;
    for i=floor(negativeConcavityVibLevel*0.5):floor(negativeConcavityVibLevel*0.9)
        numDerivLoop = numDerivLoop + 1;
        repulsiveTurningPoint = r1(i);
        scale = scale + -1*(ppval(d2VInterp, repulsiveTurningPoint)/ppval(dVInterp, repulsiveTurningPoint))*repulsiveTurningPoint - 1;
    end
    avgScale = scale/numDerivLoop;
    integerPower = round(avgScale); %Integer power for hyperbolic fit to repulsive branch. Useful for curve smoothing.

    smoothRepulsiveFitEqn = "(A/(x^(" + integerPower + "))) + B";
    xForSmoothing = [r1(floor(negativeConcavityVibLevel*0.85)-1), r1(floor(negativeConcavityVibLevel*0.85))];
    yForSmoothing = [Te+G(floor(negativeConcavityVibLevel*0.85)-1), Te+G(floor(negativeConcavityVibLevel*0.85))];

    smoothingFitOptions =  fitoptions('Method','NonlinearLeastSquares', 'StartPoint', [1e5, -1e4]);
    %smoothFit = fit(transpose(xForSmoothing), transpose(yForSmoothing), smoothRepulsiveFitEqn, smoothingFitOptions);

    E1 = yForSmoothing(1);
    E2 = yForSmoothing(2);
    d1 = xForSmoothing(1);
    d2 = xForSmoothing(2);

    Asmooth = (E1 - E2)/((1/d1^integerPower) - (1/d2^integerPower)) ;
    Bsmooth = E2 - Asmooth/(d2^integerPower);

    disp(Asmooth + "and " + Bsmooth + "\n");

    smoothFitFn = @(r) (Asmooth/r^(integerPower)) + Bsmooth;
    disp("The required power: "+ integerPower);
    %Conditions and symbolic function
    syms rsym
    assume(rsym, "positive")
    assumeAlso(rsym, "real")
    smoothFitSymFn = smoothFitFn(rsym);

    %Calculate extrapolated inner turning points
    smoothedTurningPoints = turningPoints;
    correspondingIndex = 1;
    lastVibrationalEnergy = -inf;
    for v=vmin:space:vmax
        if v >= floor(negativeConcavityVibLevel*0.85)
            innerTurningPoint = double(solve(smoothFitSymFn == Te+G(v)));
            if UseSingularityMethod == "Watson"
                distance = 2*fWatson(v);
            elseif UseSingularityMethod == "Fleming"
                distance=2*fFleming(v);
            else
                distance = 2*f(v);
            end
            smoothedTurningPoints(correspondingIndex) = innerTurningPoint;
            smoothedTurningPoints(correspondingIndex + 1) = innerTurningPoint + distance;
            purevibenergies(correspondingIndex) = Te+G(v);
            purevibenergies(correspondingIndex+1) = Te+G(v);
            if G(v) < lastVibrationalEnergy
                error("Maximum vibrational energy for current approximation reached at v=" + v + ". Try again with a lower vmax.")
            end
            lastVibrationalEnergy = G(v);
        end
        correspondingIndex = correspondingIndex + 2;
    end
    [smoothedTurningPointsOrdered, I2] = sort(smoothedTurningPoints);
    purevibenergiesSmoothedOrdered = purevibenergies(I2);

   useSmoothedCurveForMorseFit = true;
end

wellBottom = (distanceDomain(wellBottomIndex(1)) + distanceDomain(wellBottomIndex(2)))/2;

%Find dissociation energy using Morse potential fit
k = (d2V(wellBottomIndex(1))+d2V(wellBottomIndex(2)))/2; %This is in cm^-1 / Ang^2
kOut = k * 1.986e-23 * 1e10 * 1e8 * 1e5; %Convert to dyn/cm for output

fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',0,...
    'StartPoint',V(end)-Te);

atxt = "sqrt(" + k + "/(2*D))";
morse = "(D*(1-exp(-" + atxt + "*(x-" + wellBottom + ")))^2) +" + Te;
morseFit = fit(transpose(distanceDomain),transpose(V),morse,fo);

dissociationEnergy = coeffvalues(morseFit);

disp("Curve successfully interpolated. Time elapsed: " + toc + " seconds.")

disp("Plotting data...")
%Plotting data for fitted curve:
if InwardCurvatureDetected
xfit = linspace(min(smoothedTurningPointsOrdered),max(smoothedTurningPointsOrdered),1e6);
else
xfit = linspace(min(turningPointsOrdered),max(turningPointsOrdered),1e6);
end
yfit = feval(morseFit,xfit);

%Get vibrational levels properly ordered
hLines=[];
hLines(end+1)=Te+G(0); disp("Making zero-point energy at E = " + hLines(end))
dissociationReached = false;
for i = 1:vmax
    if G(i) >= dissociationEnergy
        dissociationReached = true;
        maxVibLevel = i-1;
        break
    else
        hLines(end+1) = Te+G(i); disp("Making vibrational energy level at E = " + hLines(end) + " for vibrational level " + i)
    end
end

%COMPARISON PLOT
figure

if UseQuasiMorse
    plot(turningPointsOrdered, energiesOrdered, 'DisplayName', 'Morse from Full Constants')
    hold on
end

plot(turningPointsOrdered, purevibenergiesOrdered, 'DisplayName', 'RKR Raw')
hold on
plot(xfit,yfit, 'DisplayName', 'Morse Fit to RKR')
if firstNegativeConcavityDistance ~= 0
    hold on
    plot(smoothedTurningPointsOrdered, purevibenergiesSmoothedOrdered, 'DisplayName', "Smoothed + Extrapolated RKR Past v="+floor(negativeConcavityVibLevel*0.85))
end
legend
title("Comparison Graph")

%PURE RKR PLOT (No ladder since this can vary unpredictably)
figure
plot(turningPointsOrdered, purevibenergiesOrdered, 'DisplayName', 'RKR Raw')
title("Raw RKR Data")
legend
%LADDER MORSE PLOT FROM EXPLICITLY DEFINED VALUES
if UseQuasiMorse
    figure
    plot(turningPointsOrdered, energiesOrdered, 'DisplayName', 'Morse from Full Constants')
    l = legend;
    l.AutoUpdate = 'off';
    hold on
    %Always plot zero-point energy
    vibLevel = hLines(1);
    intercepts = solve(Vsym == vibLevel, rsym);
    right = double(vpa(intercepts(1),5));
    left = double(vpa(intercepts(2),5));
    %display(left + "," + right)
    if isreal(left) && isreal(right)
        plot([left right],[vibLevel vibLevel], "red")
        text(left-0.075,vibLevel,"v =" + 0);
    end
    hold on
    for i=2:length(hLines)
        %Plot vibrational level every plotVibLevelEvery increments (stupid name, I
        %know)
        if mod(i-1, plotVibLevelEvery) == 0
            vibLevel = hLines(i);
            intercepts = solve(Vsym == vibLevel, rsym);
            right = double(vpa(intercepts(1),5));
            left = double(vpa(intercepts(2),5));
            %display(left + "," + right)
            if isreal(left) && isreal(right)
                plot([left right],[vibLevel vibLevel], "blue")
                text(left-0.075,vibLevel,"v =" + (i-1));
            end
        end
    end
    %Always plot top vibrational level
    vibLevel = hLines(end);
    intercepts = solve(Vsym == vibLevel, rsym);
    right = double(vpa(intercepts(1),5));
    left = double(vpa(intercepts(2),5));
    %display(left + "," + right)
    if isreal(left) && isreal(right)
        plot([left right],[vibLevel vibLevel], "red")
        text(left-0.075,vibLevel,"v =" + (length(hLines)-1));
    end
    title("Explicit Morse")
end

%LADDER MORSE PLOT FOR FITTED MORSE
figure
plot(xfit,yfit, 'DisplayName', 'Morse Fit to RKR')
l =legend;
l.AutoUpdate = 'off';
hold on
vibLevel = hLines(1);
%Always plot zero point energy
[~, intercepts] = mink(abs(yfit-vibLevel),50);
right = xfit(max(intercepts));
left = xfit(min(intercepts));
%display(left + "," + right)
if isreal(left) && isreal(right)
    plot([left right],[vibLevel vibLevel], "red")
    text(left-0.075,vibLevel,"v =" + 0);
end
hold on
for i=2:length(hLines)
    %Plot vibrational level every plotVibLevelEvery increments (stupid name, I
    %know)
    if mod(i-1, plotVibLevelEvery) == 0
        vibLevel = hLines(i);
        [~, intercepts] = mink(abs(yfit-vibLevel),50);
        right = xfit(max(intercepts));
        left = xfit(min(intercepts));
        %display(left + "," + right)
        if isreal(left) && isreal(right)
            plot([left right],[vibLevel vibLevel], "blue")
            text(left-0.075,vibLevel,"v =" + (i-1));
        end
    end
end
%Always plot top vibrational level
vibLevel = hLines(end);
[~, intercepts] = mink(abs(yfit-vibLevel),50);
right = xfit(max(intercepts));
left = xfit(min(intercepts));
%display(left + "," + right)
if isreal(left) && isreal(right) && left < wellBottom && right > wellBottom
    plot([left right],[vibLevel vibLevel], "red")
    text(left-0.075,vibLevel,"v =" + (length(hLines)-1));
end
title("Calculated Morse Potential")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(InwardCurvatureDetected== true)
    data = [smoothedTurningPointsOrdered(:), purevibenergiesSmoothedOrdered(:)];
    dlmwrite('C:\msys64\home\dadel\Evsr_data_matlab.dat', data, 'delimiter', '\t');

else 
    data = [turningPointsOrdered(:), purevibenergiesOrdered(:)];
    dlmwrite('C:\msys64\home\dadel\Evsr_data_matlab.dat', data, 'delimiter', '\t');
end
%LADDER GRAPH FROM SMOOTHED RKR CURVE
if firstNegativeConcavityDistance ~= 0
    figure
    plot(smoothedTurningPointsOrdered, purevibenergiesSmoothedOrdered, 'DisplayName', "Smoothed + Extrapolated RKR Past v="+floor(negativeConcavityVibLevel*0.85))
    leg = legend;
    leg.AutoUpdate = "off";
    title("Smoothed RKR Curve")
    hold on
    for i=1:length(hLines)
        if mod(i-1, plotVibLevelEvery) == 0
            vibLevel = Te+G(i-1);
            if (i-1 < negativeConcavityVibLevel)
                left = r1(i-1);
                right = r2(i-1);
            else
                left = double(solve(smoothFitSymFn == vibLevel));
                right = left + 2*f(i-1);
            end

            plot([left right],[vibLevel vibLevel], "--r")
            text(left-0.075,vibLevel,"v =" + (i-1));
        end
    end
end

%OUTPUTS
if UseQuasiMorse
    disp("Potential Asymptote from Explicit Morse: " + double(asy))
end

disp("Asymptote for computed Morse potential: " + (Te+dissociationEnergy) + " cm^-1")
disp("Maxmimum energy for attractive branch in graphed range for raw RKR: " + purevibenergiesOrdered(end) + " cm^-1")
disp("Zero-point energy: " + hLines(1) + " cm^-1")
disp("Force constant at well minimum: " + kOut + " dyn/cm")

if dissociationReached
    disp("Maxmimum vibrational level: " + maxVibLevel)
end

disp("Equilibrium internuclear distance: " + wellBottom + " Angstrom")

%Build and output data

%Raw RKR Data
outRawMat = [transpose(round(turningPointsOrdered,8)) transpose(purevibenergiesOrdered)];
writematrix(outRawMat, outDir + "\OUT_RAW_"+name+"_"+state+".txt", 'Delimiter', 'tab')
%Smoothed data as well as smoothed data relative to the minimum energy for
%the ground state if applicable. If the state in question is a ground state
%already, smoothed and smoothed-zeroed data is identical.
if firstNegativeConcavityDistance ~= 0
    outSmoothMat = [transpose(round(smoothedTurningPointsOrdered,8)) transpose(purevibenergiesSmoothedOrdered)];
    writematrix(outSmoothMat, outDir + "\OUT_SMOOTH_"+name+"_"+state+".txt", 'Delimiter', 'tab')

    outZeroedSmoothMat = [transpose(round(smoothedTurningPointsOrdered,8)) transpose(purevibenergiesSmoothedOrdered-Te)];
    writematrix(outZeroedSmoothMat, outDir + "\OUT_SMOOTH_ZEROED_"+name+"_"+state+".txt", 'Delimiter', 'tab')
end

outMorseCompMat = [transpose(xfit) yfit];
%writematrix(outMorseCompMat, outDir + "\OUT_COMP_MORSE_"+name+"_"+state+".txt", 'Delimiter', 'tab')

