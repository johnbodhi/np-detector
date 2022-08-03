clear all; close all; clc

figure(1);

X = (-10:0.004:10-0.004)';

% H0...

y0 = normpdf(X,-1,sqrt(12)); plot(y0); hold on;

% H1...

y1 = 0.4 * normpdf(X,4,sqrt(0.25)) + 0.6 * normpdf(X,8,sqrt(0.25)); plot(y1);

% H0 and H1 realizations...

x0 = load('lab3_dataH0.txt'); x1 = load('lab3_dataH1.txt'); 

figure(2)

y0 = histogram(x0,100); y1 = histogram(x1,100);


% Find PD and PFA for H0 and H1...

N = size(x0,1);

t = 3.25; 
 
% Find PFA by deciding for H0 in X0

PFA = size( find(x0 > t),1) / N;

% Find PD by deciding for H1 in X1

PD = size( find(x1 > t),1) / N;


% We can use qfunc() to find analytical answers for PD and PFA for any
% threshold.

t = 3.25; 

% Count...

PFA = size( find(x0 > t),1) / N; 

PD = size( find(x1 > t),1) / N;

% Use Q function...

% FInd PFA0 by deciding for H0 wirh Q function...
 
PFA = qfunc( ( t - ( -1 ) ) / sqrt( 12 ) );

% Find PD1 by deciding for H1 with Q function...

PD = 0.4 * qfunc( ( t - 4 ) / sqrt( 0.25 ) ) + 0.6 * qfunc( ( t - 8 ) / sqrt( 0.25 ) ) ; 


% Compute PD and PFA for a range of thresholds...

t = (2:0.25:8.5)';

%t = (0:10/N:10-10/N)';

for i = 1:size(t,1)
    
    PFA(i) = size( find(x0 > t( i ) ),1) / N;

    PD(i) = size( find(x1 > t( i ) ),1) / N; 

end
PFA = PFA'; PD = PD';

figure(3);

ROC1 = [PFA PD]; 

plot( ROC1(:,1), ROC1(:,2), 'k');

% Build NP detector...

for i = 1:N    
    X11(i) = 386 *x1(i) - 47 * x1(i)^2;
end
X11 = X11';

for i = 1:N    
    X01(i) = 386 *x0(i) - 47 * x0(i)^2;   
end
X01 = X01';

t = 24*log( (5/2) * sqrt(pi / 2) * (1 / sqrt(24*pi)) * 3.25) + 767;

PFA1 =  size(find(X01 > t),1) / N; PD1 =  size(find(X11 > t),1) / N;

for i = 1:N    
    X12(i) = 770 * x1(i) - 47 * x1(i)^2;
end
X12 = X12';

for i = 1:N    
    X02(i) = 770 * x0(i) - 47 * x0(i)^2;   
end
X02 = X02';

t = 24*log( (5/3) * sqrt(pi / 2) * (1 / sqrt(24*pi)) * 3.25) + 3071;

PFA2 =  size(find(X02 > t),1) / N; PD2  =  size(find(X12 > t),1) / N;

PD = PD1 + PD2; PFA = PFA1 + PFA2;


% Find new threshold values...

T = (0.005:10/N:10.005 -10/N)';

t = (2:0.25:8.5)';

for i = 1:N    
    T1(i) = 24*log( (5/2) * sqrt(pi / 2) * (1 / sqrt(24*pi)) * T(i) ) + 767;
    T2(i) = 24*log( (5/3) * sqrt(pi / 2) * (1 / sqrt(24*pi)) * T(i) ) + 3071;     
end
T1 = T1';T2 = T2';


% Find PFA and PD for a range of new thresholds...

for i = 1:size(T,1)

    % Find PD by deciding for H1 with x1 > t...

    PD1(i) = size( find( X11 > T1(i) ), 1) / N;
    PD2(i) = size( find( X12 > T2(i) ), 1) / N;

    % Find PFA by deciding for H0 with x0 > t...

    PFA1(i) = size( find( X01 > T1(i) ), 1) / N;
    PFA2(i) = size( find( X02 > T2(i) ), 1) / N;

end
PD  = PD1 + PD2; PFA = PFA1 + PFA2;
PD = PD'; PFA = PFA';

ROC2 = [PFA PD]; 

figure(4); plot(ROC2(:,1),ROC2(:,2),'k');

