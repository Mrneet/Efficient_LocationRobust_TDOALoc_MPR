% Generation of Figs. 10 in the paper
% Y. Sun, K. C. Ho, G. Wang. J. Chen, Y. Yang, L. Chen, and Q. Wan, 
% "Computationally attractive and location robust estimator for IoT device positioning," 
% IEEE Internet Things J., Nov. 2021.
%
% Yimao Sun and K. C. Ho   04-08-2022
%
%       Copyright (C) 2022
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA.
%       hod@missouri.edu
%

clear all;
% close all;
% clc

rng('default');
% warning off

clor = [0, 114, 189;217, 83, 25;237, 177, 32;126, 47, 142;119, 172, 48;77, 190, 238;162, 20, 47]/256;
rng('default');

senPos = [
        % minimum number of sensors
        10.23    38.38    16.29
        46.64   -87.12    62.94
        124.02  -7.98     81.16
        105.02  -51.72    26.47
       -81.56    104.48  -80.49
        ]';

thetaN = (-180:10:180)*pi/180;
phiN = [-89,-85:5:85,89]*pi/180;

[N,M] = size(senPos);

% setting
sigma_sqr = -10;  % 10log(m^2)
range = 500;         % m
% range = 10000;
mon = 5000;

R = length(phiN);
K = length(thetaN);

% Monte-Carlo Simulation

aveNse = 0;
for l=1:mon
    aveNse = aveNse + randn(M,1); 
end
aveNse = aveNse/mon/sqrt(2);
PP = aveNse(2:end) - aveNse(1);

disp('Simulation is running ...');
senPosTmp = senPos;

for t = 1:R
    phi = phiN(t);

    for k = 1:K
        disp(['phi: ',num2str(phiN(t)*180/pi),'(deg), ',num2str(t),'/',num2str(R),'; ','theta: ',num2str(thetaN(k)*180/pi),'(deg), ',num2str(k),'/',num2str(K),' ...']);
        theta = thetaN(k);
        phiTmp = phi;
        
        souLoc = range * [cos(theta)*cos(phi); sin(theta)*cos(phi); sin(phi)] + senPos(:,1);
        souLocTmp = souLoc;
        r = sqrt(sum((souLoc-senPos).^2,1))';
        rd = r(2:end) - r(1);
        rTmp = r;
            
        Q = 10^(sigma_sqr/10) * (ones(M-1, M-1)+eye(M-1))/2;
        
        % calculate CRLB
        CRB = ConsCRLB( senPos, souLoc, Q );
        CRLB_a(k,t) = CRB(1,1)+CRB(2,2);
        CRLB_g(k,t) = CRB(3,3);
        
        % Position and DOA estimate
        nsePwr = 10^(sigma_sqr/10);
        
        % SCO-MPR Method
        rng('default');
        for i = 1:mon
            % measured TDOAs
            tmp=randn(M,1);
            rdNse = sqrt(nsePwr) * ((tmp(2:M)-tmp(1))/sqrt(2)-PP);
            rd_m = rd + rdNse;
            
            % SCO
            [mprSol,~] = TDOA_SCO_MPR( senPos, rd_m, Q );
            Th1 = mprSol(1);Ph1 = mprSol(2);g1 = mprSol(3);
            if abs(theta - Th1) > pi
                eTh1(k,i) = (2*pi-abs(theta-Th1))^2;
            else
                eTh1(k,i) = (theta - Th1)^2;
            end
            ePh1(k,i) = (phi - Ph1)^2;
            eg1(k,i) = (1/r(1) - g1)^2;
            uTh1(k,i) = Th1;
            uPh1(k,i) = Ph1;
            ug1(k,i) = g1;
            
            %% SUM-MPR Method
            [Th2, Ph2, g2, ~] = TDOA_SUM_MPR( senPos, rd_m, Q );
            if abs(theta - Th2) > pi
                eTh2(k,i) = (2*pi-abs(theta-Th2))^2;
            else
                eTh2(k,i) = (theta - Th2)^2;
            end
            ePh2(k,i) = abs(phi - Ph2)^2;
            eg2(k,i) = (1/r(1) - g2)^2;
            uTh2(k,i) = Th2;
            uPh2(k,i) = Ph2;
            ug2(k,i) = g2;

            %% GTRS-MPR Method
            [Th3, Ph3, g3, ~] = TDOA_GTRS_MPR( senPos, rd_m, Q );
            if abs(theta - Th3) > pi
                eTh3(k,i) = (2*pi-abs(theta-Th3))^2;
            else
                eTh3(k,i) = (theta - Th3)^2;
            end
            ePh3(k,i) = (phi - Ph3)^2;
            eg3(k,i) = (1/r(1) - g3)^2;
            uTh3(k,i) = Th3;
            uPh3(k,i) = Ph3;
            ug3(k,i) = g3;
        end
    end
    
    % calculate MSE
    % MSE of angle
    mse_a1(:,t) = mean(eTh1+ePh1,2);
    mse_a2(:,t) = mean(eTh2+ePh2,2);
    mse_a3(:,t) = mean(eTh3+ePh3,2);
    
    % MSE of g
    mse_g1(:,t) = mean(eg1,2);
    mse_g2(:,t) = mean(eg2,2);
    mse_g3(:,t) = mean(eg3,2);
end

figure;
subplot(2,1,1);
plot(thetaN*180/pi, 10*log10(mean(mse_a1./CRLB_a,2)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SCO-MPR'); hold on; grid on;
plot(thetaN*180/pi, 10*log10(mean(mse_a2./CRLB_a,2)), 'v', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');
plot(thetaN*180/pi, 10*log10(mean(mse_a3./CRLB_a,2)), 's', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
xlabel('\theta^o (deg)', 'FontSize', 13);ylabel('$\bar{R}_a(\theta^o)$','Interpreter','Latex', 'FontSize', 13);
xlim([min(thetaN),max(thetaN)]*180/pi);set(gca,'XTick',-180:60:180);
ylim([-0.5,3]);
lgd = legend('show');
set(lgd,'FontSize',11, 'Location', 'Northeast');
subplot(2,1,2);
plot(thetaN*180/pi, 10*log10(mean(mse_g1./CRLB_g,2)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SCO-MPR'); hold on; grid on;
plot(thetaN*180/pi, 10*log10(mean(mse_g2./CRLB_g,2)), 'v', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');
plot(thetaN*180/pi, 10*log10(mean(mse_g3./CRLB_g,2)), 's', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
xlabel('\theta^o (deg)', 'FontSize', 13);ylabel('$\bar{R}_g(\theta^o)$','Interpreter','Latex', 'FontSize', 13);
xlim([min(thetaN),max(thetaN)]*180/pi);set(gca,'XTick',-180:60:180);
ylim([-0.5,3]);
lgd = legend('show');
set(lgd,'FontSize',11, 'Location', 'Northeast');

figure;
subplot(2,1,1);
plot(phiN*180/pi, 10*log10(mean(mse_a1./CRLB_a,1)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SCO-MPR'); hold on; grid on;
plot(phiN*180/pi, 10*log10(mean(mse_a2./CRLB_a,1)), 'v', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');
plot(phiN*180/pi, 10*log10(mean(mse_a3./CRLB_a,1)), 's', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
xlabel('\phi^o (deg)', 'FontSize', 13);ylabel('$\bar{R}_a(\phi^o)$','Interpreter','Latex', 'FontSize', 13);
xlim([-90,90]);set(gca,'XTick',-90:30:90);
ylim([-0.5,3]);
lgd = legend('show');
set(lgd,'FontSize',11, 'Location', 'Northeast');
subplot(2,1,2);
plot(phiN*180/pi, 10*log10(mean(mse_g1./CRLB_g,1)), 'o', 'LineWidth', 1.5, 'DisplayName', 'SCO-MPR'); hold on; grid on;
plot(phiN*180/pi, 10*log10(mean(mse_g2./CRLB_g,1)), 'v', 'LineWidth', 1.5, 'DisplayName', 'SUM-MPR');
plot(phiN*180/pi, 10*log10(mean(mse_g3./CRLB_g,1)), 's', 'LineWidth', 1.5, 'DisplayName', 'GTRS-MPR');
xlabel('\phi^o (deg)', 'FontSize', 13);ylabel('$\bar{R}_g(\phi^o)$','Interpreter','Latex', 'FontSize', 13);
xlim([-90,90]);set(gca,'XTick',-90:30:90);
ylim([-0.5,3]);
lgd = legend('show');
set(lgd,'FontSize',11, 'Location', 'Northeast');
