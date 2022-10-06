% Generation of Figs. 20-23 in the paper
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

senPos = [
    % minimum number of sensors
    10.23    38.38    16.29
    46.64   -87.12    62.94
    124.02  -7.98     81.16
    105.02  -51.72    26.47
    -81.56    104.48  -80.49
    ]';

% source direction
theta = 22.13*pi/180;
phi = 14.41*pi/180;

[N,M] = size(senPos);

mon = 1000;

% figure(1);plot(souLoc(1), souLoc(2), 'x','markersize',8,'linewidth',1);hold on;grid on;

% Monte-Carlo Simulation

aveNse = 0;
for l = 1:mon
    nse(:,l) = randn(M-1,1);
    err(:,:,l) = randn(N,M);
end
nse = nse - mean(nse,2);
err = err - mean(err,3);

aa = [1,3,7,10,4,1,9,7,2,1,3];
SS = kron(diag(aa(1:M)),eye(N));

disp('Simulation is running ...');

nsePwr = -10;
errLvl = -40:5:20;
souRange = 15*1e2;


NumAlg = 3; % number of compared algorithms
tProc = zeros(1,NumAlg);

for ir = 1:length(souRange)
    disp(['Range: ',num2str(souRange(ir)),'m, ',num2str(ir),'/',num2str(length(souRange)),' ...']);

    %******** Generate Data ********
    % source location
    souLoc = souRange(ir) * [cos(theta)*cos(phi); sin(theta)*cos(phi); sin(phi)] + senPos(:,1);
    % true range
    r = sqrt(sum((souLoc-senPos).^2,1))';
    % true TDOAs
    rd = r(2:end) - r(1);
    g = 1/r(1);
    u0 = (souLoc-senPos(:,1))/r(1);

    [uTh,uPh,ug,ep,er] = deal(zeros(length(nsePwr),mon,NumAlg));
    for is = 1:length(errLvl)
        disp(['10log(Error Level): ',num2str(errLvl(is)),', ',num2str(is),'/',num2str(length(errLvl)),' ...']);

        %     Q = eye(M-1)*sigma(k)^2;
        Qr = 10^(nsePwr/10) * (ones(M-1, M-1)+eye(M-1))/2;
        Qs = 10^(errLvl(is)/10) * SS;
        Qsm = 10^(errLvl(is)/10) * diag(aa(1:M));

        b = sqrt(sum((u0+g*(senPos(:,1)-senPos(:,2:end))).^2,1))';
        B = -diag(b);
        for i = 1:M-1
            C(i,1:N) = -u0';
            C(i,(1:N)+N*i) = (u0 - (senPos(:,i+1)-senPos(:,1))*g)';
        end
        Q = Qr + B\C*Qs*C'/B;

        % calculate CRLB
        CRB = ConsCRLB( senPos, souLoc, Q );
        CRLB_a(is,ir) = CRB(1,1)+CRB(2,2);
        CRLB_g(is,ir) = CRB(3,3);

        % position and DOA estimation
        rng('default');

        bia_p = zeros(N,mon,NumAlg);
        for i = 1:mon
            % measured TDOAs
            rd_m = rd + sqrtm(Qr)*nse(:,i);
            senPos_m = senPos + err(:,:,i)*sqrtm(Qsm);

            nAg = 0;
            % SCO-MPR Method
            nAg = nAg + 1;
            [mprSol, ~] = TDOA_SCO_MPR( senPos_m, rd_m, Qr, Qs );
            uTh(i,nAg) = mprSol(1);
            uPh(i,nAg) = mprSol(2);
            ug(i,nAg) = mprSol(3);

            % SUM-MPR Method
            nAg = nAg + 1;
            [Th2, Ph2, g2, ~] = TDOA_SUM_MPR( senPos_m, rd_m, Qr, Qs );
            uTh(i,nAg) = Th2;
            uPh(i,nAg) = Ph2;
            ug(i,nAg) = g2;

            % GTRS-MPR Method
            nAg = nAg + 1;
            [Th3, Ph3, g3, ~] = TDOA_GTRS_MPR( senPos_m, rd_m, Qr, Qs );
            uTh(i,nAg) = Th3;
            uPh(i,nAg) = Ph3;
            ug(i,nAg) = g3;
        end

        % calculate MSE and bias
        for ia = 1:nAg
            mse_a(is,ir,ia) = mean((uTh(:,ia)-theta).^2+(uPh(:,ia)-phi).^2);
            mse_g(is,ir,ia) = mean((ug(:,ia)-g).^2);

            avBia_a(is,ir,ia) = sqrt(abs(mean(uTh(:,ia))-theta).^2+abs(mean(uPh(:,ia))-phi).^2);
            avBia_g(is,ir,ia) = abs(mean(ug(:,ia))-g);
        end
    end
end

symbs = ['o','v','s','*','^','+','x'];
name = {'SCO-MPR','SUM-MPR','GTRS-MPR'};

xlabtext = '10log(\eta^2(m^2))';
xdata = errLvl;
yl_mse = [-50, 10;
    -90,-20];
yl_bias = [-100,-20;
    -180,-40];
% mstr = 'Err';


% plot results
figure;
for ia = 1:nAg
    plot(xdata, 10*log10(mse_a(:,:,ia)), symbs(ia), 'LineWidth', 1.5, 'DisplayName', name{ia});hold on;grid on;
end
plot(xdata, 10*log10(CRLB_a), '-', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
xlabel(xlabtext, 'FontSize', 13);
ylabel('10log(MSE(\theta,\phi)(rad^2))', 'FontSize', 13);
lgd11 = legend('Show');
set(lgd11, 'FontSize',11, 'Location', 'Northwest');
ylim(yl_mse(1,:));
set(gcf,'Position',[404 310 560 300]);

figure;
for ia = 1:nAg
    plot(xdata, 10*log10(mse_g(:,:,ia)), symbs(ia), 'LineWidth', 1.5, 'DisplayName', name{ia});hold on;grid on;
end
plot(xdata, 10*log10(CRLB_g), '-', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
xlabel(xlabtext, 'FontSize', 13);
ylabel('10log(MSE(g)(1/m^2))', 'FontSize', 13);
lgd2 = legend('Show');
set(lgd2, 'FontSize',11, 'Location', 'Northwest');
ylim(yl_mse(2,:));
set(gcf,'Position',[404 310 560 300]);

% Bias
figure;
for ia = 1:nAg
    plot(xdata, 20*log10(avBia_a(:,:,ia)), symbs(ia), 'LineWidth', 1.5, 'DisplayName', name{ia});hold on;grid on;
end
xlabel(xlabtext, 'FontSize', 13);
ylabel('20log(Bias(\theta,\phi)(rad))', 'FontSize', 13);
h3 = legend('Show');
set(h3, 'FontSize',11, 'Location', 'Northwest');
ylim([-100 20]);
set(gcf,'Position',[404 310 560 300]);

figure;
for ia = 1:nAg
    plot(xdata, 20*log10(avBia_g(:,:,ia)), symbs(ia), 'LineWidth', 1.5, 'DisplayName', name{ia});hold on;grid on;
end
xlabel(xlabtext, 'FontSize', 13);
ylabel('20log(Bias(g)(m^{-1}))', 'FontSize', 13);
h3 = legend('Show');
set(h3, 'FontSize',11, 'Location', 'Southeast');
ylim(yl_bias(2,:));
set(gcf,'Position',[404 310 560 300]);
