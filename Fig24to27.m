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

% target direction
theta = 22.13*pi/180;
phi = 14.41*pi/180;

[N,M] = size(senPos);

mon = 1000;
for l=1:mon
    tmpNsed(:,l) = randn(M,1);
end
aveNse = mean(tmpNsed,2)/sqrt(2);
PP = aveNse(2:end) - aveNse(1);

for l=1:mon
    tmpMp(:,l) = randn(M,1);
end


%% Monte-Carlo Simulation

disp('Simulation is running ...');

% ******* vs. noise power config *******
nsePwr = -70:10:20;
souRange = 15*1e2;

NumAlg = 1; % number of compared algorithms

for ir = 1:length(souRange)
    disp(['Range: ',num2str(souRange(ir)),', ',num2str(ir),'/',num2str(length(souRange)),' ...']);
    %******** Generate Data ********
    % source location
    souLoc = souRange(ir) * [cos(theta)*cos(phi); sin(theta)*cos(phi); sin(phi)] + senPos(:,1);
    % true range
    r = sqrt(sum((souLoc-senPos).^2,1))';
    % true TDOAs
    rd = r(2:end) - r(1);
    g = 1/r(1);

    [uTh,uPh,ug,ep,er] = deal(zeros(length(nsePwr),mon,NumAlg));
    for in = 1:length(nsePwr)
        disp(['Noise power: ',num2str(nsePwr(in)),', ',num2str(in),'/',num2str(length(nsePwr)),' ...']);

        % Q = eye(M-1)*sigma(k)^2;
        Q = 10^(nsePwr(in)/10) * (ones(M-1, M-1)+eye(M-1))/2;

        for ig = 1:3
            gamma = 1 - 0.2*(ig-1);
            % position and DOA estimation
            rng('default');

            bia_p = zeros(N,mon,NumAlg);
            for i = 1:mon
                % measured TDOAs
                mu = rand(M,1)*20*sqrt(10^(nsePwr(in)/10));
                gg = ones(M,1);
                gg((rand(M,1)+gamma-1)>0) = 0;
                rMp = sqrt(10^(nsePwr(in)/10))*(tmpMp(:,i)+mu).*gg;
                rdMp = (rMp(2:M)-rMp(1))/sqrt(2);
                rdNse = sqrt(10^(nsePwr(in)/10)) * ((tmpNsed(2:M,i)-tmpNsed(1,i))/sqrt(2)-PP);
                rd_m = rd + rdNse + rdMp;

                nAg = 0;
                % SCO-MPR Method
                nAg = nAg + 1;
                [mprSol, ~] = TDOA_SCO_MPR( senPos, rd_m, Q );
                uTh(i,nAg) = mprSol(1);
                uPh(i,nAg) = mprSol(2);
                ug(i,nAg) = mprSol(3);
            end

            % calculate MSE and bias
            for ia = 1:nAg
                mse_a(in,ir,ig,ia) = mean((uTh(:,ia)-theta).^2+(uPh(:,ia)-phi).^2);
                mse_g(in,ir,ig,ia) = mean((ug(:,ia)-g).^2);

                avBia_a(in,ir,ig,ia) = sqrt(abs(mean(uTh(:,ia))-theta).^2+abs(mean(uPh(:,ia))-phi).^2);
                avBia_g(in,ir,ig,ia) = abs(mean(ug(:,ia))-g);
            end
        end
    end
end

symbs = ['o','^','*'];

xlabtext = '10log(\sigma^2(m^2))';
xdata = nsePwr;
fileName = ['mat_Ref_TDOA_Noise_',datestr(now,'yyyymmmdd_HHMM'),'.mat'];
yl_mse = [-60, 20;
    -100,0];
yl_bias = [-120,20;
    -160,0];


%% plot results
figure;
for ia = 1:nAg
    for ig = 1:3
        plot(xdata, 10*log10(mse_a(:,:,ig,ia)), symbs(ig), 'LineWidth', 1.5, ...
            'DisplayName', ['SCO-MPR, \gamma_i = ',num2str(1 - 0.2*(ig-1))]);hold on;grid on;
    end
end
xlabel(xlabtext, 'FontSize', 13);
ylabel('10log(MSE(\theta,\phi)(rad^2))', 'FontSize', 13);
lgd11 = legend('Show');
set(lgd11, 'FontSize',11, 'Location', 'Northwest');
set(gcf,'Position',[404 310 560 300]);

figure;
for ia = 1:nAg
    for ig = 1:3
        plot(xdata, 10*log10(mse_g(:,:,ig,ia)), symbs(ig), 'LineWidth', 1.5, ...
            'DisplayName', ['SCO-MPR, \gamma_i = ',num2str(1 - 0.2*(ig-1))]);hold on;grid on;
    end
end
xlabel(xlabtext, 'FontSize', 13);
ylabel('10log(MSE(g)(1/m^2))', 'FontSize', 13);
lgd2 = legend('Show');
set(lgd2, 'FontSize',11, 'Location', 'Northwest');
set(gcf,'Position',[404 310 560 300]);

% bias
figure;
for ia = 1:nAg
    for ig = 1:3
        plot(xdata, 20*log10(avBia_a(:,:,ig,ia)), symbs(ig), 'LineWidth', 1.5, ...
            'DisplayName', ['SCO-MPR, \gamma_i = ',num2str(1 - 0.2*(ig-1))]);hold on;grid on;
    end
end
xlabel(xlabtext, 'FontSize', 13);
ylabel('20log(Bias(\theta,\phi)(rad))', 'FontSize', 13);
h3 = legend('Show');
set(h3, 'FontSize',11, 'Location', 'Northwest');
xlim([min(xdata),max(xdata)]); ylim([-220 10]);
set(gcf,'Position',[404 310 560 300]);

figure;
for ia = 1:nAg
    for ig = 1:3
        plot(xdata, 20*log10(avBia_g(:,:,ig,ia)), symbs(ig), 'LineWidth', 1.5, ...
            'DisplayName', ['SCO-MPR, \gamma_i = ',num2str(1 - 0.2*(ig-1))]);hold on;grid on;
    end
end
xlabel(xlabtext, 'FontSize', 13);
ylabel('20log(Bias(g)(m^{-1}))', 'FontSize', 13);
h3 = legend('Show');
set(h3, 'FontSize',11, 'Location', 'Northwest');
set(gcf,'Position',[404 310 560 300]);
