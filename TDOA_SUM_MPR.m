function [ varargout ] = TDOA_SUM_MPR( senPos, rd, varargin )
%  [ varargout ] = TDOA_SUM_MPR( senPos, rd, varargin )
%
% Estimation of of the source position by MPR using successive unconstrained 
% minimization, s1 not necessarily in the origin.
%
% Input:
%   senPos:	        (Dim x M), postions of reciveing sensors, each column is a sensor position 
%               and the first column is the reference sensor location for TDOA.
%   rd:         ((M-1)x1), TDOA measurement vector.
%   Q:          ((M-1)x(M-1)), covariance matrix of TDOAs.
%
% Output:
%   varargout: containing
%     **Dim == 2:
%       theta:	(1x2), DOA estimation, including 1st & 2nd stage.
%       g:      (1x2), inverse-range (g) estimation, including 1st & 2nd stage.
%       pos:	(1x2), source position, only 2nd stage.
%     **Dim == 3:
%       theta:	(1x2), azimuth estimation, including 1st & 2nd stage.
%       phi:	(1x2), elevation estimation, including 1st & 2nd stage.
%       g:      (1x2), inverse-range (g) estimation, including 1st & 2nd stage.
%       pos:	(1x3), source position, only 2nd stage.
%
% Reference: Y. Sun, K. C. Ho, G. Wang. J. Chen, Y. Yang, L. Chen, and Q. Wan, 
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

[N,M] = size(senPos);

Qr = varargin{1};
Qs = 0;
if length(varargin) == 2
    Qs = varargin{2};
end

h1 = -rd;
l = sqrt(sum((senPos(:,2:end) - senPos(:,1)).^2, 1))';
s_bar = (senPos(:,2:end)- senPos(:,1))./l';

% first stage
G1 = [(senPos(:,2:end)-senPos(:,1))', 0.5*(rd.^2-sum((senPos(:,2:end)-senPos(:,1))'.^2,2))];
Phi1 = (G1'/Qr*G1)\(G1')/Qr*h1;   % initial estimation

b = 1 + l.^2*Phi1(end)^2 - 2*Phi1(end)*l.*(s_bar'*Phi1(1:N));
B1 = -diag(sqrt(b));
for i = 1:M-1
    T(i,(1:N)+N*(i-1)) = (Phi1(1:N) - (senPos(:,i+1)-senPos(:,1))*Phi1(end))';
end
C1 = -T*kron([-ones(M-1,1),eye(M-1)],eye(N));
W1_inv = B1*Qr*B1 + C1*Qs*C1';
Phi20 = ((G1'/W1_inv*G1)\G1'/W1_inv*h1); % solution of 1st stage
Phi2 = sign(real(Phi20)).*abs(Phi20);   % to keep it real

% second stage
B2 = diag([2*Phi2(1:N);1]);
h2 = [Phi2(1:N-1).^2;Phi2(N)^2-1;Phi2(end)];
G2 = [eye(N-1),zeros(N-1,1);-ones(1,N-1),0;zeros(1,N-1),1];
W2 = B2\(G1'/W1_inv*G1)/B2;
Psi0 = ((G2'*W2*G2)\(G2')*W2*h2);  % solution of 2nd stage
Psi = sign(real(Psi0)).*abs(Psi0);  % to keep it real


if N == 2
    theta = atan2( abs(sqrt(1-Psi(1)))*sign(Phi2(2)), abs(sqrt(Psi(1)))*sign(Phi2(1)) );
    g = Psi(2);
    pos = [cos(theta);sin(theta)]/g + senPos(:,1);

    varargout{1} = theta;
    varargout{2} = g;
    varargout{3} = pos;
elseif N == 3
    theta = atan2( abs(sqrt(Psi(2)))*sign(Phi2(2)), abs(sqrt(Psi(1)))*sign(Phi2(1)) );
    phi = atan2( abs(sqrt(1-sum(Psi(1:2))))*sign(Phi2(3)), abs(sqrt(sum(Psi(1:2)))) );
    g = Psi(end);
    pos = [cos(theta)*cos(phi);sin(theta)*cos(phi);sin(phi)]/g + senPos(:,1);
    
    varargout{1} = theta;
    varargout{2} = phi;
    varargout{3} = g;
    varargout{4} = pos;
else
    error('Please check your input format of sensor positions');
end

