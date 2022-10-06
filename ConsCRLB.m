function CRB = ConsCRLB( senPos, srcLoc, Q )
% ConsCRLB( senPos, srcLoc, Q )
% 
% Evaluation of the CRLB.
%
% Input:
%   senPos:	    (Dim x M), postions of reciveing sensors, each column is a sensor position 
%               and the first column is the reference sensor location for TDOA.
%   srcLoc  :     (Dim x 1), source location.
%   Q:          ((M-1)x(M-1)), covariance matrix of TDOAs.
%
% Output:
%   CRB:        CRB matrix     
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
r = sqrt(sum((senPos-srcLoc).^2,1))';
u0 = (srcLoc-senPos(:,1))/r(1);
g0 = 1/r(1);

p = (u0 - g0*(senPos(:,2:end)-senPos(:,1)))./sqrt(sum((u0 - g0*(senPos(:,2:end)-senPos(:,1))).^2,1));
% p = u0 - (senPos(:,2:end)-senPos(:,1))*g0;
P = p./sqrt(sum(p.^2,1));

DF(1:N,:) = -(senPos(:,2:end)-senPos(:,1))./sqrt(sum((u0 - g0*(senPos(:,2:end)-senPos(:,1))).^2,1));
% DF(1:N,:) = -(senPos(:,2:end)-senPos(:,1))./r(2:end)'*r(1);
DF(N+1,:) = (-P'*u0+ones(M-1,1))*r(1)^2;
CRB1 = inv(DF/Q*DF');

S = diag([ones(N,1),;0]);
Phio = [u0;g0];
F = Phio'*S;
CRB2 = CRB1 - CRB1*F'/(F*CRB1*F')*F*CRB1;

if N == 2 % 2D
    doa = atan2(srcLoc(2),srcLoc(1));
    D1 = [-sin(doa),cos(doa),0;
          0,0,1];
elseif N == 3 % 3D
    theta = atan2(srcLoc(2)-senPos(2,1), srcLoc(1)-senPos(1,1));
    phi = atan2(srcLoc(3)-senPos(3,1), norm(srcLoc(1:2)-senPos(1:2,1),'fro'));
    D1 = [-sin(theta)/cos(phi), cos(theta)/cos(phi),    0,          0;
          -cos(theta)*sin(phi), -sin(theta)*sin(phi),   cos(phi),   0;
          0,                    0,                      0,          1];
else
    error('Please check your input format of sensor positions');
end

CRB = D1*CRB2*D1';