function [ varargout ] = TDOA_GTRS_MPR( senPos, rd, varargin )
% 
% Estimation of the source location in MPR by the GTRS method, s1 not 
% necessarily in the origin.
%
% Input:
%   senPos:	    (Dim x M), postions of reciveing sensors, each column is a sensor position
%               and the first column is the reference sensor location for TDOA.
%   rd:         ((M-1)x1), TDOA measurement vector.
%   Q:          ((M-1)x(M-1)), covariance matrix of TDOAs.
%
% Output:
%   varargout: containing
%     **Dim == 2:
%       theta:	(1x2), DOA estimation, including 1st & 2nd stages.
%       g:      (1x2), inverse-range (g) estimation, including 1st & 2nd stages.
%       pos:	(1x2), source position, only 2nd stage.
%     **Dim == 3:
%       theta:	(1x2), azimuth estimation, including 1st & 2nd stage.
%       phi:	(1x2), elevation estimation, including 1st & 2nd stage.
%       g:      (1x2), g estimation, include 1st & 2nd stage.
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
B1 = eye(M-1);
W = eye(M-1)/(B1*Qr*B1);
% W = eye(M-1);

for iter = 1:2
    A = G1(:,1:N);
    a = G1(:,N+1);
    O = W-W*a/(a'*W*a)*a'*W;
    AOA = A'*O*A;
    S = eye(N);
    
    [U,D,~] = svd(AOA);
    r = diag(D);
    k = U'*A'*O*h1;
    
    if N == 2
        x(1) = 1;
        x(2) = 2*r(1)+2*r(2);
        x(3) = r(1)^2+2*r(1)*r(2)+r(2)^2-k(1)^2-k(2)^2;
        x(4) = 2*r(1)^2*r(2)+2*r(1)*r(2)^2-2*k(1)^2*r(2)-2*k(2)^2*r(1);
        x(5) = r(1)^2*r(2)^2-k(1)^2*r(2)^2-k(2)^2*r(1)^2;
    elseif N == 3
        x(1) = 1;
        x(2) = 2*r(3)+2*r(2)+2*r(1);
        x(3) = r(3)^2+4*r(2)*r(3)+4*r(1)*r(3)+r(2)^2+4*r(1)*r(2)+r(1)^2-k(3)^2-k(2)^2-k(1)^2;
        x(4) = 2*r(2)*r(3)^2+2*r(1)*r(3)^2+2*r(2)^2*r(3)+8*r(1)*r(2)*r(3)+2*r(1)^2*r(3)-2*k(2)^2*r(3)-2*k(1)^2*r(3)+2*r(1)*r(2)^2+2*r(1)^2*r(2)-2*k(3)^2*r(2)-2*k(1)^2*r(2)-2*k(3)^2*r(1)-2*k(2)^2*r(1);
        x(5) = r(2)^2*r(3)^2+4*r(1)*r(2)*r(3)^2+r(1)^2*r(3)^2-k(2)^2*r(3)^2-k(1)^2*r(3)^2+4*r(1)*r(2)^2*r(3)+4*r(1)^2*r(2)*r(3)-4*k(1)^2*r(2)*r(3)-4*k(2)^2*r(1)*r(3)+r(1)^2*r(2)^2-k(3)^2*r(2)^2-k(1)^2*r(2)^2-4*k(3)^2*r(1)*r(2)-k(3)^2*r(1)^2-k(2)^2*r(1)^2;
        x(6) = 2*r(1)*r(2)^2*r(3)^2+2*r(1)^2*r(2)*r(3)^2-2*k(1)^2*r(2)*r(3)^2-2*k(2)^2*r(1)*r(3)^2+2*r(1)^2*r(2)^2*r(3)-2*k(1)^2*r(2)^2*r(3)-2*k(2)^2*r(1)^2*r(3)-2*k(3)^2*r(1)*r(2)^2-2*k(3)^2*r(1)^2*r(2);
        x(7) = r(1)^2*r(2)^2*r(3)^2-k(1)^2*r(2)^2*r(3)^2-k(2)^2*r(1)^2*r(3)^2-k(3)^2*r(1)^2*r(2)^2;
    end
    root = roots(x);
    
    % delete complex roots
    reRoot = root(imag(root)==0);
    L = length(reRoot);
    if L == 0 % to insure that Y is not empty
        [~,I] = min(imag(root));
        reRoot = real(root(I));
        L = 1;
    end
    
    Y = zeros(N+1,L);J = zeros(1,L);
    for i = 1:L
        Y(1:N,i) = (AOA+reRoot(i)*S)\A'*O*h1;
        Y(N+1,i) = (a'*W*a)\a'*W*(h1-A*Y(1:N,i));
        J(i) = (h1-G1*Y(:,i))'*W*(h1-G1*Y(:,i));
    end
    [~,ind] = min(J);
    Phi0 = Y(:,ind);
    Phi = sign(real(Phi0)).*abs(Phi0);   % to keep it real
    
    b = 1 + l.^2*Phi(end)^2 - 2*Phi(end)*l.*(s_bar'*Phi(1:N));
    B1 = -diag(sqrt(b));
    for i = 1:M-1
        T(i,(1:N)+N*(i-1)) = (Phi(1:N) - (senPos(:,i+1)-senPos(:,1))*Phi(end))';
    end
    C1 = -T*kron([-ones(M-1,1),eye(M-1)],eye(N));
    W = eye(M-1)/(B1*Qr*B1 + C1*Qs*C1');
end

if N == 2
    theta = atan2(Phi(2),Phi(1));
    g = Phi(3);
    pos = [cos(theta);sin(theta)]/g + senPos(:,1);
    
    varargout{1} = theta;
    varargout{2} = g;
    varargout{3} = pos;
elseif N == 3
    theta = atan2(Phi(2),Phi(1));
    phi = atan2(Phi(3),sqrt(Phi(1)^2+Phi(2)^2));
    g = Phi(4);
    pos = [cos(theta)*cos(phi);sin(theta)*cos(phi);sin(phi)]/g + senPos(:,1);
    
    varargout{1} = theta;
    varargout{2} = phi;
    varargout{3} = g;
    varargout{4} = pos;
else
    error('Please check your input format of sensor positions');
end

end
