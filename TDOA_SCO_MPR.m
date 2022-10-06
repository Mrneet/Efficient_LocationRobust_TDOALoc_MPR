function [mprSol,pos] = TDOA_SCO_MPR(senPos, rd, varargin)
% [mprSol,pos] = TDOA_SCO_MPR(senPos, rd, varargin)
%
% Closed-form method for localization in MPR using successive constrained 
% optimization, s1 not necessarily in the origin.
%
% Input:
%   senPos:     (Dim x M), postions of reciveing sensors, each column is a sensor position 
%               and the first column is the reference sensor location for TDOA.
%   rd:         ((M-1)x1), TDOA measurement vector.
%   Q:          ((M-1)x(M-1)), covariance matrix of TDOAs.
%
% Output:
%   varargout: including
%     **Dim == 2:
%       theta:	(1x2), DOA estimation, including 1st & 2nd stage.
%       g:      (1x2), inverse-range (g) estimation, including 1st & 2nd stage.
%     **Dim == 3:
%       theta:	(1x2), azimuth estimation, including 1st & 2nd stage.
%       phi:	(1x2), elevation estimation, including 1st & 2nd stage.
%       g:      (1x2), g estimation, including 1st & 2nd stage.
%       pos:    (1xN), source position, only 2nd stage.
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
v1 = 0.5*(rd.^2-sum((senPos(:,2:end)-senPos(:,1))'.^2,2));

% first stage
G1 = (senPos(:,2:end)-senPos(:,1))';
W1 = inv(Qr);
for itNum = 1:2
    iG1 = (G1'*W1*G1)\(G1')*W1;
    alpha = iG1*h1;
    beta = -iG1*v1;
    
    x(1) = beta'*beta;%(beta-senPos(:,1))'*(beta-senPos(:,1));
    x(2) = 2*beta'*alpha;%2*(beta-senPos(:,1))'*alpha;
    x(3) = alpha'*alpha - 1;
    root = [-x(2)+sqrt(x(2)^2-4*x(1)*x(3)),-x(2)-sqrt(x(2)^2-4*x(1)*x(3))]/(2*x(1));
    
    % delete complex roots
    reRoot = root(imag(root)==0);
    L = length(reRoot);
    % guarantee that Y is not empty
    if L == 0
        [~,I] = min(imag(root));
        reRoot = real(root(I));
        L = 1;
    end
    
    % find property u_bar and g
    u_tmp = zeros(N,L);J=zeros(L,1);
    for i = 1:L
        u_tmp(:,i) = alpha + beta*reRoot(i);
        r_rec = sqrt(sum((u_tmp(:,i)-reRoot(i)*(senPos-senPos(:,1))).^2,1))'/reRoot(i);
        rd_rec = r_rec(2:end) - r_rec(1);
        J(i) = (rd-rd_rec)'/Qr*(rd-rd_rec);
    end
    [~,ind] = min(J);
    ubar1 = u_tmp(:,ind);
    g1 = reRoot(ind);
        
    h2 = -rd - G1*ubar1 - v1*g1;
    G3 = -G1;
    v2 = -v1;
    O = W1 - W1*v2/(v2'*W1*v2)*v2'*W1;
    iG2 = eye(N)/(G3'*O*G3);
    lambda = - (ubar1'*iG2*G3'*O*h2) / (ubar1'*iG2*ubar1);
    delta_u = iG2*(G3'*O*h2+lambda*ubar1);
    delta_g = (v2'*W1*v2)\v2'*W1*(h2-G3*delta_u);

    psi2 = [ubar1;g1] - [delta_u;delta_g];

    % update weighting matrix
    b = sqrt(sum((psi2(1:N)+psi2(end)*(senPos(:,1)-senPos(:,2:end))).^2,1))';
    B1 = -diag(b);
    for i = 1:M-1
        C(i,1:N) = -psi2(1:N)';
        C1(i,(1:N)+N*i) = (psi2(1:N) - (senPos(:,i+1)-senPos(:,1))*psi2(end))';
    end
    W1 = inv(B1*Qr*B1 + C1*Qs*C1');
end

psi2 = sign(real(psi2)).*abs(psi2);
mprSol = [atan2( psi2(2), psi2(1) );
          atan2( psi2(3), norm(psi2(1:2),2));
          psi2(N+1)];
pos = psi2(1:N)/psi2(N+1) + senPos(:,1);


% B2 = diag([2*psi2(1:N);1]);
% h2 = [psi2(1:N-1).^2;psi2(N)^2-1;psi2(end)];
% G2 = [eye(N-1),zeros(N-1,1);-ones(1,N-1),0;zeros(1,N-1),1];
% W2 = B2\(G3'*W1*G3)/B2;
% Psi0 = ((G2'*W2*G2)\(G2')*W2*h2);  % solution of 2nd stage
% Psi2 = sign(real(Psi0)).*abs(Psi0);  % to keep it real
% 
% mprSol(1,1) = atan2( abs(sqrt(Psi2(2)))*sign(psi2(2)), abs(sqrt(Psi2(1)))*sign(psi2(1)) );
% mprSol(2,1) = atan2( abs(sqrt(1-sum(Psi2(1:2))))*sign(psi2(3)), abs(sqrt(sum(Psi2(1:2)))) );
% mprSol(3,1) = Psi2(end);
% pos = [cos(mprSol(1))*cos(mprSol(2));sin(mprSol(1))*cos(mprSol(2));sin(mprSol(2))]/mprSol(3) + senPos(:,1);

end

