function [tsig,teps] = rotmat(theta)
% Usage: [tsig,teps] = rotmat(theta)
%
% Purpose: Compute rotation matrices for stress and strain vectors
%
% Input: theta -- rotation angle (degrees)
%
% Output: tsig -- rotation matrix for stress vectors
%         teps -- rotation matrix for strain vectors

    theta = theta*pi/180;
    m = cos(theta);
    n = sin(theta);
    
    tsig = [m*m n*n  2*m*n;...
            n*n m*m -2*m*n;...
           -m*n m*n   m*m-n*n];
        
    teps = [m*m     n*n  m*n ;...
            n*n     m*m -m*n;...
           -2*m*n 2*m*n  m*m-n*n];
end 