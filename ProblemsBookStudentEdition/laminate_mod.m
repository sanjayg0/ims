function [SF] = laminate_mod(moduli,layup,filename,varargin)
% Usage: [SF] = laminate_mod(moduli,layup,filename,<loads>,<failure>)
%
% Purpose: Compute moduli and compliance matrices for
%          a laminate, evaluate failure criteria
%
% Input: moduli   -- structure with lamina moduli (SI)
%        filename -- output file name
%        layup    -- structure with ply thickness (SI) and code (degrees)
%        loads    -- (optional) structure with in-plane and bending loads (SI)
%        failure  -- (optional) structure with failure criteria (SI)
%
% Output: computations written to filename
%         SF      -- Tsai-Wu safety factor

    % Unpack laminate information
    % Elastic moduli
    E1   = moduli.e1;
    E2   = moduli.e2;
    NU12 = moduli.nu12;
    G12  = moduli.g12;  
    
    % Laminate layout details
    h0   = layup.h0;     % Ply thickness
    code = layup.code;   % Laminate code (in degrees)
    N    = length(code); % Number of plies
    
    if nargin > 3
        % Set load details
        Nx  = varargin{1}.nx;
        Ny  = varargin{1}.ny;
        Nxy = varargin{1}.nxy;
        Mx  = varargin{1}.mx;
        My  = varargin{1}.my;
        Mxy = varargin{1}.mxy;
    end

    if nargin == 5
        % Failure parameters
        F1t = varargin{2}.f1t;
        F1c = varargin{2}.f1c;
        F2t = varargin{2}.f2t;
        F2c = varargin{2}.f2c;
        F12 = varargin{2}.f12;
    end

    % Open output file
    fid      = fopen(filename,'w');

    % Write input data to file
    sep(1:63) = '*'; sep2(1:63)='-';

    fprintf(fid,'INPUT\n%s\n',sep);
    
    fprintf(fid,'Lamina properties (Moduli in Pa. Thickness in m):\n%s\n',sep);
    fprintf(fid,'E11=%5.3e\t E22=%5.3e\t NU12=%5.3e\t G12=%5.3e\n\n',...
        E1,E2,NU12,G12);
    
    if (nargin == 5)
        fprintf(fid,'F1t=%5.3e\t F1c=%5.3e\t F2t=%5.3e\t F2c=%5.3e\n',...
            F1t,F1c,F2t,F2c);
        fprintf(fid,'F12=%5.3e\n\n',F12);
    end
    
    fprintf(fid,'Lamina thickness: H0 = %5.3e\n',h0);
    fprintf(fid,'Number of plies: N = %2i\n',N);
    fprintf(fid,'Laminate code:\t');
    fprintf(fid,'%+2i,',code);
    fprintf(fid,'\n');
    
    if (nargin > 3)
        fprintf(fid,'%s\nLaminate load N_{ij} in N/m, M_{ij} in N-m/m:\n%s\n',...
            sep,sep);
        fprintf(fid,'Nx = % 5.3e\t Ny = % 5.3e\t Nxy = % 5.3e\n',Nx,Ny,Nxy);
        fprintf(fid,'Mx = % 5.3e\t My = % 5.3e\t Mxy = % 5.3e\n',Mx,My,Mxy);
    end
    fprintf(fid,'%s\n',sep);

    % Start calculation
    % Laminate thickness
    h = h0*N;

    % On axis compliance and stiffness matrices from properties
    Scomp(:,:) = zeros(3,3);
    Scomp(1,1) = 1.0/E1;
    Scomp(2,2) = 1.0/E2;
    Scomp(3,3) = 1.0/G12;
    Scomp(2,1) = -NU12/E1;
    Scomp(1,2) = Scomp(2,1);
    Q          = inv(Scomp);

    % Compute z-coordinates for various plies in the laminate
    z(1) = -0.5*h;
    for k = 1:N
        z(k+1) = z(k) + h0;
    end

    % Intitialize laminate stiffnesses:   
    %  [A] -- In-plane stiffness matrix
    %  [D] -- Bending stiffness matrix
    %  [B] -- Bending-stretching coupling matrix
    A(:,:) =  zeros(3,3);
    B(:,:) =  zeros(3,3);
    D(:,:) =  zeros(3,3);

    % Loop over plies and sum contributions to the stiffnesses
    for k = 1:N    
        % Compute lamina stiffness in (x,y) coordinate frame
        [Tsig Teps] = rotmat(code(k));
        Qbar        = inv(Tsig)*Q*Teps;

        % Compute components of [A],[B],[D]  
        A = A + Qbar*(z(k+1)-z(k));
        B = B + Qbar*(z(k+1)^2 - z(k)^2);
        D = D + Qbar*(z(k+1)^3 - z(k)^3);
    end

    % Adjust pre-factors
    B = 0.5*B;
    D = D/3.0;

    % Print Output
    fprintf(fid,'\nOUTPUT');

    % Write lamina stiffness matrices  
    out3(A*1e-6,'Extensional Stiffness Matrix [A], MN/m',fid);
    out3(B*1e-3,'Coupling Matrix [B], kN',fid);
    out3(D,     'Bending Stiffness Matrix [D], N-m',fid);

    % Check if [B] is zero-valued
    isymm = sum(sum( abs(B) > 1e-9 ));
    fprintf(fid,'\n\n%s\n',sep2);
    switch (isymm)
        case 0
            fprintf(fid,'%s\n','LAMINATE IS SYMMETRIC');
        otherwise    
            fprintf(fid,'%s\n','NOT A SYMMETRIC LAMINATE');
    end
    fprintf(fid,'%s\n',sep2);

    % Compute and write compliances useful for symmetric laminates
    Ainv = inv(A);
    Dinv = inv(D);
    out3(Ainv*1e9,'Extensional Compliance Matrix [A]^{-1}, (GN/m)^{-1}',fid);
    out3(Dinv*1e3,'Bending Compliance Matrix [D]^{-1}, (kNm)^{-1}',fid);

    % Write in-plane laminate stiffness and compliance
    out3(1e-9*A/h,'In-plane Laminate Stiffness Matrix [A]/h, (GPa)',fid);
    out3(1e12*h*inv(A),...
        'In-plane Laminate Compliance Matrix h[A]^{-1}, (TPa)^{-1}',fid);

    % Determine in-plane laminate engineering constants
    Astari   =  h*Ainv;
    Ebarx    =  1.0e-9/Astari(1,1);
    nubarxy  = -Astari(2,1)/Astari(1,1);
    etabarxs =  Astari(3,1)/Astari(1,1);

    Ebary    =  1.0e-9/Astari(2,2);
    nubaryx  = -Astari(1,2)/Astari(2,2);
    etabarys =  Astari(3,2)/Astari(2,2);

    Gbarxy   =  1.0e-9/Astari(3,3);
    etabarsx =  Astari(1,3)/Astari(3,3);
    etabarsy =  Astari(2,3)/Astari(3,3);

    fprintf(fid,'\n%s\nIn-plane Laminate Engineering Constants\n%s\n',sep,sep);

    fprintf(fid,'Ebarx    = % 4.3e GPa\n', Ebarx);
    fprintf(fid,'nubarxy  = % 4.3e\n',   nubarxy);
    fprintf(fid,'etabarxs = % 4.3e\n\n',etabarxs);

    fprintf(fid,'Ebary    = % 4.3e GPa\n', Ebary);
    fprintf(fid,'nubaryx  = % 4.3e\n',   nubaryx);
    fprintf(fid,'etabarys = % 4.3e\n\n',etabarys);

    fprintf(fid,'Gbarxy   = % 4.3e GPa\n',Gbarxy);
    fprintf(fid,'etabarsx = % 4.3e\n',  etabarsx);
    fprintf(fid,'etabarsy = % 4.3e\n',  etabarsy);


    % Write laminate face bending stiffness and compliance
    out3(1e-9*D*12/h^3,'Effective Laminate Face Stiffness (12/h^3)[D], (GPa)',fid);
    out3(1e12*Dinv*h^3/12,...
        'Effective Laminate Face Compliance (h^3/12)[D]^{-1}, (TPa)^{-1}',fid);

    Dstari = Dinv*h^3/12; 
    Ebarxf    =  1.0e-9/Dstari(1,1);
    nubarxyf  = -Dstari(2,1)/Dstari(1,1);
    etabarxsf =  Dstari(3,1)/Dstari(1,1);

    Ebaryf    =  1.0e-9/Dstari(2,2);
    nubaryxf  = -Dstari(1,2)/Dstari(2,2);
    etabarysf =  Dstari(3,2)/Dstari(2,2);

    Gbarxyf   =  1.e-9/Dstari(3,3);
    etabarsxf =  Dstari(1,3)/Dstari(3,3);
    etabarsyf =  Dstari(2,3)/Dstari(3,3);

    fprintf(fid,'\n%s\nLaminate Flexural Engineering Constants\n%s\n',sep,sep);

    fprintf(fid,'Ebarxf    = % 4.3e GPa\n', Ebarxf);
    fprintf(fid,'nubarxyf  = % 4.3e\n',   nubarxyf);
    fprintf(fid,'etabarxsf = % 4.3e\n\n',etabarxsf);

    fprintf(fid,'Ebaryf    = % 4.3e GPa\n', Ebaryf);
    fprintf(fid,'nubaryxf  = % 4.3e\n',   nubaryxf);
    fprintf(fid,'etabarysf = % 4.3e\n\n',etabarysf);

    fprintf(fid,'Gbarxyf   = % 4.3e GPa\n',Gbarxyf);
    fprintf(fid,'etabarsxf = % 4.3e\n',  etabarsxf);
    fprintf(fid,'etabarsyf = % 4.3e\n\n',etabarsyf);

    % Compute and write laminate compliances
    switch (isymm)
        case 0     
            Aprime = Ainv;       % [A']        
            Bprime = zeros(3,3); % [B']       
            Cprime = zeros(3,3); % [C']   
            Dprime = Dinv;       % [D']
        otherwise                  
            AUX1 = Ainv*B;                       
            AUX2 = B*Ainv;                      
            Dprime =  inv(D-B*Ainv*B);       % [D']    
            Aprime =  Ainv+AUX1*Dprime*AUX2; % [A']    
            Bprime = -AUX1*Dprime;           % [B']    
            Cprime = -Dprime*AUX2;           % [C']    
    end

    out3(1e9*Aprime,"Aprime matrix [A'], m/GN" ,fid);
    out3(1e6*Bprime,"Bprime matrix [B'], 1/MN" ,fid);
    out3(1e6*Cprime,"Cprime matrix [C'], 1/MN" ,fid);
    out3(1e3*Dprime,"Dprime matrix [D'], 1/kNm",fid);

    if (nargin > 3)
        % Laminate mid-plane strain and curvature
        Nvec = [Nx Ny Nxy]';
        Mvec = [Mx My Mxy]';

        epsv  = Aprime*Nvec + Bprime*Mvec;
        kappa = Cprime*Nvec + Dprime*Mvec;

        fprintf(fid,'\n\n%s\nMid-plane strain components\n%s\n',sep,sep);
        fprintf(fid,'eps_x    = % 4.3e\n',  epsv(1));
        fprintf(fid,'eps_y    = % 4.3e\n',  epsv(2));
        fprintf(fid,'gamma_xy = % 4.3e\n\n',epsv(3));

        fprintf(fid,'%s\nMid-plane curvature components\n%s\n',sep,sep);
        fprintf(fid,'kappa_x  = % 4.3e\n',  kappa(1));
        fprintf(fid,'kappa_y  = % 4.3e\n',  kappa(2));
        fprintf(fid,'kappa_xy = % 4.3e\n\n',kappa(3));

        % Compute ply level stresses
        sigply = zeros(3,N);
        for k=1:N
            % Compute layer stiffness in (x,y) coordinate frame
            [Tsig,Teps] = rotmat(code(k));
            Qbar = inv(Tsig)*Q*Teps;

            % Compute middle of ply stresses
            zmean    = (z(k)+z(k+1))/2;
            epsmid   = epsv + kappa*zmean;
            siglocal = Qbar*epsmid;

            % Compute stress and strain in ply coordinates
            epslocal_ply = Teps*epsmid;
            sigply(:,k)  = Q*epslocal_ply;
        end


        % Compute stresses and failure criteria ply by ply
        fprintf(fid,'\n%s\n',sep);
        fprintf(fid,...
            'Compute lamina stresses, (opt) failure criterion:\n%s\n',sep);
        Sfactor = zeros(N,1);
        for k=1:N
            % Extract ply stresses in lamina coordinate frame
            sig1  = sigply(1,k);
            sig2  = sigply(2,k);
            sig12 = sigply(3,k);

            fprintf(fid,'Layer %u:\n',k);
            fprintf(fid,'sigma_1\t\t sigma_2\t sigma_12 in MPa:\n');
            fprintf(fid,'%+4.3e\t %+4.3e\t %+4.3e\n',...
                sig1/1.e6,sig2/1.e6,sig12/1.e6);

            if (nargin == 5)
                % Check individual criteria
                if  sig1/F1t  >=1
                    fprintf(fid,'sigma1/f1t >= 1, lamina fiber failure\n');
                end
                if abs(sig1)/F1c >=1     
                    fprintf(fid,'abs(sigma1)/f1c >= 1, lamina fiber failure\n');
                end
                if  sig2/F2t >=1 
                    fprintf(fid,'sigma2/f2t >= 1, lamina matrix failure\n');
                end
                if abs(sig2)/F2c >=1 
                    fprintf(fid,'abs(sigma2)/f2c >= 1, lamina matrix failure\n');
                end
                if abs(sig12)/F12 >=1 
                 fprintf(fid,'abs(sigma12)/f12 >= 1, lamina shear failure\n'); 
                end

                % Azzi-Tsai criteria
                AT = (sig1/F1t)^2+(sig2/F2t)^2 +(sig12/F12)^2 - sig1*sig2/(F1t^2);
                if AT >=1
                    fprintf(fid,...
                        'Azzi-Tsai parameter = % 4.3e >=1, lamina has failed\n',...
                        AT);
                end

                % Tsai-Wu criteria
                f1   = (1/F1t) - (1/F1c);
                f2   = (1/F2t) - (1/F2c);
                f11  = 1/(F1t*F1c);
                f22  = 1/(F2t*F2c);
                f66  = 1/(F12*F12);
                f12  = -0.5/sqrt(F1t*F1c*F2t*F2c);
                afac = f11*sig1^2+ f22*sig2^2 + f66*sig12^2 + 2*f12*sig1*sig2;
                bfac = f1*sig1 + f2*sig2;
                TW   = afac + bfac; 
                if TW >=1
                    fprintf(fid,...
                        'Tsai-Wu parameter = % 4.3e >=1, lamina has failed\n',...
                        TW);       
                end

                % Tsai-Wu strength safety factor 
                Sfactor(k) = (-bfac + sqrt(bfac^2 + 4*afac))/(2*afac); 

                fprintf(fid,'%s\n',sep2);
            end

            % Compute the minimum value of the safety factor for the plies
            SF = min(Sfactor);
            fprintf(fid,'\n%s\nTsai-Wu Safety Factor SF = % 4.3e\n%s\n',...
                sep,SF,sep);
        end
    end
    % Close output file
    fclose(fid);
end



