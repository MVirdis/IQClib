classdef IQC_NBN < IQCDescriptor
    %IQC_NBN IQC for norm-bounded nonlinearities
    %   
    
    properties
        alpha
    end
    
    methods
        function obj = IQC_NBN(alpha, nz, nv)
            %IQC_NBN(alpha, nz, nv)
            
            mustBePositive(alpha);
            mustBePositive(nz);
            mustBePositive(nv);

            % Construct and assign params
            obj = obj@IQCDescriptor();
            obj.alpha = alpha;
            obj.nv = nv;
            obj.nz = nz;
            obj.nu = 0; nu = 0;
            obj.rho = 0; rho = 0;

            % Assign decision variables
            obj.P11 = eye(nz);
            obj.P12 = zeros(nz,nv);
            obj.P21 = zeros(nv,nz);
            obj.P22 = -eye(nv);

            % Instantiate Psi
            [psi,~,~,~,~] = basisTF(nu,rho,nz);
            obj.psi11 = alpha*psi;
            [psi,~,~,~,~] = basisTF(nu,rho,nv);
            obj.psi22 = psi;

            % Construct additional constraints
            obj.constr = [];
        end
    end
end

