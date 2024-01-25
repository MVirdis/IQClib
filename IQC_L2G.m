classdef IQC_L2G < IQCDescriptor
    %IQC_L2G IQC for L2 gain bound
    %   
    
    properties
        gamma
    end
    
    methods
        function obj = IQC_L2G(gamma, ne, nw)
            %IQC_L2G(gamma, ne, nw)
            
            % Unpack args
            if nargin < 2
                error('IQC_L2G constructor takes 3 params (gamma, ne, nw)');
            end
            mustBePositive(gamma);
            mustBePositive(ne);
            mustBePositive(nw);

            % Construct and assign params
            obj = obj@IQCDescriptor();
            obj.gamma = gamma;
            obj.nv = nw;
            obj.nz = ne;
            obj.nu = 0; nu = 0;
            obj.rho = 0; rho = 0;

            % Assign decision variables
            obj.P11 = 1/gamma*eye(ne);
            obj.P12 = zeros(ne,nw);
            obj.P21 = obj.P12';
            obj.P22 = -gamma*eye(nw);

            % Instantiate Psi
            [psi,~,~,~,~] = basisTF(nu,rho,ne);
            obj.psi11 = psi;
            [psi,~,~,~,~] = basisTF(nu,rho,nw);
            obj.psi22 = psi;

            % Construct additional constraints
            obj.constr = [];
        end
    end
end

