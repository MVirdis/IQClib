classdef IQC_L2G < IQCDescriptor
    %IQC_L2G IQC for L2 gain bound
    %   
    
    properties
        gamma
    end
    
    methods
        function obj = IQC_L2G(varargin)
            %IQC_L2G(gamma, ne, nw) or IQC_L2G(gamma, ne, nw, Ts)
            
            % Unpack args
            if nargin < 3
                error('IQC_L2G not enough inputs');
            end
            gamma = varargin{1};
            ne = varargin{2};
            nw = varargin{3};
            mustBePositive(gamma);
            mustBePositive(ne);
            mustBePositive(nw);
            if nargin >= 4
                Ts = varargin{4};
                mustBeNonnegative(Ts);
            else
                Ts = 0;
            end

            % Construct and assign params
            obj = obj@IQCDescriptor();
            obj.gamma = gamma;
            obj.nv = nw;
            obj.nz = ne;
            obj.nu = 0; nu = 0;
            obj.rho = 0; rho = 0;
            obj.Ts = Ts;

            % Assign decision variables
            obj.P11 = 1/gamma*eye(ne);
            obj.P12 = zeros(ne,nw);
            obj.P21 = obj.P12';
            obj.P22 = -gamma*eye(nw);

            % Instantiate Psi
            [psi,~,~,~,~] = basisTF(nu,rho,ne,Ts);
            obj.psi11 = psi;
            [psi,~,~,~,~] = basisTF(nu,rho,nw,Ts);
            obj.psi22 = psi;

            % Construct additional constraints
            obj.constr = [];
        end
    end
end

