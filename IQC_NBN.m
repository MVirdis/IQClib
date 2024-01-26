classdef IQC_NBN < IQCDescriptor
    %IQC_NBN IQC for norm-bounded nonlinearities
    %   
    
    properties
        alpha
    end
    
    methods
        function obj = IQC_NBN(varargin)
            %IQC_NBN(alpha, nz, nv) or IQC_NBN(alpha, nz, nv, Ts)

            % Unpack args
            if nargin < 3
                error('IQC_NBN not enough inputs');
            end
            
            alpha = varargin{1};
            nz = varargin{2};
            nv = varargin{3};
            mustBePositive(alpha);
            mustBePositive(nz);
            mustBePositive(nv);
            if nargin >= 4
                Ts = varargin{4};
                mustBeNonnegative(Ts);
            else
                Ts = 0;
            end

            % Construct and assign params
            obj = obj@IQCDescriptor();
            obj.alpha = alpha;
            obj.nv = nv;
            obj.nz = nz;
            obj.nu = 0; nu = 0;
            obj.rho = 0; rho = 0;
            obj.Ts = Ts;

            % Assign decision variables
            obj.P11 = eye(nz);
            obj.P12 = zeros(nz,nv);
            obj.P21 = zeros(nv,nz);
            obj.P22 = -eye(nv);

            % Instantiate Psi
            [psi,~,~,~,~] = basisTF(nu,rho,nz,Ts);
            obj.psi11 = alpha*psi;
            [psi,~,~,~,~] = basisTF(nu,rho,nv,Ts);
            obj.psi22 = psi;

            % Construct additional constraints
            obj.constr = [];
        end
    end
end

