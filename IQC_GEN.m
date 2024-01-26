classdef IQC_GEN < IQCDescriptor
    %IQC_GEN Generic positive-negative multiplier
    %   

    methods
        function obj = IQC_GEN(varargin)
            %IQC_GEN(nz, nv) or IQC_GEN(nz, nv, Ts)
            
            % Unpack args
            if nargin < 2
                error('IQC_GEN not enough inputs');
            end
            nz = varargin{1};
            nv = varargin{2};
            mustBePositive(nz);
            mustBePositive(nv);
            if nargin >= 3
                Ts = varargin{3};
                mustBeNonnegative(Ts);
            else
                Ts = 0;
            end

            % Construct and assign params
            obj = obj@IQCDescriptor();
            obj.nv = nv;
            obj.nz = nz;
            obj.nu = 0; nu = 0;
            obj.rho = 0; rho = 0;
            obj.Ts = Ts;

            % Assign decision variables
            obj.P11 = sdpvar(nz,nz,'symmetric','real');
            obj.P12 = sdpvar(nz,nv,'full','real');
            obj.P21 = obj.P12';
            obj.P22 = sdpvar(nv,nv,'symmetric','real');

            % Instantiate Psi
            [psi,~,~,~,~] = basisTF(nu,rho,nz,Ts);
            obj.psi11 = psi;
            [psi,~,~,~,~] = basisTF(nu,rho,nv,Ts);
            obj.psi22 = psi;

            % Construct additional constraints
            obj.constr = [(obj.P11 >= 0), (obj.P22 <= 0)];
        end
    end
end

