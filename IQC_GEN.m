classdef IQC_GEN < IQCDescriptor
    %IQC_GEN Generic positive-negative multiplier
    %   

    methods
        function obj = IQC_GEN(nz, nv)
            %IQC_GEN Construct an instance of this class
            
            mustBePositive(nz);
            mustBePositive(nv);

            % Construct and assign params
            obj = obj@IQCDescriptor();
            obj.nv = nv;
            obj.nz = nz;
            obj.nu = 0; nu = 0;
            obj.rho = 0; rho = 0;

            % Assign decision variables
            obj.P11 = sdpvar(nz,nz,'symmetric','real');
            obj.P12 = sdpvar(nz,nv,'full','real');
            obj.P21 = obj.P12';
            obj.P22 = sdpvar(nv,nv,'symmetric','real');

            % Instantiate Psi
            [psi,~,~,~,~] = basisTF(nu,rho,nz);
            obj.psi11 = psi;
            [psi,~,~,~,~] = basisTF(nu,rho,nv);
            obj.psi22 = psi;

            % Construct additional constraints
            obj.constr = [(obj.P11 >= 0), (obj.P22 <= 0)];
        end
    end
end

