classdef IQCDescriptor
    %IQCDESCRIPTOR A single IQC
    %   
    
    properties
        % Size of the described operator
        nv % Delta output dimension
        nz % Delta input dimension

        % P matrix of the IQC multiplier
        P11 % 11 block
        P12 % 12 block
        P21 % 21 block
        P22 % 22 block
        
        % Dynamics of the multiplier
        nu      % McMillan degree
        rho     % pole location
        psi11   % tf object
        psi22   % tf object

        % SDP constraints
        constr
    end
    
    methods
        function obj = IQCDescriptor()
        end
    end
end

