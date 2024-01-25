classdef Delta
    %DELTA IQC description of uncertain operator Delta
    %   
    
    properties
        % Size of the described operator
        nv % Delta output dimension
        nz % Delta input dimension
        num_operators

        % P matrix of the IQC multiplier
        P11 % 11 block
        P12 % 12 block
        P21 % 21 block
        P22 % 22 block
        
        % Dynamics of the multiplier
        psi11
        psi22

        % SDP constraints
        constr

        IQCDescriptors
    end
    
    methods
        function obj = Delta(nz, nv, num_operators)
            mustBePositive(nz);
            mustBePositive(nv);
            mustBePositive(num_operators);
            obj.nv = nv;
            obj.nz = nz;
            obj.num_operators = num_operators;
            obj.IQCDescriptors = cell(1,num_operators);
        end

        function obj = attachIQC(obj, IQC, opId)
            mustBePositive(opId);
            if isempty(obj.IQCDescriptors{opId})
                obj.IQCDescriptors{opId} = {IQC};
            else
                nmults = length(obj.IQCDescriptors{opId});
                obj.IQCDescriptors{opId}{nmults+1} = IQC;
            end
            
            % Launch consistency check
            obj.consistencyCheck(opId);
        end

        function consistencyCheck(obj, opId)
            mustBePositive(opId);
            if isempty(obj.IQCDescriptors{opId})
                warning('Called consistency check on empty multiplier channel');
                return;
            end
            nmults = length(obj.IQCDescriptors{opId});
            
            IQC1 = obj.IQCDescriptors{opId}{1};
            nv = IQC1.nv;
            nz = IQC1.nz;
            for i=2:nmults
                IQCi = obj.IQCDescriptors{opId}{i};
                if nv ~= IQCi.nv || nz ~= IQCi.nz
                    error('Incosistent IQC dimensions across channel');
                end
            end
        end

        function obj = compile(obj)
            % Check whether entire operator is covered
            for i=1:obj.num_operators
                if isempty(obj.IQCDescriptors{i})
                    error('IQC coverage of operator is not complete. Add missing IQCs.');
                end
            end

            for i=1:obj.num_operators
                nmults = length(obj.IQCDescriptors{i});
                if nmults > 1
                    error('Error multiple IQC for single channel not implemented');
                end

                % Select current IQC
                IQCi = obj.IQCDescriptors{i}{1};
                obj.constr = [obj.constr, IQCi.constr];

                if isempty(obj.psi11)
                    obj.psi11 = IQCi.psi11;
                else
                    obj.psi11 = blkdiagtf(obj.psi11, IQCi.psi11);
                end

                if isempty(obj.psi22)
                    obj.psi22 = IQCi.psi22;
                else
                    obj.psi22 = blkdiagtf(obj.psi22, IQCi.psi22);
                end

                if isempty(obj.P11)
                    obj.P11 = IQCi.P11;
                else
                    obj.P11 = blkdiag(obj.P11, IQCi.P11);
                end

                if isempty(obj.P12)
                    obj.P12 = IQCi.P12;
                else
                    obj.P12 = blkdiag(obj.P12, IQCi.P12);
                end

                if isempty(obj.P21)
                    obj.P21 = IQCi.P21;
                else
                    obj.P21 = blkdiag(obj.P21, IQCi.P21);
                end

                if isempty(obj.P22)
                    obj.P22 = IQCi.P22;
                else
                    obj.P22 = blkdiag(obj.P22, IQCi.P22);
                end
            end
        end
    end
end

