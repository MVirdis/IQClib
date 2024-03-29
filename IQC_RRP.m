classdef IQC_RRP < IQCDescriptor
    %IQC_RRP Real Repeated Parameter
    %   
    
    properties
        alpha % bound such that |delta| <= alpha
    end
    
    methods
        function obj = IQC_RRP(varargin)
            %IQC_RRP(alpha, nq) or IQC_RRP(alpha, nq, Ts) or IQC_RRP(alpha, nq, nu, rho) or IQC_RRP(alpha, nq, nu, rho, Ts)

            % Unpack args
            if nargin < 2
                error('IQC_RRP constructor takes at least 2 params (alpha, size)');
            end
            
            alpha = varargin{1};
            nq = varargin{2};
            mustBePositive(alpha);
            mustBePositive(nq);

            Ts = 0;
            nu = 0;
            rho = 0;
            if nargin == 3
                Ts = varargin{3};
            end
            if nargin >= 4
                nu = varargin{3};
                rho = varargin{4};
            end
            if nargin >= 5
                Ts = varargin{5};
            end
            mustBeNonnegative(nu);
            mustBeNegative(rho);
            mustBeNonnegative(Ts);

            % Construct and assign params
            obj = obj@IQCDescriptor();
            obj.alpha = alpha;
            obj.nv = nq;
            obj.nz = nq;
            obj.nu = nu;
            obj.rho = rho;
            obj.Ts = Ts;

            % Assign decision variables
            obj.P11 = sdpvar(nq*(nu+1),nq*(nu+1),'symmetric','real');
            obj.P12 = sdpvar(nq*(nu+1),nq*(nu+1),'full','real');
            obj.P21 = obj.P12';
            obj.P22 = -obj.P11;

            % Instantiate Psi
            [psi,Anu,Bnu,Cnu,Dnu] = basisTF(nu,rho,nq,Ts);
            obj.psi11 = alpha*psi;
            obj.psi22 = psi;

            % Construct additional constraints
            M = vstack_IO(Anu, Bnu, Cnu, Dnu);
            nAnu = length(Anu);
            Xnu = sdpvar(nAnu, nAnu, 'symmetric', 'real');
            if Ts > 0
                T = [-Xnu, zeros(nAnu); zeros(nAnu), Xnu];
            else
                T = [zeros(nAnu), Xnu; Xnu, zeros(nAnu)];
            end
            T = blkdiag(T,obj.P11);
            obj.constr = [(obj.P12+obj.P12' == 0), (M'*T*M >= 0)];
        end
    end
end

