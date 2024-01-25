classdef IQCProblem
    %IQCPROBLEM Analysis problem
    %   
    
    properties
        % Partitioned Linear System
        Gvz
        Gve
        Gwz
        Gwe

        % Dimensions
        nv
        nz
        nw
        ne

        % IQC Descriptions
        Delta
        DeltaP

        % Whether performance is included
        perf
        compiled

        constr
    end
    
    methods
        function obj = IQCProblem(varargin)
            %IQCProblem(Gvz, Delta) or IQCProblem(Gvz, Delta, Gve, Gwz, Gwe, DeltaP)

            obj.perf = false;
            obj.compiled = false;

            % Unpack input arguments
            % Gvz, Delta
            % Gvz, Delta, Gve, Gwz, Gwe, DeltaP
            if nargin == 2
                Gvz = varargin{1};
                obj.Delta = varargin{2};
            elseif nargin == 6
                Gvz = varargin{1};
                obj.Delta = varargin{2};
                Gve = varargin{3};
                Gwz = varargin{4};
                Gwe = varargin{5};
                obj.DeltaP = varargin{6};
                obj.perf = true;
            elseif nargin < 2 || (nargin > 2 && nargin < 6)
                error('Not enough input arguments');
            end

            % Make sure everything is in ss
            Gvz_ = ss(Gvz);
            if obj.perf
                Gve_ = ss(Gve);
                Gwz_ = ss(Gwz);
                Gwe_ = ss(Gwe);
            end

            % Infer dimensions
            [nz,nv] = size(Gvz_);
            if obj.perf
                [ne,nw] = size(Gwe_);
            end

            % Check whether the system is strictly stable
            if obj.perf
                G = [Gvz_,Gwz_;
                     Gve_,Gwe_];
            else
                G = Gvz_;
            end
            eA = eig(G.A);
            if ~all(real(eA) < 0)
                error('G must be strictly stable for IQC analysis');
            end

            % Assign
            obj.Gvz = Gvz_;
            if obj.perf
                obj.Gve = Gve_;
                obj.Gwe = Gwe_;
                obj.Gwz = Gwz_;
            end
            obj.nv = nv;
            obj.nz = nz;
            if obj.perf
                obj.ne = ne;
                obj.nw = nw;
            end

        end

        function obj = compile(obj)
            if ~obj.perf
                obj = obj.compileStab();
            else
                obj = obj.compilePerf();
%                 error('Not yet implemented');
            end
        end

        function obj = compileStab(obj)
            obj.constr = obj.Delta.constr;

            % Make ss of Psi and [Gvz;I]
            Gvz_ = [obj.Gvz; ss(eye(size(obj.Gvz,2)))];
            psi = blkdiagtf(obj.Delta.psi11,obj.Delta.psi22);
            G = series(Gvz_, ss(psi));

            nA = length(G.A);
            M = vstack_IO(G.A,G.B,G.C,G.D);
            X = sdpvar(nA,nA,'symmetric','real');
            T = blkdiag([zeros(nA),X;X,zeros(nA)], ...
                        [obj.Delta.P11,obj.Delta.P12; ...
                         obj.Delta.P21,obj.Delta.P22]);
            
            obj.constr = [obj.constr, (M'*T*M <= 0)];

            obj.compiled = true;
        end

        function obj = compilePerf(obj)
            obj.constr = [obj.Delta.constr, obj.DeltaP.constr];

            % Make ss of Psi and [Gvz;I]
            G_ = [obj.Gvz, obj.Gwz;
                  ss(eye(obj.nv)), ss(zeros(obj.nv, obj.nw))];
            G_ = [G_; obj.Gve, obj.Gwe; ss(zeros(obj.nw, obj.nv)), ss(eye(obj.nw))];
            psi = blkdiagtf(obj.Delta.psi11,obj.Delta.psi22, ...
                            obj.DeltaP.psi11,obj.DeltaP.psi22);
            G = series(G_, ss(psi));

            nA = length(G.A);
            M = vstack_IO(G.A,G.B,G.C,G.D);
            X = sdpvar(nA,nA,'symmetric','real');
            T = blkdiag([zeros(nA),X;X,zeros(nA)], ...
                        [obj.Delta.P11,obj.Delta.P12; ...
                         obj.Delta.P21,obj.Delta.P22],...
                        [obj.DeltaP.P11,obj.DeltaP.P12; ...
                         obj.DeltaP.P21,obj.DeltaP.P22]);
            
            obj.constr = [obj.constr, (M'*T*M <= 0)];

            obj.compiled = true;
        end

        function solve(obj)
            if ~obj.compiled
                error('Before solving it the problem needs to be compiled.');
            end

            options = sdpsettings('verbose',1,'solver','mosek');
            optimize(obj.constr,[],options);
        end
    end
end

