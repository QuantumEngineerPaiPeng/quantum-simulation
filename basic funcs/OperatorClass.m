% 20180619
% class for operators
% created by Pai Peng
classdef OperatorClass < matlab.mixin.Copyable
    properties
        matrix % Matrix cell. Element in the cell represents the matrix in a symmetry sector.
        sym % Symmtry type, eg. 'kpz'
        L % System size
        eigsys % EigSys class cell. Element in the cell contains the eigenvalues and eigenvectors. See also EigSys.
    end
    
    methods
        function obj=OperatorClass(L, model, pre, bc, cplist, direction)
            % create an operator, with argument L, model, prefactor, boundary codition, coupling strength list, direction
            % nargin can be 0, 1, 3, 5 or 6.
            switch nargin
                case 0
                    
                case 1
                    obj.L=L;
                case 3
                    obj.L=L;
                    temp=model-'x'+1;
                    if temp~=1 && temp~=2 &&temp~=3
                        error('When only 2 arguments present, 1st is L, 2nd is x, y, or z')
                    end
                    Sigma_i=IndivPauliSparse(L);
                    matr=Sigma_i{temp}{1};
                    for p=2:L
                        matr=matr+Sigma_i{temp}{p};
                    end
                    obj.matrix={pre*matr};
                case 6
                    obj.matrix={pre*Hamiltonian(L,bc,cplist,model,direction)};
                    obj.L=L;
                case 5
                    obj.matrix={pre*Hamiltonian(L,bc,cplist,model)};
                    obj.L=L;
                otherwise
                    error('invalid number of inputs')
            end
        end
        
        function symmetrize(self,symcode)
            % Divide the matrix in symmetry sectors. Type of symmetry
            % specified by the argument, eg. 'kpz' or 'kp'.
            if isempty(symcode)
                return
            end
            if ~isempty(self.sym)
                warning([inputname(1),' is already symmetrized'])
            end
            if or(strcmp(symcode,'kpz'),strcmp(symcode,'kp'))
                if isunix
                    filename=sprintf('~/Dropbox (MIT)/grad/research/codes/basic funcs/AllkpzSym_%d.mat',self.L);
                else
                    filename=sprintf('C:\\Users\\Pai\\Dropbox (MIT)\\grad\\research\\codes\\basic funcs\\AllkpzSym_%d.mat',self.L);
                end
                if exist(filename,'file')
                    load(filename)
                else
                    pt=AllkpzProject(self.L);
                    save(filename,'pt')
                end
                switch symcode
                    case 'kpz'
                        self.sym=symcode;
                        temp=pt;
                    case 'kp'
                        self.sym=symcode;
                        temp={};
                        for p=1:length(pt)/2
                            temp{end+1}=[pt{p},pt{length(pt)/2+p}];
                        end
                end
                matr={};
                for p=1:length(temp)
                    if isempty(temp{p})
                        matr{end+1}=[];
                    else
                        matr{end+1}=temp{p}'*self.matrix{1}*temp{p};
                    end
                end
                self.matrix=matr;
                self.eigsys=[];
            else
                if strcmp(symcode,'Mz')
                    Mz=MzProject(self.L);
                    matr=cell(self.L+1,1);
                    for p=1:(self.L+1)
                        matr{p}=self.matrix{1}(Mz{p},Mz{p});
                    end
                    self.sym=symcode;
                    self.matrix=matr;
                    self.eigsys=[];
                else
                    error('Invalid symmetry type')
                end
            end
        end
        
        function r=plus(o1,o2)
            % Add the matrix of two operators, the two must have the same
            % system size and symmetry type.
            if IsSameSym({o1,o2})
                matr=cell(1,length(o1.matrix));
                for p=1:length(o1.matrix)
                    matr{p}=o1.matrix{p}+o2.matrix{p};
                end
                r=copy(o1);
                r.eigsys=[];
                r.matrix=matr;
            else
                error('Adding Hamiltonian with different L or symmetry')
            end
        end
        
        function r=minus(o1,o2)
            % Subtract the matrix of two operators, the two must have the same
            % system size and symmetry type.
            if IsSameSym({o1,o2})
                matr=cell(1,length(o1.matrix));
                for p=1:length(o1.matrix)
                    matr{p}=o1.matrix{p}-o2.matrix{p};
                end
                r=copy(o1);
                r.eigsys=[];
                r.matrix=matr;
            else
                error('Subtracting Hamiltonian with different L or symmetry')
            end
        end
        
        function r=mpower(o1,b)
            % Power of the matrix of an operator
            if isempty(o1.eigsys)
                o1.diagonalize
            end
            if isscalar(b)
                r=copy(o1);
                for p=1:length(o1.matrix)
                    r.eigsys{p}.D=(o1.eigsys{p}.D).^b;
                    r.matrix{p}=r.eigsys{p}.V*sparse(diag(r.eigsys{p}.D))*r.eigsys{p}.V';
                end
                
            else
                error('The power must be a scalar')
            end
        end
        
        function r=ctranspose(o1)
            % Hermittian conjugate of the matrix of an operator
            r=copy(o1);
            for p=1:length(o1.matrix)
                if ~isempty(r.eigsys)
                    r.eigsys{p}.D=conj(o1.eigsys{p}.D);
                end
                r.matrix{p}=o1.matrix{p}';
            end
            
        end
        
        function r=mtimes(o1,o2)
            % Multiply the matrices of two operators or a scalar and matrices of two operators
            tag=0; % 0 for invalid input
            if isnumeric(o1) && isscalar(o1) && isa(o2,'OperatorClass')
                r=copy(o2);
                s=o1;
                tag=1; % 1 for a scalar times an operator
            else
                if isnumeric(o2) && isscalar(o2) && isa(o1,'OperatorClass')
                    r=copy(o1);
                    s=o2;
                    tag=1;
                end
            end
            
            
            
            if isa(o1,'OperatorClass') && isa(o2,'OperatorClass')
                if IsSameSym({o1,o2})
                    tag=2; % 2 for multiplying two operators
                else
                    error('Multiplying Hamiltonian with different L or symmetry')
                end
            end
            
            switch tag
                case 1
                    for p=1:length(r.matrix)
                        r.matrix{p}=r.matrix{p}*s;
                        if ~isempty(r.eigsys)
                            r.eigsys{p}.D=r.eigsys{p}.D*s;
                        end
                    end
                case 2
                    matr=cell(1,length(o1.matrix));
                    for p=1:length(o1.matrix)
                        matr{p}=o1.matrix{p}*o2.matrix{p};
                    end
                    r=copy(o1);
                    r.eigsys=[];
                    r.matrix=matr;
                otherwise
                    error('Cannot multiply. Argument must be (scalar, OperatorClass) or (OperatorClass, OperatorClass)')
            end
            
        end
        
        function diagonalize(self)
            % diagonalize the all matrices in the matrix property.
            % Eigenvalues and eigenvectors stored in the eigsys property.
            % The argument is optional, indicating whether the operator is
            % Hermittian. Default value is True (Hermittian).
            if ~isempty(self.eigsys)
                warning([inputname(1),' is already diagonalized.'])
            end
            eigcell=cell(1,length(self.matrix));
            for p=1:length(self.matrix)
                if max(max(abs(self.matrix{p}-self.matrix{p}')))>1e-13
                    % non-Hermitian ase
                    [V,D]=schur(full(self.matrix{p}));
                    eigcell{p}=EigSys(V,diag(D));
                else % Hermitian case
                    [V,D]=schur(full((self.matrix{p}+self.matrix{p}')/2));
                    eigcell{p}=EigSys(V,diag(D));
                end
                
            end
            self.eigsys=eigcell;
        end
        
        function r=trace(self)
            r=0;
            for p=1:length(self.matrix)
                r=r+trace(self.matrix{p});
            end
        end
        
        function normalize(self)
            t=self.trace;
            if abs(t)<1e-13
                error('Normalizing a operator with zero trace')
            else
                if t<0
                    warning('Normalizing a operator with negative trace')
                else
                    no=self*(1/t);
                    self.matrix=no.matrix;
                    self.eigsys=no.eigsys;
                end
            end
        end
    end
end