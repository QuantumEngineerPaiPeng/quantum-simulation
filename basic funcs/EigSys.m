% 20180620
% eigensystem class
% property:
% D: eigenvalues
% V: eigenvectors
% created by Pai Peng
classdef EigSys% < matlab.mixin.Copyable
    properties
        D
        V
    end
    
    methods
        function obj=EigSys(V,D)
            % test
            if nargin==2
                if size(D,2)==1 || isempty(D)
                    obj.D=D;
                else
                    if size(D,1)==1
                        obj.D=D.';
                    else
                        error('D must be a vector')
                    end
                end
                obj.V=V;
            else
                error('Number of arguments must be 2')
            end
        end
    end
end