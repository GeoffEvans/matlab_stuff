function [X,Y1,Y2,Z] = Test()

A = 1:100;
B = 1:300;

X = mapping();
Y1 = forloopMeshCell();
Y1 = forloopCell();
Y2 = forloopArray();
Z = multidim();

    function f = fun(a,b)
        f = [cos(b) .* sin(a), sin(b) .* sin(a), cos(a)];
    end
    function X = mapping()
        [aMesh,bMesh] = meshgrid(A,B);
        X = arrayfun(@fun, aMesh, bMesh,'uniformoutput', false);
    end
    function Y2 = forloopArray()
        Y2 = zeros(3,length(B),length(A));
        for k = B
            for l = A
                Y2(:,k,l) = fun(l,k);
            end
        end
    end
    function Y1 = forloopCell()
        Y1 = cell(length(B),length(A));
        for k = B
            for l = A
                Y1{k,l} = fun(l,k);
            end
        end
    end
    function Z = multidim()
        Z1 = cos(B)' * sin(A);
        Z2 = sin(B)' * sin(A);
        Z3 = ones(length(B),1) * cos(A);
        Z = zeros(3,length(B),length(A));
        Z(1,:,:) = Z1;
        Z(2,:,:) = Z2;
        Z(3,:,:) = Z3;
    end
end