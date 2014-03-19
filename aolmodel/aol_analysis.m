function [ output_args ] = aolAnalysis( input_args )
%AOLANALYSIS Summary of this function goes here
%   Detailed explanation goes here



        function PlotRays(plotRays,zFocusModel)
            if plotRays
                figure()
                hold on;
                for m=1:numOfAods+1
                    zValAsArray = repmat(zPlanesAod(m),1,4);
                    fill3([max(x(:)) max(x(:)) min(x(:)) min(x(:))],[max(y(:)) min(y(:)) min(y(:)) max(y(:))], zValAsArray, zValAsArray);
                end
                zValAsArray = repmat(zFocusModel,1,4);
                fill3([max(x(:)) max(x(:)) min(x(:)) min(x(:))],[max(y(:)) min(y(:)) min(y(:)) max(y(:))], zValAsArray, zValAsArray);
                alpha(0.1)
                for q = 1:numOfPerturbations
                    indicesForQthPerturbation = (1:numOfRaysPerPerturbation)+(q-1)*numOfRaysPerPerturbation;
                    plot3(x(:,indicesForQthPerturbation),y(:,indicesForQthPerturbation),z(:,indicesForQthPerturbation),'r');
                end
                grid on;
                grid minor;
                axis square;
                xlabel('x')
                ylabel('y')
                zlabel('z')
                hold off;
            end
        end

function [ prodEffDeflection ] = AnalysePerformance(eff, numToOptimise, xFocus, yFocus, theta, phi)
            prodEffSingleRay = prod(eff(1:numToOptimise,:),1);
            prodEffDeflectionMat = reshape(prodEffSingleRay,numOfTimes*numOfPositions,numOfDeflections*numOfPerturbations); 
            prodEffDeflection = mean(prodEffDeflectionMat,1); % average rays for each deflection
            maxFracAngleError = max(fractionalAngleErrorMax);
            
            % REMOVE THIS AFTER DOING SCANNING WORK
            %prodEffDeflection = mean(reshape(prodEffSingleRay,numOfTimes,numOfPositions*numOfDeflections*numOfPerturbations),2); 
        end

                function CompareModelToSimpleAngles(kIn, kOut, localFreq)
                    modelAngle = acos( dot(kOut,kIn) ./ (mag(kOut).*mag(kIn)) );
                    isotropicAngle = wavelengthVac * localFreq / V;
                    fractionalAngleError = abs((isotropicAngle - modelAngle)./modelAngle);
                    fractionalAngleErrorMax = fractionalAngleErrorMax + (fractionalAngleError > fractionalAngleErrorMax) .* (fractionalAngleError - fractionalAngleErrorMax);  
                end

end

