classdef FTIRexperiment
    properties
        data double = 0
        freqAxis (:,1) double = 0
        volume (:,1) double = 1 % must be in microliters
        pathLength (:,1) double = 12 % must be in micrometers
        radius (:,1) double % not part of constructor method
        gasFactor (:,1) double = 0
        timeInterval (:,1) double = 0 % must be in seconds
        sample string = ""
        diffusionCoeff double = 0
        finalSpectrum double = 0 % to get the concentration at saturation
        finalConc double = 0 % concentration at saturation
    end
    methods
        function f = FTIRexperiment(data,freqAxis,volume,pathLength,timeInterval,sample)
            if nargin > 0
                f.data = data;
                f.freqAxis = freqAxis;
                f.volume = volume;
                f.pathLength = pathLength;
                f.sample = sample;
                f.timeInterval = timeInterval;
                f.radius = sqrt(f.volume*1e9/(pi*f.pathLength)); % returns radius in micrometers
            end
        end
        function plts = plotSpectra(f,specNum,gasFactor)
            % plots the spectra for data f for specNum indices, w optional
            % gasFactor
            %syntax: plotSpectra(f,specNum,gasFactor)
            if numel(specNum) == 1
                if nargin == 3
                    %find the indicies for the amount of spectra desired
                    spectraIndicies = zeros(1,specNum);
                    interval = floor(size(f.data,2)/specNum);
                    for ii = 0:specNum-1
                        spectraIndicies(ii+1) = 1+(ii*interval);
                    end
                    if spectraIndicies(end) ~= size(f.data,2)
                        spectraIndicies = [spectraIndicies size(f.data,2)];
                    end
                    %baseline correction
                    baseline = f.data - f.data(28375,:);
                    %gas line correction
                    %cd("Users/matthewliberatore/Library/CloudStorage/OneDrive-UniversityofPittsburgh/data/ftir_data/Matt/Gas Lines Ref")
                    wha = load("CO2_gas_lines.mat",'gasLines');
                    fixed = baseline-gasFactor.*wha.gasLines;
                elseif nargin == 2
                    %find the indicies for the amount of spectra desired
                    spectraIndicies = zeros(1,specNum);
                    interval = floor(size(f.data,2)/specNum);
                    for ii = 0:specNum-1
                        spectraIndicies(ii+1) = 1+(ii*interval);
                    end
                    if spectraIndicies(end) ~= size(f.data,2)
                        spectraIndicies = [spectraIndicies size(f.data,2)];
                    end
                    %baseline correction
                    baseline = f.data - f.data(28375,:);
                    %gas line correction
                    %cd("Users/matthewliberatore/Library/CloudStorage/OneDrive-UniversityofPittsburgh/data/ftir_data/Matt/Gas Lines Ref")
                    wha = load("CO2_gas_lines.mat",'gasLines');
                    fixed = baseline-f.gasFactor.*wha.gasLines;
                end
            elseif numel(specNum) > 1
                if nargin == 3
                    spectraIndicies = specNum;
                    %baseline correction
                    baseline = f.data - f.data(28375,:);
                    %gas line correction
                    %cd("Users/matthewliberatore/Library/CloudStorage/OneDrive-UniversityofPittsburgh/data/ftir_data/Matt/Gas Lines Ref")
                    wha = load("CO2_gas_lines.mat",'gasLines');
                    fixed = baseline-gasFactor.*wha.gasLines;
                elseif nargin == 2
                    %find the indicies for the amount of spectra desired
                    spectraIndicies = specNum;
                    %baseline correction
                    baseline = f.data - f.data(28375,:);
                    %gas line correction
                    %cd("Users/matthewliberatore/Library/CloudStorage/OneDrive-UniversityofPittsburgh/data/ftir_data/Matt/Gas Lines Ref")
                    wha = load("CO2_gas_lines.mat",'gasLines');
                    fixed = baseline-f.gasFactor.*wha.gasLines;
                end
            end
            plts = plot(f.freqAxis(:,1),fixed(:,spectraIndicies));

            hold on
            xlim([2290 2390])
            xlabel("Frequency (cm^{-1})")
            ylabel("Absorbance (AU)")
            legend(string(spectraIndicies))
        end
        function conc = concOverTime(f)
                %baseline correction
                baseline = f.data - f.data(28375,:);
                %gas line correction
                %cd("Users/matthewliberatore/Library/CloudStorage/OneDrive-UniversityofPittsburgh/data/ftir_data/Matt/Gas Lines Ref")
                wha = load("CO2_gas_lines.mat",'gasLines');
                fixed = baseline-f.gasFactor.*wha.gasLines;
                
                CO2band = fixed(f.freqAxis(:,1) > 2290 & f.freqAxis(:,1) < 2390,:);
                conc = max(CO2band)./(1000*f.pathLength*1e-4);
                conc = conc-conc(1); % DOES THIS NEED TO BE HERE?
        end
        function axis = timeAxis(f)
            % generates time axis for data in f depending on time interval
            % specified
            %syntax: timeAxis(f)
            axis = (0:(size(f.data,2)-1)).*10;
        end
        function plts = plotConcOverTime(f)
            %converts data in f to vector of concentration values of CO2
            %for each spectrum, generates time axis, returns a plot of
            %concentration vs time
            %syntax: plotConcOverTime(f)
            plts = plot(timeAxis(f),concOverTime(f));
            hold on
            xlabel("Time (s)")
            ylabel("Concentration (M)")
            hold off
        end
        function model = uptakeModel(f,a0,k1,k2)
            %y = a0*(1-k1/(k2-k1)*(exp(-k1*x)-exp(-k2*x))-exp(-k1*x));
            model = a0*(1-k1/(k2-k1)*(exp(-k1*timeAxis(f))-exp(-k2*timeAxis(f)))-exp(-k1*timeAxis(f)));
        end
        function u = getDiffusionEquation(f,nmax)
            %returns a symbolic function for the complex diffusion model
            %approximated to nmax terms. The resulting function depends on
            %r,t, and D.
            %syntax: getDiffusionEquation(f,nmax)
            syms a x u b C
            
            A = f.radius; % radius of disk
            %D = 1; % diffusion coeff?
            C = f.finalConc; % initial conc on outside
            
            j0 = besselzero(0,nmax,1); % get nmax zeros of 0th order bessel function of first kind
            
            % constants defined by boundary condition
            c0 = zeros(1,nmax);
            for ii = 1:nmax
                c0(ii) = (C*2)/(j0(ii)*besselj(1,j0(ii)));
            end
            
            u(a,b,x)=0;
            for ii = 1:nmax
                u = u + c0(ii)*besselj(0,j0(ii)/A*a)*exp(-(j0(ii)/A)^2*b*x);
            end
            u = -u+C;
        end
        function f = getFinalConc(f)
            %baseline correction
            baseline = f.finalSpectrum - f.finalSpectrum(28375);
            %gas line correction
            %cd("Users/matthewliberatore/Library/CloudStorage/OneDrive-UniversityofPittsburgh/data/ftir_data/Matt/Gas Lines Ref")
            wha = load("CO2_gas_lines.mat",'gasLines');
            fixed = baseline-f.gasFactor.*wha.gasLines;
            
            CO2band = fixed(f.freqAxis(:,1) > 2290 & f.freqAxis(:,1) < 2390,:);
            f.finalConc = max(CO2band)./(1000*f.pathLength*1e-4);
        end
    end
end