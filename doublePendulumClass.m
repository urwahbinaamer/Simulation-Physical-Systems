% To run, enter in command window: 
%    pendulum_sim = doublePendulumClass;

% Ref: "Physically Based Modeling: Principles and Practice
% Constrained Dynamics" http://www.cs.cmu.edu/~baraff/sigcourse/notesf.pdf
% accessed Apr 2024

classdef doublePendulumClass
    
    properties
        numberOfBodies (1, 1) double {mustBeInteger, mustBePositive} = 4
        numberOfConstraints (1, 1) double {mustBeInteger, mustBePositive} = 6
        
        % frequency of physics calculation (per second)
        computeFrequency = 5000
        
        % update frequency of plot
        framerate = 50
        
        % duration of simulation in seconds
        simDuration = 10;
        
        % defines masses of particles
        massesArray = [0.1 0.1 0.1 0.1]
        
        % defines intial orientation of links, positive is clockwise from
        % east
        thetaArray = [-30 -0]
        
        % defines length of links
        lengthArray = [0.3 0.3]
        
        % gravity
        g = -9.81;
        
        Ks = 50; % spring constant 
        Kd = 10; % damping constant
        
        moveFrequency = 0.9
        
        numbStateParams
        
        F1 % standard function 1
        timestep
        frameTimestep
        J
        Jdot
        W
        Q
        q
        qDot
        C
        t = 0
        %animLine1
        %animLine2
        axisLimits
    end
    
    methods
        
        function obj = run(obj)
            obj = obj.computeConstants;
            obj = obj.calculateJandF1;
            % pendulum_sim.createFigure;
            obj = obj.solver;
        end
        
        function [obj] = doublePendulumClass
            % createObject; constructs an instance of this class, and runs
            % code
            
%             obj.numberOfConstraints = numberOfConstraints;
%             obj.numberOfBodies = numberOfBodies;
%             obj.computeFrequency = computeFrequency;
%             obj.simDuration = simDuration;
%             obj.massesArray = massesArray;
%             obj.thetaArray = thetaArray;
%             obj.lengthArray = lengthArray;
%             obj.axisLimits = axisLimits;
            obj = computeConstants(obj);
            obj = calculateJandF1(obj);
            obj = solver(obj);
            
        end
        
        function obj = computeConstants(obj)
            %% computes / initializes constants
            obj.numbStateParams = 2 * obj.numberOfBodies;% is this actually the max numb possible?
            obj.timestep = 1/obj.computeFrequency;
            obj.frameTimestep = 1/obj.framerate;
            obj.Jdot = zeros(obj.numberOfConstraints, obj.numbStateParams);
            obj.qDot = zeros(obj.numbStateParams, 1);
            obj.Q = zeros(obj.numbStateParams, 1);
            obj.C = zeros(obj.numberOfConstraints, 1);
            
            %% Define vector q
            obj.q = [0; 0; obj.lengthArray(1) * cosd(obj.thetaArray(1)); obj.lengthArray(1) * -sind(obj.thetaArray(1))];
            obj.q(5, 1) = obj.q(3, 1);
            obj.q(6, 1) = obj.q(4, 1);
            obj.q(7, 1) = obj.q(5) + obj.lengthArray(2) * cosd(obj.thetaArray(2));
            obj.q(8, 1) = obj.q(6) + obj.lengthArray(2) * -sind(obj.thetaArray(2));
                       
            %% Define W
            obj.W = zeros(obj.numbStateParams);
            obj.W(1, 1) = 1/obj.massesArray(1);
            obj.W(2, 2) = 1/obj.massesArray(1);
            obj.W(3, 3) = 1/obj.massesArray(2);
            obj.W(4, 4) = 1/obj.massesArray(2);
            obj.W(5, 5) = 1/obj.massesArray(3);
            obj.W(6, 6) = 1/obj.massesArray(3);
            obj.W(7, 7) = 1/obj.massesArray(4);
            obj.W(8, 8) = 1/obj.massesArray(4);
            
            %% Define Q
            %obj.Q(1) = 10*sin(2*pi*obj.t*obj.moveFrequency);
            obj.Q(2) = obj.massesArray(1) * obj.g;
            obj.Q(4) = obj.massesArray(2) * obj.g;
            obj.Q(6) = obj.massesArray(3) * obj.g;
            obj.Q(8) = obj.massesArray(4) * obj.g;
        end
        
        function obj = calculateJandF1(obj)
            %% F1 function (denominator)
            obj.F1(1) = sqrt((obj.q(3, 1) - obj.q(1, 1))^2 + (obj.q(4, 1) - obj.q(2, 1))^2); %sqrt((x2-x1)^2 + (y2-y1)^2);
            %obj.F1(2) = sqrt((obj.q(5, 1) - obj.q(3, 1))^2 + (obj.q(6, 1) - obj.q(4, 1))^2);
            obj.F1(2) = sqrt((obj.q(7, 1) - obj.q(5, 1))^2 + (obj.q(8, 1) - obj.q(6, 1))^2); %sqrt((x4-x3)^2 + (y4-y3)^2);
            
            %% Define J
            obj.J = zeros(obj.numberOfConstraints, obj.numbStateParams);
            obj.J(1, 1) = 1;
            obj.J(2, 2) = 1;
            obj.J(3, 1) = (obj.q(1) - obj.q(3))/obj.F1(1);
            obj.J(3, 2) = (obj.q(2) - obj.q(4))/obj.F1(1);
            obj.J(3, 3) = (obj.q(3) - obj.q(1))/obj.F1(1);
            obj.J(3, 4) = (obj.q(4) - obj.q(2))/obj.F1(1);
            
            obj.J(4, 3) = -1;
            obj.J(4, 5) = 1;
            obj.J(5, 4) = -1;
            obj.J(5, 6) = 1;
            
            obj.J(6, 5) = (obj.q(5) - obj.q(7))/obj.F1(2);
            obj.J(6, 6) = (obj.q(6) - obj.q(8))/obj.F1(2);
            obj.J(6, 7) = (obj.q(7) - obj.q(5))/obj.F1(2);
            obj.J(6, 8) = (obj.q(8) - obj.q(6))/obj.F1(2);
            
        end
        
        function obj = updateJ(obj)
            %% Define J
            %obj.J = zeros(obj.numberOfConstraints, obj.numbStateParams);
            %obj.J(1, 1) = 1;
            %obj.J(2, 2) = 1;
            obj.J(3, 1) = (obj.q(1) - obj.q(3))/obj.F1(1);
            obj.J(3, 2) = (obj.q(2) - obj.q(4))/obj.F1(1);
            obj.J(3, 3) = (obj.q(3) - obj.q(1))/obj.F1(1);
            obj.J(3, 4) = (obj.q(4) - obj.q(2))/obj.F1(1);
            
            %obj.J(4, 3) = -1;
            %obj.J(4, 5) = 1;
            %obj.J(5, 4) = -1;
            %obj.J(5, 6) = 1;
            
            obj.J(6, 5) = (obj.q(5) - obj.q(7))/obj.F1(2);
            obj.J(6, 6) = (obj.q(6) - obj.q(8))/obj.F1(2);
            obj.J(6, 7) = (obj.q(7) - obj.q(5))/obj.F1(2);
            obj.J(6, 8) = (obj.q(8) - obj.q(6))/obj.F1(2);
            
        end
        
        function [Cdot, obj] = calculateCandCdot(obj)
            C_old = obj.C;
            
            obj.C(1, 1) = obj.q(1);%-0.1*sin(2*pi*obj.t*obj.moveFrequency);%obj.q(1);
            obj.C(2, 1) = obj.q(2);
            obj.C(3, 1) = sqrt( (obj.q(3) - obj.q(1))^2+(obj.q(4) - obj.q(2))^2)-0.3;
            obj.C(4, 1) = obj.q(5) - obj.q(3);
            obj.C(5, 1) = obj.q(6) - obj.q(4);
            obj.C(6, 1) = sqrt((obj.q(7) - obj.q(5))^2+(obj.q(8) - obj.q(6))^2)-0.3;
            
            Cdot = (obj.C - C_old) / obj.timestep;
            
        end
        
        function obj = solver(obj)
            %
            %% Creates figure
            % creates lines
            figure('Position', [10 50 1500 700]);
%             animLine1 = animatedline([obj.q(1) obj.q(3) ],...
%                 [obj.q(2) obj.q(4)], 'LineWidth', 10, ...
%                 'Marker', 'o', 'MarkerSize', 9, 'MarkerFaceColor', 'w');
            
            animLine1 = animatedline([obj.q(1) obj.q(3) obj.q(5) obj.q(7)],...
                [obj.q(2) obj.q(4) obj.q(6) obj.q(8)], 'LineWidth', 10, ...
                'Marker', 'o', 'MarkerSize', 9, 'MarkerFaceColor', 'w');
            
            
            %animLine2 = animatedline([)], [], 'LineWidth', 10, 'Marker', 'o', 'MarkerSize', 9, 'MarkerFaceColor', 'w');
            axis equal
            xlim([-0.6 0.6])
            ylim([-0.8 0.4])
            pause(0.5)
            
            %% Solver
            
            i = 1;
            totalIterationCount = obj.computeFrequency * obj.simDuration;
            %intialize time to zero
            obj.t = 0;
            %lamda_old = -obj.Jdot * obj.qDot - obj.J * obj.W * obj.Q - (obj.Ks * obj.C) - (obj.Kd * 0);
            
            %loopRateControl = rateControl(obj.computeFrequency);
            %reset(loopRateControl)
            loopTimeReduction = 0.67; % loop takes too long, hence this factor
            totalTime = tic;
            while i <= totalIterationCount
                loopTime = tic;
                i = i+1;
                obj.t = obj.t + obj.timestep;
                %q_old = obj.q;
                [Cdot, obj] = calculateCandCdot(obj);
                
                A = obj.J * obj.W * obj.J';
                b = -obj.Jdot * obj.qDot - obj.J * obj.W * obj.Q - (obj.Ks * obj.C) - (obj.Kd * Cdot);
                
                % solves for lamda
                lamda = A \ b; % gaussian elim time 0.246 sec vs 0.366 sec by conjgrad for 30k calls...
                %lamda = conjgrad(A,b, lamda_old);
                
                %lamda_old = lamda;
                fHat = obj.J'*lamda;
                
                qDoubleDot = obj.W * (obj.Q + fHat); % acceleration ; try the runge-kutta method instead of this?
                obj.qDot = obj.qDot + (qDoubleDot .* obj.timestep);
                obj.q = obj.q + (obj.qDot .* obj.timestep);
                
                J_old = obj.J;
                
                obj = updateJ(obj);
                
                obj.Jdot = (obj.J - J_old) ./ obj.timestep;
                
                if mod(obj.t, obj.frameTimestep) <= (obj.timestep)
                    clearpoints(animLine1)
                    addpoints(animLine1, [obj.q(1) obj.q(3) obj.q(5) obj.q(7)], [obj.q(2) obj.q(4) obj.q(6) obj.q(8)])
                    %[obj.q(1) obj.q(3)], [obj.q(2) obj.q(4)]
                    drawnow
                    potentialEnergy = obj.massesArray(1) * -obj.g * (obj.q(2) + obj.q(4) + obj.q(6) + obj.q(8));
                    kineticEnergy = 0.5 * obj.massesArray(1)*( obj.qDot(1)^2+obj.qDot(2)^2 + obj.qDot(3)^2+obj.qDot(4)^2 + obj.qDot(5)^2+obj.qDot(6)^2 + obj.qDot(7)^2+obj.qDot(8)^2 );
                    
                    titleText = sprintf('Time: %.2f s, Energy: %.3f J', obj.t, (potentialEnergy+kineticEnergy));
                    title(titleText)
                    %ylim([(q(4)-0.1) (length+0.1)])
                    %toc(startTime)
                end
                
                %waitfor(loopRateControl);
                pause((obj.timestep - toc(loopTime))*loopTimeReduction)
                %toc(loopTime)
            end
            
            toc(totalTime)
            
        end
        
    end
    
end
