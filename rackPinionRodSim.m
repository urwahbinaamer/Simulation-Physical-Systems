% To run, enter in command window: 
%    pendulum_sim = rackPinionRodSim;
%    pendulum_sim = pendulum_sim.run;

% Ref: "Physically Based Modeling: Principles and Practice
% Constrained Dynamics" http://www.cs.cmu.edu/~baraff/sigcourse/notesf.pdf
% accessed Apr 2024

classdef rackPinionRodSim
    
    properties
        numberOfBodies (1, 1) double {mustBeInteger, mustBePositive} = 4
        numberOfConstraints (1, 1) double {mustBeInteger, mustBePositive}  = 6
        
        % frequency of physics calculation (per second)
        computeFrequency = 2000
        
        % update frequency of plot
        framerate = 50
        
        % duration of simulation in seconds
        simDuration = 30;
        
        % defines masses of particles, kg
        massesArray = [0.1 0.1 0.1 0.1]
        
        rackMass = 0.5;
        rackLength = 0.5
        
        gearMass = 0.1;
        radGyrationGear = 0.025;
        gearRadius = 0.025;
        I_gear = 0.1
        gearCentre = [0.5 0.025]
        
        gearRotationCount = 0 % used to adjust output of rolling constraint
        % function when gear undergoes full rotation
        
        % defines intial orientation of links, positive is anticlockwise from
        % east
        pendulumAngle = -90
        
        gearParticleAngle = -90.1 % initial orientation, jacobian of rolling 
        % constraint formula blows up if this value is initialized at exactly 90 and -90 deg
        % doesn't seem to attain these values on its own during the simulation
        
        % defines length of links
        lengthArray = [0.3 0.3]
        
        % gravity
        g = -9.81;
        
        Ks = 50; % spring constant 
        Kd = 20; % damping constant
        
        %moveFrequency = 0.9
        timeArray
        
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
        axisLimits
        gearDisp
        rackDisp
        displacementPlot
        gearVertices
        torque
        angleDeltaFromVertical = 0
        angleDeltaModifier % the purpose of this is to return balanced rod to x  = zero position
        lastImpulseForce = 0
        stopMechanism = false
        positionalDelta = 0
        q3Array
        pendulumAngleArr
    end
    
    methods
        
        function obj = run(obj)
            obj = obj.computeConstants;
            obj = obj.calculateJandF1;
            %pendulum_sim.createFigure;
            obj = obj.solver;
        end
        
        function [obj] = createObject(simDuration)
            % createObject; constructs an instance of this class, and runs
            % code
            %clear; close all;
            obj.simDuration = simDuration;
            obj = obj.run;
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
            obj.I_gear = 0.5 * obj.gearMass * obj.radGyrationGear;
            obj.q3Array = zeros(1, obj.simDuration / obj.timestep);
            obj.pendulumAngleArr = zeros(1, obj.simDuration / obj.timestep);
            
            %% Define vector q
            %second particle is pinned to first, first is pinned to rack ???
            %they both make up the rod, third particle is belt system/rack, fourth
            %represents the gear
            obj.q = [obj.lengthArray(1) * cosd(obj.pendulumAngle(1)); obj.lengthArray(1) * sind(obj.pendulumAngle(1)); 0; 0];
            obj.q(5, 1) = obj.q(3, 1) + obj.rackLength; % rack
            obj.q(6, 1) = obj.q(4, 1);
            obj.q(7, 1) = obj.q(3, 1) + obj.rackLength+obj.gearRadius*cosd(obj.gearParticleAngle); % gear
            obj.q(8, 1) = 0+obj.gearRadius+obj.gearRadius*sind(obj.gearParticleAngle); % initial position of gear particle is at ? deg 
                       
            %% Define W
            obj.W = zeros(obj.numbStateParams);
            obj.W(1, 1) = 1/obj.massesArray(2);
            obj.W(2, 2) = 1/obj.massesArray(2);
            obj.W(3, 3) = 1/obj.massesArray(1);
            obj.W(4, 4) = 1/obj.massesArray(1);
            obj.W(5, 5) = 1/obj.rackMass;
            obj.W(6, 6) = 1/obj.rackMass;
            obj.W(7, 7) = 1/obj.gearMass;
            obj.W(8, 8) = 1/obj.gearMass;
            
            %% Define Q

            obj.Q(2) = obj.massesArray(1) * obj.g;
            obj.Q(4) = obj.massesArray(2) * obj.g;
            obj.Q(6) = obj.massesArray(3) * obj.g;
            
            % force on gear particle, zero for now
            %obj.Q(7)
            %obj.Q(8) = 10*sin(2*pi*obj.t*obj.moveFrequency);
        end
        
        function obj = calculateJandF1(obj)
            %% F1 function (denominator)
            obj.F1(1) = sqrt((obj.q(3, 1) - obj.q(1, 1))^2 + (obj.q(4, 1) - obj.q(2, 1))^2); %sqrt((x2-x1)^2 + (y2-y1)^2);
            
            obj.F1(2) = 2*sqrt((obj.q(7, 1) - 0.5)^2 + (obj.q(8, 1) - 0.025)^2); %sqrt((x4-0.5)^2 + (y4-0.025)^2);
            % could've automated the whole differentiation business breh
            
            %% Define J
            obj.J = zeros(obj.numberOfConstraints, obj.numbStateParams);

            obj.J(1, 1) = (obj.q(1) - obj.q(3))/obj.F1(1);
            obj.J(1, 2) = (obj.q(2) - obj.q(4))/obj.F1(1);
            obj.J(1, 3) = (obj.q(3) - obj.q(1))/obj.F1(1);
            obj.J(1, 4) = (obj.q(4) - obj.q(2))/obj.F1(1);
            
            obj.J(2, 3) = -1;
            obj.J(2, 5) = 1;
            
            obj.J(3, 4) = 1;
            
            obj.J(4, 5) = -1; %rolling constraint row in matrix
            
            obj.J(4, 7) = (obj.gearRadius.*(obj.gearRadius - obj.q(8)))./((obj.rackLength - obj.q(7)).^2.*((obj.gearRadius - obj.q(8)).^2./(obj.rackLength - obj.q(7)).^2 + 1));
            obj.J(4, 8) = -obj.gearRadius./((obj.rackLength - obj.q(7)).*((obj.gearRadius - obj.q(8)).^2./(obj.rackLength - obj.q(7)).^2 + 1));
            
            obj.J(5, 6) = 1;
          
            obj.J(6, 7) = (2*obj.q(7) - 1)/obj.F1(2);
            obj.J(6, 8) = (2*obj.q(8) - 0.05)/obj.F1(2);
            
        end
        
        function obj = updateJ(obj)
            %% Define J
            %obj.J = zeros(obj.numberOfConstraints, obj.numbStateParams);
            obj.J(1, 1) = (obj.q(1) - obj.q(3))/obj.F1(1);
            obj.J(1, 2) = (obj.q(2) - obj.q(4))/obj.F1(1);
            obj.J(1, 3) = (obj.q(3) - obj.q(1))/obj.F1(1);
            obj.J(1, 4) = (obj.q(4) - obj.q(2))/obj.F1(1);
            
%             obj.J(2, 3) = -1;
%             obj.J(2, 5) = 1;
%             
%             obj.J(3, 4) = 1;
            
            obj.J(4, 5) = -1; %rolling constraint row in matrix
            
            obj.J(4, 7) = (obj.gearRadius.*(obj.gearRadius - obj.q(8)))./((obj.rackLength - obj.q(7)).^2.*((obj.gearRadius - obj.q(8)).^2./(obj.rackLength - obj.q(7)).^2 + 1));
            obj.J(4, 8) = -obj.gearRadius./((obj.rackLength - obj.q(7)).*((obj.gearRadius - obj.q(8)).^2./(obj.rackLength - obj.q(7)).^2 + 1));
            
            
%             obj.J(5, 6) = 1;
          
            obj.J(6, 7) = (2*obj.q(7) - 1)/obj.F1(2);
            obj.J(6, 8) = (2*obj.q(8) - 0.05)/obj.F1(2);
            if isnan(obj.J(4,8))
                pause(0.5)
                %this seems to happen if pendulum Angle becomes exactly
                %90 or -90
            end 
        end
        
        function obj = angleDeltaModifierPD(obj)
            
            % the angle delta modifier is intended to modifiy the target
            % angle of the inverted rod slightly from the vertical to help
            % return its base (which is particle 2 (q(3) and q(4) are its x and
            % y coordinates respectively) to x = 0; This can be
            % probably easily expanded to bring the rod's position to any
            % desired point. 
            % The mechanism to achieve this is: if rod is in +ve x, tilt the
            % rod slighly to the left so the balancing controller brings it
            % back to zero. The tilt angle of the rod is controlled
            % using a PD controller
            
            % will only implement here when belt/rack is slow
            % (which are particles 2 and 3), to avoid interference from
            % this controller in balancing the rod
            
            if  abs(obj.qDot(3)) < 1 
                proportionalGain = -0.37;
                derivativeGain = -0.3; 
                obj.angleDeltaModifier = proportionalGain * obj.q(1)...
                + derivativeGain * obj.qDot(1);
            else
                %proportionalGain = 0;
                %derivativeGain = 0;
                obj.angleDeltaModifier = 0;
            end
            
            
        % Old cmnts
        % not yet very accurate in the X positioning, tends to overshoot
        % cannot increase derivative gain because it causes jitter in the
        % controller output because it is in the above conditional.
        % Derivative gain appears to rapidly oscillate between zero and its
        % given value.
        
%             if obj.angleDeltaModifier > 0.1  && abs(obj.qDot(3)) < 0.8 
%                 obj.angleDeltaModifier = 0.1;
%             elseif obj.angleDeltaModifier < 0.1 && abs(obj.qDot(3)) < 0.8
%                 obj.angleDeltaModifier = 0.1;
%             end

        end
        
        function obj = motorTorque(obj)
            
            proportionalGain = -2; % needs to be negative?
            derivativeGain = -0.2; % gains have been tuned by a trial-and-error sort of approach
            
            if abs(90-obj.pendulumAngle) <= 12
                % angle deviation is angle difference of pendulum link from
                % the vertical
                obj = angleDeltaModifierPD(obj);

                oldAngleDelta = obj.angleDeltaFromVertical;
                obj.angleDeltaFromVertical = -(pi/2) + atan2( (obj.q(2)-obj.q(4))...
                    , ((obj.q(1)-obj.q(3))) ) + obj.angleDeltaModifier;

                omegaOfRod = (obj.angleDeltaFromVertical - oldAngleDelta)/obj.timestep;

                obj.torque = proportionalGain * obj.angleDeltaFromVertical ...
                    + derivativeGain * omegaOfRod;
            else
                if obj.t <= 5.12
                    freq = 0.78;
                elseif obj.t > 5.12
                    freq = 0.5;
                end
                idleTime = 2;
                targetPos = 0.12 * sin(2*pi*freq*(obj.t-idleTime)); % code to swing up pendulum
                
                proportionalGain = +2; 
                derivativeGain = +0.2;
                
                oldPostionalDelta = obj.positionalDelta;
                obj.positionalDelta = targetPos-obj.q(3);
                errorDot = (obj.positionalDelta - oldPostionalDelta)/obj.timestep;
                obj.torque = proportionalGain * (obj.positionalDelta)...
                    + derivativeGain * errorDot;
                
                if obj.t < idleTime
                    obj.torque = 0;
                end
                
            end
            
            if obj.t >= 24.9
                obj.stopMechanism = true;
            end
            % old comments
            % if pendulum's angle deviation from vertical is greater than
            % 12 deg, controller will cease attempting balance, and attempt
            % to damp out and stop the linear motion of belt /rack system,
            % also condition will be acitvated if controller output motor 
            % torque is greater than 75... (70 Nm) is the intial torque if
            % pendulum angle deviation from the vertical is 10 degrees
            
            % all this tuning and values coded in will likely need to be
            % adjusted for a physical system with a different configuration
            
%             if abs(90-obj.pendulumAngle) > 12 || abs(obj.torque) > 75
%                 obj.stopMechanism = true;
%             end
            
            if obj.stopMechanism == true
                obj.torque = 0;
            end
            
            %stopMechanismBoolean = false;
%             if abs(obj.torque) > 10
%                 obj.torque = 0;
%                 stopMechanismBoolean = true;
%                 warning('Motor output maxed out.')
%             end
            
%             if stopMechanismBoolean == true
%                 obj.torque = derivativeGain * obj.qDot(3);
%             end
        end
       
        function obj = updateGearForces(obj)
            obj = motorTorque(obj);
            %obj.torque = 0.05 * sin(2*pi*2*obj.t);
            forceOnParticle = obj.torque / obj.gearRadius;
            
            % the sign convention for angle here is +ve anti-clockwise from
            % 'east' position of gear particle :)
            
            deltaY = obj.q(8)-obj.gearCentre(2);
            deltaX = obj.q(7)-obj.gearCentre(1);
            gearAngle = atan2(deltaY,deltaX); % 4 quadrant arctan function
            
            obj.Q(7) = forceOnParticle * cos((pi/2)+gearAngle);
            obj.Q(8) = forceOnParticle * sin((pi/2)+gearAngle);
        end
        
        function [Cdot, obj] = calculateCandCdot(obj)
            C_old = obj.C;
            %parameterise these formulae or automate diffing breh
            obj.C(1, 1) = sqrt( (obj.q(3) - obj.q(1))^2+(obj.q(4) - obj.q(2))^2)-0.3;
            obj.C(2, 1) = obj.q(5)-obj.q(3)-obj.rackLength;%-0.1*sin(2*pi*obj.t*obj.moveFrequency);%obj.q(1);
            obj.C(3, 1) = obj.q(4);
            
            % Rolling constraint function stuff
            oldGearParticleAngle = obj.gearParticleAngle;
            obj.gearParticleAngle = (180/pi)* atan2(obj.q(8)-obj.gearCentre(2), obj.q(7)-obj.gearCentre(1));
            obj.pendulumAngle = (180/pi)* atan2(obj.q(2)-obj.q(4), obj.q(1)-obj.q(3));
            
            
            if ((oldGearParticleAngle <= 180 && oldGearParticleAngle >= 90 ) && (obj.gearParticleAngle >= -180 && obj.gearParticleAngle <= -90)) % if gearParticle goes from 2nd to 3rd quadrant
                obj.gearRotationCount = obj.gearRotationCount + 1;
            elseif ((obj.gearParticleAngle <= 180 && obj.gearParticleAngle >= 90 ) && (oldGearParticleAngle >= -180 && oldGearParticleAngle <= -90)) % if gearParticle goes from 3rd to 2nd quadrant
                obj.gearRotationCount = obj.gearRotationCount - 1;
            end
            
            obj.C(4, 1) = obj.gearRadius * (atan2((obj.q(8)-obj.gearCentre(2)),(obj.q(7)-obj.gearCentre(1))) + pi/2 + (obj.gearRotationCount*2*pi)) +obj.rackLength-obj.q(5);
%             if obj.t >= 2.1855
%                 disp('hold')
%             end
            % debugging code
%             if abs(obj.C(4,1)) > 0.001
%                 disp('this formula right above here sir')
%             end
            
            obj.C(5, 1) = obj.q(6);
            obj.C(6, 1) = sqrt( (obj.q(7) - obj.rackLength)^2+(obj.q(8) - obj.gearRadius)^2)-obj.gearRadius;
            
            Cdot = (obj.C - C_old) / obj.timestep;
            
        end
        
        function obj = createDisplacementArrays(obj)   
            obj.gearDisp = 0;%zeros(1, 200);
            obj.rackDisp = 0;%zeros(1, 200);
            obj.timeArray= 0;%zeros(1, 200);
            figure(2)
            obj.displacementPlot = plot(obj.timeArray, obj.gearDisp, 'o-.g', obj.timeArray, obj.rackDisp, '*-.r');
            legend('Gear particle disp','Rack particle disp')
            
            obj.gearDisp = [];
            obj.rackDisp = [];
            obj.timeArray= [];% this twister game has to be played so that the legend gets displayed
            % earlier a value was stored in these arrays to get the legend
            % to display when they are plotted
            % note that this plotting business will slow down the code
            % execution
        end
        
        
        function obj = plotDisplacements(obj)
            
            if size(obj.timeArray,1) >= 200
                obj.gearDisp = circshift(obj.gearDisp, -1);
                obj.rackDisp = circshift(obj.rackDisp, -1);
                obj.timeArray= circshift(obj.timeArray, -1);
                
                %all values get shifted back and last element gets replaced
                %with current value
                obj.gearDisp(end) = obj.gearRadius * (obj.gearParticleAngle*(pi/180) -(atan(obj.q(8)-obj.gearCentre(2)/obj.q(7)-obj.gearCentre(1))));
                obj.rackDisp(end) = obj.q(5)-obj.rackLength;
                obj.timeArray(end) = obj.t;
            else
                obj.gearDisp = [obj.gearDisp (obj.gearRadius * (obj.gearParticleAngle*(pi/180) -(atan(obj.q(8)-obj.gearCentre(2)/obj.q(7)-obj.gearCentre(1)))))];
                obj.rackDisp = [obj.rackDisp (obj.q(5)-obj.rackLength)];
                obj.timeArray = [obj.timeArray obj.t];
            end
            %gearDisp contains arc length of gear particle from vertically
            %down position
            % rackDisp contains displacement of rack particle from x=0.5
            set(obj.displacementPlot(1),'XData',obj.timeArray,'YData',obj.gearDisp);
            set(obj.displacementPlot(2),'XData',obj.timeArray,'YData',obj.gearDisp);
            
        end
        
        function obj = computeGearVertices(obj)
            var = obj.gearParticleAngle;
            %gear is depicted as a hexagon
            gearVerticesAngles = [var+60 var+120 var+180 var+240 var+300 var];
            %[var+120 var+240 var];
            
            %first row contains x values, second contains y values
            obj.gearVertices(1, :) = obj.gearRadius * cosd(gearVerticesAngles) + obj.gearCentre(1);
            obj.gearVertices(2, :) = obj.gearRadius * sind(gearVerticesAngles) + obj.gearCentre(2);
        end
        
        function obj = randomImpulseGenerator(obj)
            % old comments: every x seconds, a random impulse is calculated to be imparted to the
            % pendulum
            % impulse will be zero if stopMechanism is true or if time (t)
            % is zero
            % if false
             
%             if obj.t >= 15
%                 obj.t
%             end
                maxImpulseForce = 30;
                if (obj.t >= 10 && obj.t <= (10+0.5*obj.timestep)) && obj.stopMechanism == false
                    obj.lastImpulseForce = -50 - maxImpulseForce * rand; % random number in range
                    % of [-maxImpulse, maxImpulse] newtons

                    %add last impulse to x force on top of pendulum
                    obj.Q(1) = obj.Q(1) + obj.lastImpulseForce;
                    
                    % below condition will impart impulse at 20, 20.6 and
                    % 21.2 seconds. Have to implement the in btween
                    % condition because obj.t value seems to have some
                    % small error in it, it is not precisely the expected
                    % value e.g.15.00000000000000000000000
                elseif ((obj.t >= 15 && obj.t <= (15+0.5*obj.timestep))...
                        || (obj.t >= 15.8 && obj.t <= (15.8+0.5*obj.timestep))...
                        || (obj.t >= 16.6 && obj.t <= (16.6+0.5*obj.timestep)))...
                        && obj.stopMechanism == false
                    obj.lastImpulseForce = 50+ 0.8 * maxImpulseForce * rand;
                    obj.Q(1) = obj.Q(1) + obj.lastImpulseForce;
                else
                    obj.Q(1) = 0;% resets x force on particle 1

                end
            %end
        end
        
        function obj = frictionForceGenerator(obj)
            % for a damping effect
            dynamicFrictionCoefficient = 0.02;
            obj.Q(3) = obj.qDot(3) * -dynamicFrictionCoefficient;
            
            % friction gets added to any forces from impulse generator
            obj.Q(1) = obj.Q(1) + obj.qDot(1) * -dynamicFrictionCoefficient;
            obj.Q(2) = obj.Q(2) + obj.qDot(2) * -dynamicFrictionCoefficient;
        end
        
        function obj = reset_Q_forces(obj)
            % used to reset forces on particles to prevent accumulation in
            % between timesteps
            obj.Q = zeros(obj.numbStateParams, 1);
            obj.Q(2) = obj.massesArray(1) * obj.g;
            obj.Q(4) = obj.massesArray(2) * obj.g;
            obj.Q(6) = obj.massesArray(3) * obj.g;
        end
        
        function obj = solver(obj)
            %
            %% Creates figure
            % creates lines
            obj = computeGearVertices(obj);
            var = obj.gearVertices;
            
            figure('Position', [10 50 1510 730]);
            animLine1 = animatedline([obj.q(1) obj.q(3) obj.q(5)],...
                [obj.q(2) obj.q(4) obj.q(6)], 'LineWidth', 5, ...
                'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'w');
            gearLine = animatedline([obj.rackLength obj.q(7) var(1,:)],...
                [obj.gearRadius obj.q(8) var(2,:)], 'LineWidth', 2, ...
                'Marker', 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'w');
            
            %animLine2 = animatedline([)], [], 'LineWidth', 10, 'Marker', 'o', 'MarkerSize', 9, 'MarkerFaceColor', 'w');
            axis equal
            xlim([-0.8 1.5])
            ylim([-0.5 0.5])
            pause(0.5)
            
            %% Solver
            
            %obj = createDisplacementArrays(obj);
            
            i = 0;
            totalIterationCount = obj.computeFrequency * obj.simDuration +1;
            %intialize time to zero
            obj.t = -obj.timestep; % to start t from zero in first iteration
            %lamda_old = -obj.Jdot * obj.qDot - obj.J * obj.W * obj.Q - (obj.Ks * obj.C) - (obj.Kd * 0);
            
            %loopRateControl = rateControl(obj.computeFrequency);
            %reset(loopRateControl)
            loopTimeReduction = 0.5453*0.67; % loop takes too long, hence this factor
            totalTime = tic;
            while i <= totalIterationCount
                loopTime = tic;
                i = i+1;
                obj.t = obj.t + obj.timestep;
                %q_old = obj.q;
                [Cdot, obj] = calculateCandCdot(obj);
                
                %obj.Q(1) = 10*sin(2*pi*obj.t*obj.moveFrequency);
                obj = updateGearForces(obj);
                if obj.t==15
                    pause(4)
                end
                obj = randomImpulseGenerator(obj);
                obj = frictionForceGenerator(obj);
                                
                if obj.t > 19.99
                    obj.stopMechanism = true;
                end
                
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
                
                %storing position of base of pendulum with time
                obj.q3Array(i) = obj.q(3);
                %storing value in array
                obj.pendulumAngleArr(i) = obj.pendulumAngle;
                
                J_old = obj.J;
                
                obj = updateJ(obj);
                if isnan(obj.J(4,8))
                    pause(0.5)
                    %J goes NaN
                end %debugging code
                
                obj.Jdot = (obj.J - J_old) ./ obj.timestep;
                
                obj = reset_Q_forces(obj);
                
                if mod(obj.t, obj.frameTimestep) <= (obj.timestep)
                    clearpoints(animLine1)
                    addpoints(animLine1, [obj.q(1) obj.q(3) obj.q(5)], [obj.q(2) obj.q(4) obj.q(6)])
                    clearpoints(gearLine)
                    obj = computeGearVertices(obj);
                    var = obj.gearVertices;
                    addpoints(gearLine, [obj.rackLength obj.q(7) var(1,:)],...
                        [obj.gearRadius obj.q(8) var(2,:)])
                    drawnow
                    %potentialEnergy = obj.massesArray(1) * -obj.g * (obj.q(2) + obj.q(4))+ obj.rackMass * -obj.g * obj.q(6);% + obj.q(8));
                    %kineticEnergy = 0.5 * obj.massesArray(1)*( obj.qDot(1)^2+obj.qDot(2)^2 + obj.qDot(3)^2+obj.qDot(4)^2) + 0.5 * obj.rackMass * (obj.qDot(5)^2+obj.qDot(6)^2) + 0.5 * obj.gearMass * (obj.qDot(7)^2+obj.qDot(8)^2);
                    figure(1)
                    titleText = sprintf('Time: %.2f s     Last impulse = %.4f N.s     Gear Angle = %.1f deg \n\nPendulum free end X position = %.3f m    Pendulum Angle = %.2f deg',...
                        obj.t, ...
                        obj.lastImpulseForce*obj.timestep, obj.gearParticleAngle, obj.q(1), obj.pendulumAngle);
                    %    Energy: %.3f J    Torque = %.3f Nm , (potentialEnergy+kineticEnergy), obj.torque,
                    
                    title(titleText)
                    
                    %obj = plotDisplacements(obj);
                    %disp('s')% this line is just for a breakpoint
                    
                end
                
                %waitfor(loopRateControl); not needed
                pause((obj.timestep - toc(loopTime))*loopTimeReduction)
                %toc(loopTime)
            end
            
            toc(totalTime)
            
        end
        
    end
    
end
