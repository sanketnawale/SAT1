clear all; close all; clc;

%% 🌍 Simulation Parameters
tic
Simulation_T = 1000;  % Total simulation time (seconds)
Time_Step = 10;        % Time step duration (seconds)
MonteCarlo = 5000;    % Monte Carlo simulation (number of packets per node)
Nodes = 3;            % Number of ground nodes (Rome, Milan, NodeRM)
Pkct_Per_Hour = 100;   % Packets per hour for each node

%% 🌍 Ground Nodes (Rome, Milan, NodeRM)
Node_Coordinates = [
    41.9028, 12.4964;  % Rome
    45.4642, 9.1900;   % Milan
    41.9, 12.5         % NodeRM (close to Rome)
];

%% 🛰️ Satellite Constellation (Walker)
Num_Satellites = 2;
Num_Planes = 6;
Orbital_Inclination = deg2rad(87); % Inclination in Radians
H = 1200e3; % Satellite altitude (meters)
Earth_Radius = 6378e3; % Earth radius (meters)
Time_Vector = 0:Time_Step:Simulation_T; % Time simulation vector

% 🛰️ Generate Walker Delta Constellation
oev = walker_delta(Num_Satellites, Num_Planes, 1, pi, Earth_Radius + H, Orbital_Inclination);

%% 📡 Call Updated Satellite Geometry Function
[Distances, Elevation_Angles, Ground_Distances, Visibility, Num_Visible_Sats] = ...
    Satellite_Geometry(H, Node_Coordinates, oev, Earth_Radius, Time_Vector);

%% 📡 Random Access Logic for Packet Transmission, Reception & Collision Detection

% Initialize Success Rates and Collisions
SuccessRate = zeros(Nodes, length(Time_Vector));  % Success rate per node per time step
Collisions = zeros(Nodes, length(Time_Vector));   % Collision counts per node per time step

% Initialize Received Packets for NodeRM
Received_Packets_NodeRM = zeros(1, length(Time_Vector));

% Simulate Random Access and Reception Logic
for t = 1:length(Time_Vector)
    fprintf('\n⏳ Time %.2f sec: \n', Time_Vector(t));

    % Tracking packets received by NodeRM
    NodeRM_Packet_Times = [];

    for n = 1:Nodes-1  % Only Rome (Node 1) and Milan (Node 2) transmit
        if Num_Visible_Sats(n, t) > 1  % If the node sees at least one satellite
            % Generate Random Packet Transmission Times
            Transmission_Times = rand(1, MonteCarlo) * Time_Step; % Packets randomly sent within time step
            sorted_times = sort(Transmission_Times);

            % Select visible satellites
            Visible_Sats = find(Visibility(n, :, t));

            % 🛠 **Fix for Error: Only proceed if satellites are visible**
            if isempty(Visible_Sats)
                continue; % No visible satellites, skip transmission
            end

            % Initialize Satellite Reception Time Storage
            Sat_Receive_Times = cell(Num_Satellites, 1);
            for s = 1:Num_Satellites
                Sat_Receive_Times{s} = [];
            end

            % Assign packets to visible satellites
            for pkt = 1:MonteCarlo
                chosen_sat = Visible_Sats(randi(length(Visible_Sats))); % Random satellite
                if chosen_sat <= Num_Satellites  % 🛠 **Fix for Out of Bounds Error**
                    Sat_Receive_Times{chosen_sat} = [Sat_Receive_Times{chosen_sat}, sorted_times(pkt)];
                end
            end

            % 🔥 Check for Collisions at Each Satellite
            for s = 1:Num_Satellites
                if ~isempty(Sat_Receive_Times{s})
                    sat_tx_times = sort(Sat_Receive_Times{s});
                    
                    % Detect collisions (packets arriving too close)
                    collisions = sum(diff(sat_tx_times) < 0.01);  % Collision if packets arrive within 10ms
                    total_packets = length(sat_tx_times);
                    
                    % Store results
                    Collisions(n, t) = Collisions(n, t) + collisions; 
                    SuccessRate(n, t) = SuccessRate(n, t) + (total_packets - collisions);

                    % 🚀 If node is Rome, relay to NodeRM
                    if n == 1  % Rome transmits
                        % NodeRM receives all packets from visible satellites
                        NodeRM_Packet_Times = [NodeRM_Packet_Times, sat_tx_times];
                    end
                end
            end

            % Print Debug Info
            fprintf('Node %d transmitted %d packets, %d collisions\n', n, MonteCarlo, Collisions(n, t));
        end
    end

    % 🚀 **NodeRM Packet Reception Logic**
    if ~isempty(NodeRM_Packet_Times)
        NodeRM_Packet_Times = sort(NodeRM_Packet_Times);

        % **First packet received is successful, rest are discarded**
        Received_Packets_NodeRM(t) = 1;  % Mark as successful reception
    end

    % Print Success Rates
    fprintf('Node 1 (Rome) Success Rate: %.2f%%\n', SuccessRate(1, t) / MonteCarlo * 100);
    fprintf('Node 2 (Milan) Success Rate: %.2f%%\n', SuccessRate(2, t) / MonteCarlo * 100);
    fprintf('NodeRM Received Packet: %d\n', Received_Packets_NodeRM(t));
end

%% 📊 Results Summary
fprintf('\n==== 📊 Final Results ====\n');
fprintf('Overall Success Rate for Rome: %.2f%%\n', mean(SuccessRate(1, :)) / MonteCarlo * 100);
fprintf('Overall Success Rate for Milan: %.2f%%\n', mean(SuccessRate(2, :)) / MonteCarlo * 100);
fprintf('Total Packets Successfully Received by NodeRM: %d\n', sum(Received_Packets_NodeRM));

toc;
