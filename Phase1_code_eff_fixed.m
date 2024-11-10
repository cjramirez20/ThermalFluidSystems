%%%%%% Base Case %%%%%%%%%%%%%

% Define constants
R = 8.314; % kJ/(kmol*K), gas constant
M = 28.84; % kg/kmol, molar mass of air

% Given state variables (convert to SI units)
m_dot = 86.0464726; % kg/s, mass flow rate
T1 = 18.33 + 273.15; % Convert Celsius to Kelvin
P1 = 99.3; % kPa
P0 = 99.3; % Reference pressure (kPa)
V1 = 72.365; % m^3/s (see hand calc for details)

% Call Calculator function for State 1
[u1, h1, s1_0, s1] = thermo_properties(P1, T1);

% State 2
P2 = (99.3 - 0.996); % kPa, Pressure at state 2
T2 = T1; % Temperature remains the same
[u2, h2, s2_0, s2] = thermo_properties(P2, T2);

% State 25s
PRC1 = 6; % Compression ratio
P25s = P2 * PRC1; % kPa, Pressure at state 25s
s25s = s2; % Known entropy at state 25s
s25s0 = s25s + (R/M) * log(P25s / P0);
T25s = get_temperature(s25s0, 's0', P25s); % Solve temperature for State 25s
[u25s, h25s, s25s0, s25s] = thermo_properties(P25s, T25s);

% State 25 (Efficiency)
eta_25 = 0.82; % Efficiency of state 25
h25 = h2 - (h2 - h25s)/eta_25; % kj/kg
P25 = P25s;
T25 = get_temperature(h25, 'h', P25); % Solve temperature based on enthalpy at State 25
[u25, h25, s25_0, s25] = thermo_properties(P25, T25);
w_c1=m_dot*(h1-h25);

% State 3s
PRC2 = 4;
P3s = P25 * PRC2; % kPa, Pressure at state 3s
s3s = s25; % Known entropy at state 3s
s3s0 = s3s + (R/M) * log(P3s / P0);
T3s = get_temperature(s3s0, 's0', P3s); % Solve temperature for State 3s
[u3s, h3s, s3s_0, s3s] = thermo_properties(P3s, T3s);

% State 3 (Efficiency)
eta_3 = 0.84; % Efficiency of state 3
h3 = h25 - (h25 - h3s)/eta_3; % kj/kg
P3 = P25 * 4; % Pa, Pressure at state 3
T3 = get_temperature(h3, 'h', P3); % Solve temperature based on enthalpy at State 3
[u3, h3, s3_0, s3] = thermo_properties(P3, T3);
w_c2 = m_dot*(h25-h3);

% State 4 (After fuel addition)
m_fuel_dot = 1.7841; % kg/s, mass flow rate of fuel
LHV = 46950.31; % kj/kg, Lower Heating Value of fuel
Q_fuel = LHV * m_fuel_dot; % Heat input from fuel
m4_dot = m_dot + m_fuel_dot; % Total mass flow rate at state 4
h4 = (Q_fuel + m_dot * h3) / m4_dot; % kj/kg
P4 = P3; % Pressure stays the same
T4 = get_temperature(h4, 'h', P4); % Solve temperature based on enthalpy at State 4
[u4, h4, s4_0, s4] = thermo_properties(P4, T4);

% State 48s
P48s = 489.53; % kPa, Pressure at state 48s
s48s = s4; % Known entropy at state 48s
s48s0 = s48s + (R/M) * log(P48s / P0);
T48s = get_temperature(s48s0, 's0', P48s); % Solve temperature for State 48s
[u48s, h48s, s48s_0, s48s] = thermo_properties(P48s, T48s);

% State 48 (Efficiency)
w_t1 = abs(w_c1) + abs(w_c2);
PR_t1 = 4.824;
eta_t1 = 0.88; % Efficiency of state 3
h48 = h4 - (h4 - h48s)*eta_t1; % kj/kg
P48 = P4 * 4.824; % Pa, Pressure at state 3
T48 = get_temperature(h48, 'h', P48); % Solve temperature based on enthalpy at State 3
[u48, h48, s48_0, s48] = thermo_properties(P48, T48);

% State 6
T6=T1;
P6=P0;
[u6, h6, s6_0, s6] = thermo_properties(P6, T6);

% State 5s
P5s = P1 + (2.491); % kPa, Pressure at state 48s
s5s = s48; % Known entropy at state 48s
s5s0 = s5s + (R/M) * log(P5s / P0);
T5s = get_temperature(s5s0, 's0', P5s); % Solve temperature for State 48s
[u5s, h5s, s5s_0, s5s] = thermo_properties(P5s, T5s);

% State 5
P5 = P1 + (2.491);
PR_t2 = 4.804;
eta_t2 = 0.821; % Efficiency of state 3
h5 = h48 - (h48 - h5s)*eta_t2; % kj/kg
T5 = get_temperature(h5, 'h', P5); % Solve temperature based on enthalpy at State 3
[u5, h5, s5_0, s5] = thermo_properties(P5, T5);
w_t2 = m4_dot*abs(h5-h48);
eta_gen = w_t2/30607;


% Create a table with the relevant information
State = {'Inlet (State 1)', 'State 2', 'State 25', 'State 3', 'State 4', 'State 48', 'State 5', 'State 6'};
Temperature_K = [T1-273.15, T2-273.15, T25-273.15, T3-273.15, T4-273.15, T48-273.15, T5-273.15, T6-273.15];
Pressure_Pa = [P1, P2, P25, P3, P4, P48, P5, P6];
VolumetricFlowRate_Inlet = [V1, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
TurbineEfficiency = [NaN, NaN, NaN, NaN, NaN, eta_t1, eta_t2, NaN]; % Only for relevant turbine states

% Combine into a table
T = table(State', Temperature_K', Pressure_Pa', VolumetricFlowRate_Inlet', TurbineEfficiency', ...
          'VariableNames', {'State', 'Temperature (C)', 'Pressure (kPa)', 'Volumetric Flow Rate (m^3/s)', 'Turbine Efficiency'});

% Display the table
disp(T);

%%%%%%%%%%%%%%%%%%%%% various fuel conditions %%%%%%%%%%%%%%%%%%%%%%%%%
fuel_percentage = [.2,.4,.6,.8,1];

% temperature
temperature_fahrenheit = 65;
temperature_kelvin = 273.15 + (temperature_fahrenheit - 32) * 5 / 9;

% coresponding fuel variations
fuel_flow_lb_per_hr = fuel_percentage * 14585; % fuel mass flow based on base case
fuel_flow_kg_per_s = fuel_flow_lb_per_hr * 0.453592 / 3600; % convert to kg/s

% Given
eta_t1 = 0.88;
eta_t2 = 0.821;

% declare output variables for plotting
P_net = [NaN,NaN,NaN,NaN];
m_dot_inlet=[NaN,NaN,NaN,NaN];
m_dot_outlet=[NaN,NaN,NaN,NaN];
thermal_efficiency=[NaN,NaN,NaN,NaN];
T_outlet=[NaN,NaN,NaN,NaN];
T_inlet=[NaN,NaN,NaN,NaN];
SFC=[NaN,NaN,NaN,NaN];
Heat_rate=[NaN,NaN,NaN,NaN];

for i = 1:length(fuel_flow_kg_per_s)
    T1 = temperature_kelvin; % temperature for specific iteration
    specific_volume = (R/M)*T1/(P1); %%% update this %%%
    m_dot = V1/specific_volume; % kg/s, mass flow rate
    
    % Call Calculator function for State 1
    [u1, h1, s1_0, s1] = thermo_properties(P1, T1);
    
    % State 2
    P2 = (99.3 - 0.996); % kPa, Pressure at state 2
    T2 = T1; % Temperature remains the same
    [u2, h2, s2_0, s2] = thermo_properties(P2, T2);
    
    % State 25s
    PRC1 = 6; % Compression ratio
    P25s = P2 * PRC1; % kPa, Pressure at state 25s
    s25s = s2; % Known entropy at state 25s
    s25s0 = s25s + (R/M) * log(P25s / P0);
    T25s = get_temperature(s25s0, 's0', P25s); % Solve temperature for State 25s
    [u25s, h25s, s25s0, s25s] = thermo_properties(P25s, T25s);
    
    % State 25 (Efficiency)
    eta_25 = 0.82; % Efficiency of state 25
    h25 = h2 - (h2 - h25s)/eta_25; % kj/kg
    P25 = P25s;
    T25 = get_temperature(h25, 'h', P25); % Solve temperature based on enthalpy at State 25
    [u25, h25, s25_0, s25] = thermo_properties(P25, T25);
    w_c1=m_dot*(h1-h25);
    
    % State 3s
    PRC2 = 4;
    P3s = P25 * PRC2; % kPa, Pressure at state 3s
    s3s = s25; % Known entropy at state 3s
    s3s0 = s3s + (R/M) * log(P3s / P0);
    T3s = get_temperature(s3s0, 's0', P3s); % Solve temperature for State 3s
    [u3s, h3s, s3s_0, s3s] = thermo_properties(P3s, T3s);
    
    % State 3 (Efficiency)
    eta_3 = 0.84; % Efficiency of state 3
    h3 = h25 - (h25 - h3s)/eta_3; % kj/kg
    P3 = P25 * 4; % Pa, Pressure at state 3
    T3 = get_temperature(h3, 'h', P3); % Solve temperature based on enthalpy at State 3
    [u3, h3, s3_0, s3] = thermo_properties(P3, T3);
    w_c2 = m_dot*(h25-h3);
    
    % State 4 (After fuel addition)
    m_fuel_dot = fuel_flow_kg_per_s(i); % kg/s, mass flow rate of fuel
    LHV = 46950.31; % kj/kg, Lower Heating Value of fuel
    Q_fuel = LHV * m_fuel_dot; % Heat input from fuel
    m4_dot = m_dot + m_fuel_dot; % Total mass flow rate at state 4
    h4 = (Q_fuel + m_dot * h3) / m4_dot; % kj/kg
    P4 = P3; % Pressure stays the same
    T4 = get_temperature(h4, 'h', P4); % Solve temperature based on enthalpy at State 4
    [u4, h4, s4_0, s4] = thermo_properties(P4, T4);
    
    % State 48s
    P48s = 489.53; % kPa, Pressure at state 48s
    s48s = s4; % Known entropy at state 48s
    s48s0 = s48s + (R/M) * log(P48s / P0);
    T48s = get_temperature(s48s0, 's0', P48s); % Solve temperature for State 48s
    [u48s, h48s, s48s_0, s48s] = thermo_properties(P48s, T48s);
    
    % State 48 (Efficiency)
    w_t1 = abs(w_c1) + abs(w_c2);
    % h48 = h4 - w_t1/m4_dot; % kj/kg
    h48=h4-(eta_t1*(h4-h48s)); % uncomment this if you want to use base
    % case efficiency instead
    P48 = P48s; % Pa, given pressure
    T48 = get_temperature(h48, 'h', P48); % Solve temperature based on enthalpy at State 48
    [u48, h48, s48_0, s48] = thermo_properties(P48, T48);
    % eta_t1_array(i) = (h4-h48)/(h4-h48s);
    PR_t1 = P48/P4;
    
    % State 6
    T6=T1;
    P6=P0;
    [u6, h6, s6_0, s6] = thermo_properties(P6, T6);
    
    % State 5s
    P5s = P1 + (2.491);
    s5s = s48; % Known entropy at state 5s
    s5s0 = s5s + (R/M) * log(P5s / P0);
    T5s = get_temperature(s5s0, 's0', P5s); % Solve temperature for State 5s
    [u5s, h5s, s5s_0, s5s] = thermo_properties(P5s, T5s);

    % State 5
    P5 = P1 + (2.491);
    PR_t2 = 4.804;
    eta_5 = 0.821; % Efficiency of state 3
    h5 = h48 - (h48 - h5s)*eta_5; % kj/kg
    T5 = get_temperature(h5, 'h', P5); % Solve temperature based on enthalpy at State 3
    [u5, h5, s5_0, s5] = thermo_properties(P5, T5);
    w_t2 = m4_dot*abs(h5-h48);
    eta_gen = w_t2/30607;
        
    % Calculate requested values (in SI)
    P_net(i) = w_t2*eta_gen;
    m_dot_inlet(i) = m_dot;
    m_dot_outlet(i)=m4_dot;
    T_outlet(i) = T5;
    T_inlet(i)=T48;
    SFC(i) = m_fuel_dot / P_net(i);
    Heat_rate(i) = SFC(i) * LHV;
    thermal_efficiency(i) = P_net(i)/(m_fuel_dot*LHV);

end

m_dot_outlet_lbm_hr = m_dot_outlet * 7936.64; % kg/s to lbm/hr
m_dot_inlet_lbm_hr = m_dot_inlet * 7936.64; % kg/s to lbm/hr
fuel_flow_lbm_hr = fuel_flow_kg_per_s * 7936.64; % kg/s to lbm/hr

% Plot each variable on its own figure
figure;
plot(fuel_flow_lbm_hr, P_net, '-o');
xlabel('Fuel Flow Rate (lbm/hr)');
ylabel('Net Power (kW)');
title('Fuel Flow Rate vs Net Power');

figure;
plot(fuel_flow_lbm_hr, m_dot_inlet_lbm_hr, '-o');
xlabel('Fuel Flow Rate (lbm/hr)');
ylabel('Inlet Mass Flow Rate (lbm/hr)');
title('Fuel Flow Rate vs Inlet Mass Flow Rate');

figure;
plot(fuel_flow_lbm_hr, m_dot_outlet_lbm_hr, '-o');
xlabel('Fuel Flow Rate (lbm/hr)');
ylabel('Outlet Mass Flow Rate (lbm/hr)');
title('Fuel Flow Rate vs Outlet Mass Flow Rate');

figure;
plot(fuel_flow_lbm_hr, (T_outlet - 273.15) * 9/5 + 32, '-o');
xlabel('Fuel Flow Rate (lbm/hr)');
ylabel('Turbine Outlet Temperature (°F)');
title('Fuel Flow Rate vs Turbine Outlet Temperature (State 5)');

figure;
plot(fuel_flow_lbm_hr, (T_inlet - 273.15) * 9/5 + 32, '-o');
xlabel('Fuel Flow Rate (lbm/hr)');
ylabel('Turbine Inlet Temperature (°F)');
title('Fuel Flow Rate vs Turbine Inlet Temperature (State 48)');

figure;
plot(fuel_flow_lbm_hr, SFC*2.205*60, '-o');
xlabel('Fuel Flow Rate (lbm/hr)');
ylabel('Specific Fuel Consumption (lbm/kWh)');
title('Fuel Flow Rate vs Specific Fuel Consumption');

figure;
plot(fuel_flow_lbm_hr, Heat_rate*.9478*60, '-o'); % heat rate kJ/
xlabel('Fuel Flow Rate (lbm/hr)');
ylabel('Heat Rate (BTU/kWh)');
title('Fuel Flow Rate vs Heat Rate');

figure;
plot(fuel_flow_lbm_hr, thermal_efficiency, '-o');
xlabel('Fuel Flow Rate (lbm/hr)');
ylabel('Thermal Efficiency');
title('Fuel Flow Rate vs Thermal Efficiency');

%%%%%%%%%%%%%%%%%%%%%%%%% various temperatures %%%%%%%%%%%%%%%%%%%%%%%%%%
% temperature variations
temperatures_fahrenheit = [35, 65, 85, 105];
temperatures_kelvin = 273.15 + (temperatures_fahrenheit - 32) * 5 / 9;

% coresponding fuel variations
fuel_flow_lb_per_hr = [15981, 14585, 13518, 12410];
fuel_flow_kg_per_s = fuel_flow_lb_per_hr * 0.453592 / 3600;

% Given state variables (convert to SI units)
eta_t1 = 0.8302;
eta_t2 = 0.6592;

% declare output variables for plotting
P_net = [NaN,NaN,NaN,NaN];
m_dot_inlet=[NaN,NaN,NaN,NaN];
m_dot_outlet=[NaN,NaN,NaN,NaN];
thermal_efficiency=[NaN,NaN,NaN,NaN];
T_outlet=[NaN,NaN,NaN,NaN];
T_inlet=[NaN,NaN,NaN,NaN];
SFC=[NaN,NaN,NaN,NaN];
Heat_rate=[NaN,NaN,NaN,NaN];

for i = 1:length(temperatures_kelvin)
    T1 = temperatures_kelvin(i); % temperature for specific iteration
    specific_volume = (R/M)*T1/(P1); % using ideal gas
    m_dot = V1/specific_volume; % kg/s, mass flow rate
    
    % Call Calculator function for State 1
    [u1, h1, s1_0, s1] = thermo_properties(P1, T1);
    
    % State 2
    P2 = (99.3 - 0.996); % kPa, Pressure at state 2
    T2 = T1; % Temperature remains the same
    [u2, h2, s2_0, s2] = thermo_properties(P2, T2);
    
    % State 25s
    PRC1 = 6; % Compression ratio
    P25s = P2 * PRC1; % kPa, Pressure at state 25s
    s25s = s2; % Known entropy at state 25s
    s25s0 = s25s + (R/M) * log(P25s / P0);
    T25s = get_temperature(s25s0, 's0', P25s); % Solve temperature for State 25s
    [u25s, h25s, s25s0, s25s] = thermo_properties(P25s, T25s);
    
    % State 25 (Efficiency)
    eta_25 = 0.82; % Efficiency of state 25
    h25 = h2 - (h2 - h25s)/eta_25; % kj/kg
    P25 = P25s;
    T25 = get_temperature(h25, 'h', P25); % Solve temperature based on enthalpy at State 25
    [u25, h25, s25_0, s25] = thermo_properties(P25, T25);
    w_c1=m_dot*(h1-h25);
    
    % State 3s
    PRC2 = 4;
    P3s = P25 * PRC2; % kPa, Pressure at state 3s
    s3s = s25; % Known entropy at state 3s
    s3s0 = s3s + (R/M) * log(P3s / P0);
    T3s = get_temperature(s3s0, 's0', P3s); % Solve temperature for State 3s
    [u3s, h3s, s3s_0, s3s] = thermo_properties(P3s, T3s);
    
    % State 3 (Efficiency)
    eta_3 = 0.84; % Efficiency of state 3
    h3 = h25 - (h25 - h3s)/eta_3; % kj/kg
    P3 = P25 * 4; % Pa, Pressure at state 3
    T3 = get_temperature(h3, 'h', P3); % Solve temperature based on enthalpy at State 3
    [u3, h3, s3_0, s3] = thermo_properties(P3, T3);
    w_c2 = m_dot*(h25-h3);
    
    % State 4 (After fuel addition)
    m_fuel_dot = fuel_flow_kg_per_s(i); % kg/s, mass flow rate of fuel
    LHV = 46950.31; % kj/kg, Lower Heating Value of fuel
    Q_fuel = LHV * m_fuel_dot; % Heat input from fuel
    m4_dot = m_dot + m_fuel_dot; % Total mass flow rate at state 4
    h4 = (Q_fuel + m_dot * h3) / m4_dot; % kj/kg
    P4 = P3; % Pressure stays the same
    T4 = get_temperature(h4, 'h', P4); % Solve temperature based on enthalpy at State 4
    [u4, h4, s4_0, s4] = thermo_properties(P4, T4);
    
    % State 48s
    P48s = 489.53; % kPa, Pressure at state 48s
    s48s = s4; % Known entropy at state 48s
    s48s0 = s48s + (R/M) * log(P48s / P0);
    T48s = get_temperature(s48s0, 's0', P48s); % Solve temperature for State 48s
    [u48s, h48s, s48s_0, s48s] = thermo_properties(P48s, T48s);
    
    % State 48 (Efficiency)
    w_t1 = abs(w_c1) + abs(w_c2);
    h48 = h4 - w_t1/m4_dot; % kj/kg
    % h48=h4-(eta_t1*(h4-h48s)); % uncomment this if you want to use base
    % case efficiency instead
    P48 = P48s; % Pa, given pressure
    T48 = get_temperature(h48, 'h', P48); % Solve temperature based on enthalpy at State 48
    [u48, h48, s48_0, s48] = thermo_properties(P48, T48);
    eta_t1_array(i) = (h4-h48)/(h4-h48s);
    PR_t1 = P48/P4;
    
    % State 6
    T6=T1;
    P6=P0;
    [u6, h6, s6_0, s6] = thermo_properties(P6, T6);
    
    % State 5s
    P5s = P1 + (2.491);
    s5s = s48; % Known entropy at state 5s
    s5s0 = s5s + (R/M) * log(P5s / P0);
    T5s = get_temperature(s5s0, 's0', P5s); % Solve temperature for State 5s
    [u5s, h5s, s5s_0, s5s] = thermo_properties(P5s, T5s);

    % State 5
    eta_gen=.977;
    P5 = P1 + (2.491); % kPa, Pressure at state 5 (after pressure drop)
    h5=h48-(eta_t2*(h48-h5s));
    w_t2=m4_dot*(h48-h5);
    T5 = get_temperature(h5, 'h', P5); % Solve temperature for State 5
    [u5, h5, s5_0, s5] = thermo_properties(P5, T5);
    
    % Calculate requested values (in SI)
    P_net(i) = w_t2*eta_gen;
    m_dot_inlet(i) = m_dot;
    m_dot_outlet(i)=m4_dot;
    T_outlet(i) = T5;
    T_inlet(i)=T48;
    SFC(i) = m_fuel_dot / (P_net(i));
    Heat_rate(i) = SFC(i) * LHV;
    thermal_efficiency(i) = P_net(i)/(m_fuel_dot*LHV);

end

m_dot_outlet_lbm_hr = m_dot_outlet * 7936.64; % convert from kg/s to lbm/hr
m_dot_inlet_lbm_hr = m_dot_inlet * 7936.64; %  convert from kg/s to lbm/hr

% Plot each variable on its own figure
figure;
plot(temperatures_fahrenheit, P_net, '-o');
xlabel('Inlet Temperature (°F)');
ylabel('Net Power (kW)');
title('Inlet Temperature (°F) vs Net Power');

figure;
plot(temperatures_fahrenheit, m_dot_inlet_lbm_hr, '-o');
xlabel('Inlet Temperature (°F)');
ylabel('Inlet Mass Flow Rate (lbm/hr)');
title('Inlet Temperature (°F) vs Inlet Mass Flow Rate');

figure;
plot(temperatures_fahrenheit, m_dot_outlet_lbm_hr, '-o');
xlabel('Inlet Temperature (°F)');
ylabel('Outlet Mass Flow Rate (lbm/hr)');
title('Inlet Temperature (°F) vs Outlet Mass Flow Rate');

figure;
plot(temperatures_fahrenheit, (T_outlet - 273.15) * 9/5 + 32, '-o');
xlabel('Inlet Temperature (°F)');
ylabel('Turbine Outlet Inlet Temperature (°F)');
title('Inlet Temperature (°F) vs Turbine Outlet Temperature (State 5)');

figure;
plot(temperatures_fahrenheit, (T_inlet - 273.15) * 9/5 + 32, '-o');
xlabel('Inlet Temperature (°F)');
ylabel('Turbine Inlet Temperature (°F)');
title('Inlet Temperature (°F) vs Turbine Inlet Temperature (State 48)');

figure;
plot(temperatures_fahrenheit, SFC * 2.20462 * 60, '-o');
xlabel('Inlet Temperature (°F)');
ylabel('Specific Fuel Consumption (lbm/kWhr)');
title('Inlet Temperature (°F) vs Specific Fuel Consumption');

figure;
plot(temperatures_fahrenheit, Heat_rate * 60 * .9478, '-o');
xlabel('Inlet Temperature (°F)');
ylabel('Heat Rate (BTU/kWh)');
title('Inlet Temperature (°F) vs Heat Rate');

figure;
plot(temperatures_fahrenheit, thermal_efficiency, '-o');
xlabel('Inlet Temperature (°F)');
ylabel('Thermal Efficiency');
title('Inlet Temperature (°F) vs Thermal Efficiency');

function T = get_temperature(x, calc_type, P)
    switch calc_type

        case 'h'
            % Temperature (T) as a function of Enthalpy (H)
            T = -2E-12 * x^4 + 3E-08 * x^3 - 0.0001 * x^2 + 1.0021 * x - 15.116;
            
        case 'u'
            % Temperature (T) as a function of Internal Energy (U)
            T = -1E-11 * x^4 + 9E-08 * x^3 - 0.0003 * x^2 + 1.4249 * x - 16.618;
            
        case 's0'
            % Temperature (T) as a function of Entropy (S0)
            T = 0.00048 * x ^ 6.99128;
        
        otherwise
            error('Unknown calculation type. Please specify a valid type.');
    end
end

function [u,h,s0,s] = thermo_properties(P,T)   
   
    % Enthalpy (H) as a function of Temperature (T)
    h = 4E-12*T^4 - 5E-08*T^3 + 0.0002*T^2 + 0.9489*T + 25.787;
    
    % Internal Energy (U) as a function of Temperature (T)
    u = 4E-12 * T^4 - 5E-08 * T^3 + 0.0002 * T^2 + 0.6375 * T + 25.811;
    
    % Entropy (S0) as a function of Temperature (T)
    s0 = (T / 0.00048) ^ (1/6.99128);
    
    % Entropy (S) with pressure correction
    s = s0 - (8.314/28.84)*log(P/99.3);  
end

