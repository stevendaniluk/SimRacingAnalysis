% This script loads all analysis data for a car from Assetto Corsa Competizione,
% esimates mass and aero properties, and generates plots for the aero data.
%
% This is a sample script intended to demonstrate how to run all parts of the
% analysis. It is currenty configured to analyze the Aston Martin GT3 car, which
% there is sample data for. Note, this car does not have any front splitter
% variability, so comment out any sections for processing or plotting front aero
% element data.

init_paths();

data_dir = 'ACC_data/aston_gt3';
sim = 'ACC';

% Optionally load a pre existing analysis to avoid re running everything
% load(fullfile(data_dir, 'sample_analysis'));

analysis = SimVehicleAnalysis(data_dir, sim);
analysis.loadTrackFactsACC();

%% Process static suspension data to estimate mass and CoG properties
analysis.loadStaticSuspensionData();
analysis.processStaticSuspensionData();

%% Process coast down tests to estimate aero properties
analysis.loadAeroRideHeightData();
analysis.processAeroRideHeightData();

analysis.loadAeroElementData(false);
analysis.processAeroElementData(false);

analysis.loadAeroElementData(true);
analysis.processAeroElementData(true);

%% Optionally save the analysis and a setup with the estimated attributed set
% setup = analysis.base_setup;
% save(fullfile(data_dir, 'sample_setup'), 'setup');
% save(fullfile(data_dir, 'sample_analysis'), 'analysis');

%% Plot all the aero ride height results
analysis.plotRideHeightCoverage();
analysis.plotRideHeightAeroForcePrediction(1:12);

v_const = 180;
v_z_vals = 120:30:240;
rake_min = -0.2;
rake_max = 1.0;
rake_const = 0.4;
RH_f_min = 40;
RH_f_max = 55;
RH_f_const = 50;
RH_f_z_vals = 40:5:55;
analysis.plotAeroPropertiesRideHeight(rake_min, rake_max, rake_const, RH_f_min, RH_f_max, ...
    RH_f_const, RH_f_z_vals, v_const, v_z_vals);

%% Plot all the aero rear wing results
analysis.plotAeroElementForcePrediction(false);

element_min = 1;
element_max = 10;
element_label = 'Rear Wing Setting';
v_z_vals = 120:30:240;
RH_f_const = 50;
rake_const = 0.4;
analysis.plotAeroPropertiesElement(element_min, element_max, element_label, v_z_vals, ...
    RH_f_const, rake_const, false);

%% Plot all the aero front splitter results
analysis.plotAeroElementForcePrediction(true);

element_min = 1;
element_max = 10;
element_label = 'Front Splitter Setting';
v_z_vals = 120:30:240;
RH_f_const = 50;
rake_const = 0.4;
analysis.plotAeroPropertiesElement(element_min, element_max, element_label, v_z_vals, ...
    RH_f_const, rake_const, true);
