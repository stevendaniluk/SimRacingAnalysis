% SimVehicleAnalysis
%
% Analyzes MoTeC data produced by a simulator for a collection of scenarios to
% determine the following car parameters:
%   -Sprung mass, and CoG position
%   -Fuel density and tank CoG position
%   -Damper zero positions
%   -Front downforce, rear downforce, and drag of the chassis as a function of
%    velocity, front ride height, and rake angle
%   -Front downforce, rear downforce, and drag of front/rear aero elements
%
% This loads a spreadsheet which has all the car parameters, setup files, and
% motec files for each test performed, and will process the data to estimate the
% above parameters. The final estimated parameters will be saved in the
% base_setup.
%
% This also has plotting utilities to display all the results.
%
% This is capable of supporting any simulator, provided the appropriate
% extensions have been made to the Setup and MotecHandler classes.
classdef SimVehicleAnalysis < handle
    properties
        % Full path to directory with data legend and subdirectories with each test
        data_dir;
        % Base setup which will be populated as the tests are processed
        base_setup;
        % Type of motec handler to use
        motec_handler;

        % Data for the "Stationary Suspension" test
        test_susp;
        % Data for the "Aero Ride Height" test
        test_aero_RH;
        % Data for the "Aero Rear Elements" test
        test_aero_r;
        % Data for the "Aero Front Elements" test
        test_aero_f;
    end

    properties (Constant)
        % Window duration to smooth damper signal over
        dt_smooth_damper = 0.25;

        % Name of the spreadsheet with all the data
        spreadsheet = 'Data Legend.xlsx';
        % Names of the tests in the spreadsheet
        test_susp_name = 'Stationary Suspension';
        test_aero_RH_name = 'Aero Ride Height';
        test_aero_r_name = 'Aero Rear Elements';
        test_aero_f_name = 'Aero Front Elements';
    end

    methods
        % SimVehicleAnalysis
        %
        % INPUTS:
        %   data_dir: Directory with main data spreadsheet and subdirectories
        %             for each test
        %   sim: Name of sim {ACC, iRacing}
        function this = SimVehicleAnalysis(data_dir, sim)
            % Make sure we add the path to the test data so we can find it
            this.data_dir = data_dir;
            addpath(this.data_dir);

            if isequal(sim, 'ACC')
                this.base_setup = SetupAcc();
                this.motec_handler = MotecHandlerAcc();
            elseif isequal(sim, 'iRacing')
                error('iRacing not yet supported');
            else
                error('Unknown sim %s', sim);
            end
        end

        % loadTrackFactsACC
        %
        % Loads the track facts from the spread sheet and populates the baseline
        % setup with the dry mass and fuel consumption data.
        function loadTrackFactsACC(this)
            track_facts_table = readtable(this.spreadsheet, 'Sheet', 'Track Facts');

            this.base_setup.dry_mass_track_map = ...
                containers.Map(track_facts_table.Track, track_facts_table.Mass_kg__fromMoTeCLog_);
            this.base_setup.fuel_consumption_track_map = ...
                containers.Map(track_facts_table.Track, track_facts_table.FuelConsumption_L_s_);
        end

        % loadStaticSuspensionData
        %
        % Loads the setup and motec data defined in the test spreadsheet for
        % static suspension measurements.
        function loadStaticSuspensionData(this)
            disp('Loading static suspension data');
            this.test_susp = this.loadTestData(this.test_susp_name);

            % Override the fuel level in the second setup since it was refuelled
            this.test_susp.setups{2}.V_fuel = this.test_susp.table.Fuel_L_(2);

            % Set the dry mass as defined in the table
            for i=1:this.test_susp.runs
                this.test_susp.setups{i}.m_dry = this.test_susp.table.DryMass_kg_(i);
            end
        end

        % loadAeroRideHeightData
        %
        % Loads the setup and motec data defined in the test spreadsheet for
        % the ride height variations.
        function loadAeroRideHeightData(this)
            disp('Loading aero ride height data data');
            this.test_aero_RH = this.loadTestData(this.test_aero_RH_name);
        end

        % loadAeroElementData
        %
        % Loads the setup and motec data defined in the test spreadsheet for
        % the rear aero elements variations.
        %
        % INPUTS:
        %   is_front: True when for the front aero elements, false for rear
        function loadAeroElementData(this, is_front)
            if is_front
                disp('Loading front aero element data');
                this.test_aero_f = this.loadTestData(this.test_aero_f_name);
            else
                disp('Loading rear aero element data');
                this.test_aero_r = this.loadTestData(this.test_aero_r_name);
            end
        end

        % processStaticSuspensionData
        %
        % Processes the stationary suspension test data to estimate the
        % following:
        %   -Sprung mass
        %   -Sprung mass CoG
        %   -Fuel density
        %   -Fuel tank CoG
        %   -Damper position origins (where spring force is zero)
        function processStaticSuspensionData(this)
            disp('Processing static suspension data');
            v_min = 0;
            v_max = 1.0;
            t_trim = 0.5;

            % Get the damper positions during the first stationary period
            [t_start, t_end] = this.test_susp.motec{1}.detectCoastPeriod(v_min, v_max, t_trim);
            xf_pre_pit = mean(this.test_susp.motec{1}.getDamperPosAvgF(t_start, t_end));
            xr_pre_pit = mean(this.test_susp.motec{1}.getDamperPosAvgR(t_start, t_end));

            % Get the damper positions during the second stationary period after
            % the refuel
            [t_start, t_end] = this.test_susp.motec{2}.detectCoastPeriod(v_min, v_max, t_trim);
            xf_post_pit = mean(this.test_susp.motec{2}.getDamperPosAvgF(t_start, t_end));
            xr_post_pit = mean(this.test_susp.motec{2}.getDamperPosAvgR(t_start, t_end));

            % Compute the mass and suspension properties
            this.base_setup = computeFuelProperties(this.test_susp.setups{1}, ...
                xf_pre_pit, xr_pre_pit, xf_post_pit, xr_post_pit, this.test_susp.setups{2}.V_fuel);
            this.base_setup = computeSprungMassProperties(this.base_setup, xf_pre_pit, xr_pre_pit);

            this.base_setup = this.base_setup.setRideHeightOffset();

            disp('Done');
        end

        % processAeroRideHeightData
        %
        % Processes the aero ride height test data to extract the measured
        % downforce and drag levels, as well as to estimate the following:
        %   -Front downforce as a function of velocity, front and rear ride
        %   height
        %   -Rear downforce as a function of velocity, front and rear ride
        %   height
        %   -Drag as a function of velocity, front and rear ride height
        function processAeroRideHeightData(this)
            disp('Processing aero ride height data');
            this.test_aero_RH = this.extractCoastTestData(this.test_aero_RH);

            disp('Fitting aero force functions to ride height data');

            disp('Fitting front downforce function...');
            [this.base_setup.DF_f, this.base_setup.DF_f_params] = ...
                fitAeroForceFunction3d(this.test_aero_RH.v, this.test_aero_RH.RH_f, ...
                this.test_aero_RH.rake, this.test_aero_RH.DF_f);

            disp('Fitting rear downforce function...');
            [this.base_setup.DF_r, this.base_setup.DF_r_params] = ...
                fitAeroForceFunction3d(this.test_aero_RH.v, this.test_aero_RH.RH_f, ...
                this.test_aero_RH.rake, this.test_aero_RH.DF_r);

            disp('Fitting drag force function...');
            [this.base_setup.drag, this.base_setup.DF_drag_params] = ...
                fitAeroForceFunction3d(this.test_aero_RH.v, this.test_aero_RH.RH_f, ...
                this.test_aero_RH.rake, this.test_aero_RH.drag);

            disp('Done');
        end

        % processAeroElementData
        %
        % Processes the aero rear element test data to extract the measured
        % downforce and drag levels, computes the relative increase in
        % downforce/drag per run, and fits force delta functions.
        %
        % INPUTS:
        %   is_front: True when for the front aero elements, false for rear
        function processAeroElementData(this, is_front)
            if is_front
                disp('Processing front aero element data');
                test = this.test_aero_f;
            else
                disp('Processing rear aero element data');
                test = this.test_aero_r;
            end
            test = this.extractCoastTestData(test);

            % Compute the change in aero force introduced by the additional
            % elements by computing the predicted downforce from on the
            % velocity and ride heights, using the functions fitted in the ride
            % height test, and subtracting that from the measured values.
            test.DF_f_delta = cell(1, test.runs);
            test.DF_r_delta = cell(1, test.runs);
            test.drag_delta = cell(1, test.runs);
            element_setting = cell(1, test.runs);
            for i=1:test.runs
                DF_f_pred = this.base_setup.DF_f(test.v_runs{i}, ...
                    test.RH_f_runs{i}, test.rake_runs{i});
                DF_r_pred = this.base_setup.DF_r(test.v_runs{i}, ...
                    test.RH_f_runs{i}, test.rake_runs{i});
                drag_pred = this.base_setup.drag(test.v_runs{i}, ...
                    test.RH_f_runs{i}, test.rake_runs{i});

                test.DF_f_delta_runs{i} = test.DF_f_runs{i} - DF_f_pred;
                test.DF_r_delta_runs{i} = test.DF_r_runs{i} - DF_r_pred;
                test.drag_delta_runs{i} = test.drag_runs{i} - drag_pred;

                if is_front
                    setting = test.setups{i}.wing_front;
                else
                    setting = test.setups{i}.wing_rear;
                end
                element_setting{i} = setting * ones(1, length(test.v_runs{i}));

            end

            test.DF_f_delta = cell2mat(test.DF_f_delta_runs);
            test.DF_r_delta = cell2mat(test.DF_r_delta_runs);
            test.drag_delta = cell2mat(test.drag_delta_runs);
            element_setting = cell2mat(element_setting);

            % Fit downforce and drag functions to the delta forces
            if is_front
                disp('Fitting aero force functions to front element data');
            else
                disp('Fitting aero force functions to rear element data');
            end

            disp('Fitting front downforce function...');
            [df_f_delta_fit, p_f] = ...
                fitAeroForceFunction2d(test.v, element_setting, test.DF_f_delta);

            disp('Fitting rear downforce function...');
            [df_r_delta_fit, p_r] = ...
                fitAeroForceFunction2d(test.v, element_setting, test.DF_r_delta);

            disp('Fitting drag force function...');
            [drag_delta_fit, p_drag] = ...
                fitAeroForceFunction2d(test.v, element_setting, test.drag_delta);

            if is_front
                this.test_aero_f = test;
                this.base_setup.DF_f_front_element = df_f_delta_fit;
                this.base_setup.DF_f_front_element_params = p_f;
                this.base_setup.DF_r_front_element = df_r_delta_fit;
                this.base_setup.DF_r_front_element_params = p_r;
                this.base_setup.drag_front_element = drag_delta_fit;
                this.base_setup.drag_front_element_params = p_drag;
            else
                this.test_aero_r = test;
                this.base_setup.DF_f_rear_element = df_f_delta_fit;
                this.base_setup.DF_f_rear_element_params = p_f;
                this.base_setup.DF_r_rear_element = df_r_delta_fit;
                this.base_setup.DF_r_rear_element_params = p_r;
                this.base_setup.drag_rear_element = drag_delta_fit;
                this.base_setup.drag_rear_element_params = p_drag;
            end

            disp('Done');
        end

        % plotRideHeightCoverage
        %
        % Generates a scatter plot of all the front and rear ride heights
        % measured across all aero ride height tests. This plot is informative
        % for assessing if enought ride height variation is present in the test
        % data for a sufficient function to be fitted.
        function plotRideHeightCoverage(this)
            % Plot all the ride height levels covered throughout the tests,
            % coloured for speed level to assess that there is sufficient data
            % for the aero force fits
            v_min = min(this.test_aero_RH.v);
            v_max = max(this.test_aero_RH.v);

            figure(1);
            clf;

            subplot(121);
            scatter(this.test_aero_RH.RH_f, this.test_aero_RH.RH_r, [], ...
                this.test_aero_RH.v - v_min);
            grid on;
            xlabel('Front Ride Height [mm]');
            ylabel('Rear Ride Height [mm]');
            title(sprintf(...
                'Ride Heights Covered Over %d Tests (Purple=%.0f km/h, Yellow=%.0f km/h)', ...
                this.test_aero_RH.runs, 3.6 * v_min, 3.6 * v_max));

            subplot(122);
            scatter(this.test_aero_RH.RH_f, this.test_aero_RH.rake, [], ...
                this.test_aero_RH.v - v_min);
            grid on;
            xlabel('Front Ride Height [mm]');
            ylabel('Rake [\circ]');
            title(sprintf(...
                'Rake Angles Covered Over %d Tests (Purple=%.0f km/h, Yellow=%.0f km/h)', ...
                this.test_aero_RH.runs, 3.6 * v_min, 3.6 * v_max));
        end

        % plotRideHeightAeroForcePrediction
        %
        % Plots the measured and predicted front downforce, rear downforce, and
        % drag levels for each run in the aero ride height tests. The predicted
        % levels come from the fitted functions. This should be used to help
        % validate the quality of fit.
        %
        % INPUTS:
        %   runs: Run numbers to plot (optional, defaults to all)
        function plotRideHeightAeroForcePrediction(this, runs)
            if nargin < 2
                runs = 1:this.test_aero_RH.runs;
            end

            % Plot the measured downforce levels for the front and rear as well
            % as the levels predicted by the fitted function for each run.
            cols = 4;
            rows = ceil(length(runs) / cols);

            % Find the min and max force and time values so all the plots can
            % have consistent scales
            axis_limits = [0, -inf, 0, -inf];
            axis_limits(4) = max(max(this.test_aero_RH.DF_f), max(this.test_aero_RH.DF_r));
            for i=1:length(runs)
                run = runs(i);
                t = this.test_aero_RH.t_runs{run};
                axis_limits(2) = max(axis_limits(2), t(end) - t(1));
            end

            f = figure(2);
            f.Name = 'Measured vs. Predicted Downforce';
            % f.NumberTitle = 'off';
            clf;
            for i=1:length(runs)
                run = runs(i);

                DF_f_pred = this.base_setup.DF_f(this.test_aero_RH.v_runs{run}, ...
                    this.test_aero_RH.RH_f_runs{run}, this.test_aero_RH.rake_runs{run});
                DF_r_pred = this.base_setup.DF_r(this.test_aero_RH.v_runs{run}, ...
                    this.test_aero_RH.RH_f_runs{run}, this.test_aero_RH.rake_runs{run});
                drag_pred = this.base_setup.drag(this.test_aero_RH.v_runs{run}, ...
                    this.test_aero_RH.RH_f_runs{run}, this.test_aero_RH.rake_runs{run});
                t = this.test_aero_RH.t_runs{run};
                t = t - t(1);

                subplot(rows, cols, i);
                hold on;
                grid on;
                plot(t, this.test_aero_RH.DF_f_runs{run}, 'b-', 'DisplayName', 'DF_{F,Meas}');
                plot(t, DF_f_pred, 'b--', 'LineWidth', 1.0, 'DisplayName', 'DF_{F,Pred}');
                plot(t, this.test_aero_RH.DF_r_runs{run}, 'g-', 'DisplayName', 'DF_{R,Meas}');
                plot(t, DF_r_pred, 'g--', 'LineWidth', 1.0, 'DisplayName', 'DF_{R,Pred}');
                plot(t, this.test_aero_RH.drag_runs{run}, 'r-', 'DisplayName', 'Drag_{Meas}');
                plot(t, drag_pred, 'r--', 'LineWidth', 1.0, 'DisplayName', 'Drag_{Pred}');
                axis(axis_limits);
                title(sprintf('Run %d', run));
                xlabel('Time [s]');
                ylabel('Aerodynamic Force [N]');
                lgd = legend(gca, 'show');
                lgd.FontSize = 7;
            end
        end

        % plotAeroElementForcePrediction
        %
        % Plots the measured and predicted front downforce, rear downforce, and
        % drag levels for each run in a front/rear aero element test. This will
        % plot two pricted lines: one only from the fitted ride height
        % functions, and another that is the sum of the fitted ride height
        % function output and the fitted front/rear aero element functions.
        % This should be used to help validate the quality of fit.
        %
        % INPUTS:
        %   is_front: True when for the front aero elements, false for rear
        %   runs: Run numbers to plot (optional, defaults to all)
        function plotAeroElementForcePrediction(this, is_front, runs)
            fig = figure(6);
            clf;

            if is_front
                fig.Name = 'Measured vs. Predicted Downforce For Front Aero Elements';
                test = this.test_aero_f;
                df_delta_f = this.base_setup.DF_f_front_element;
                df_delta_r = this.base_setup.DF_r_front_element;
                drag_delta = this.base_setup.drag_front_element;
            else
                fig.Name = 'Measured vs. Predicted Downforce For Rear Aero Elements';
                test = this.test_aero_r;
                df_delta_f = this.base_setup.DF_f_rear_element;
                df_delta_r = this.base_setup.DF_r_rear_element;
                drag_delta = this.base_setup.drag_rear_element;
            end

            if nargin < 3
                runs = 1:test.runs;
            end

            % Plot the measured downforce levels for the front and rear as well
            % as the levels predicted by the fitted function for each run.
            cols = 3;
            rows = ceil(length(runs) / cols);

            % Find the min and max force and time values so all the plots can
            % have consistent scales
            axis_limits = [0, -inf, 0, -inf];
            axis_limits(4) = max(max(test.DF_f), max(test.DF_r));
            for i=1:length(runs)
                run = runs(i);
                t = test.t_runs{run};
                axis_limits(2) = max(axis_limits(2), t(end) - t(1));
            end

            for i=1:length(runs)
                run = runs(i);

                element_setting = test.setups{run}.wing_rear * ...
                    ones(1, length(test.v_runs{run}));

                DF_f_pred_0 = this.base_setup.DF_f(test.v_runs{run}, ...
                    test.RH_f_runs{run}, test.rake_runs{run});
                DF_f_delta_pred = df_delta_f(test.v_runs{run}, element_setting);
                DF_f_pred_total = DF_f_pred_0 + DF_f_delta_pred;

                DF_r_pred_0 = this.base_setup.DF_r(test.v_runs{run}, ...
                    test.RH_f_runs{run}, test.rake_runs{run});
                DF_r_delta_pred = df_delta_r(test.v_runs{run}, element_setting);
                DF_r_pred_total = DF_r_pred_0 + DF_r_delta_pred;

                drag_pred_0 = this.base_setup.drag(test.v_runs{run}, ...
                    test.RH_f_runs{run}, test.rake_runs{run});
                drag_delta_pred = drag_delta(test.v_runs{run}, element_setting);
                drag_pred_total = drag_pred_0 + drag_delta_pred;

                t = test.t_runs{run};
                t = t - t(1);

                subplot(rows, cols, i);
                hold on;
                grid on;
                plot(t, test.DF_f_runs{run}, 'b-', 'DisplayName', 'DF_{F,Meas}');
                plot(t, DF_f_pred_0, 'b-.', 'LineWidth', 1.0, 'DisplayName', 'DF_{F,Pred-0W}');
                plot(t, DF_f_pred_total, 'b--', 'LineWidth', 1.0, 'DisplayName', 'DF_{F,Pred}');

                plot(t, test.DF_r_runs{run}, 'g-', 'DisplayName', 'DF_{R,Meas}');
                plot(t, DF_r_pred_0, 'g-.', 'LineWidth', 1.0, 'DisplayName', 'DF_{R,Pred-0W}');
                plot(t, DF_r_pred_total, 'g--', 'LineWidth', 1.0, 'DisplayName', 'DF_{R,Pred}');

                plot(t, test.drag_runs{run}, 'r-', 'DisplayName', 'Drag_{R,Meas}');
                plot(t, drag_pred_0, 'r-.', 'LineWidth', 1.0, 'DisplayName', 'Drag_{R,Pred-0W}');
                plot(t, drag_pred_total, 'r--', 'LineWidth', 1.0, 'DisplayName', 'Drag_{R,Pred}');

                axis(axis_limits);
                title(sprintf('Run %d', run));
                xlabel('Time [s]');
                ylabel('Aerodynamic Force [N]');
                lgd = legend(gca, 'show');
                lgd.FontSize = 7;
            end
        end

        % plotAeroPropertiesRideHeight
        %
        % Generates the following plots:
        %   -Downforce (total)
        %   -Drag
        %   -Aero Balance
        %   -ClA
        %   -CdA
        %   -L/D
        % For the following variations:
        %   -As a function of rake and velocity
        %   -As a function of front ride height and velocity
        %   -As a function of rake and front ride height
        %
        % INPUTS:
        %   rake_min: Minimum rake angle for plots [deg]
        %   rake_max: Maximum rake angle for plots [deg]
        %   rake_const: Constant rake angle to use when not a variable [deg]
        %   RH_f_min: Minimum front ride height for plots [mm]
        %   RH_f_max: Maximum front ride height for plots [mm]
        %   RH_f_const: Constant front ride height to use when not a variable [mm]
        %   RH_f_z_vals: COnstant front ride height to use for multi line plots [mm]
        %   v_const: Constant velocity to use when not a variable [km/h]
        %   v_z_vals: Constant velocity values to use for multi line plots [km/h]
        function plotAeroPropertiesRideHeight(this, rake_min, rake_max, rake_const, ...
                RH_f_min, RH_f_max, RH_f_const, RH_f_z_vals, v_const, v_z_vals)

            rake_x_vals = linspace(rake_min, rake_max);
            RH_f_x_vals = linspace(RH_f_min, RH_f_max);

            % For each variation we'll need to extract the downforce and drag
            % levels, and feed those to the aero properties plotting utility

            % Properties vs. rake with multiple speed lines, constant front RH
            fig = figure(3);
            % fig.NumberTitle = 'off';
            fig.Name = sprintf(...
                'Aero properties as a function of rake and velocity (RH_f=%.1f mm)', ...
                RH_f_const);
            clf;
            for i=1:length(v_z_vals)
                v_vec = (v_z_vals(i) / 3.6) * ones(1, length(rake_x_vals));
                RH_f_vec = RH_f_const * ones(1, length(rake_x_vals));

                DF_f = this.base_setup.DF_f(v_vec, RH_f_vec, rake_x_vals);
                DF_r = this.base_setup.DF_r(v_vec, RH_f_vec, rake_x_vals);
                drag = this.base_setup.drag(v_vec, RH_f_vec, rake_x_vals);

                fig = this.plotAeroProperties(fig, DF_f, DF_r, drag, v_z_vals(i) / 3.6, ...
                    rake_x_vals, 'Rake [\circ]', sprintf('V=%.0f km/h', v_z_vals(i)));
            end

            % Properties vs. RH_f with multiple speed lines, constant front rake
            fig = figure(4);
            % fig.NumberTitle = 'off';
            fig.Name = sprintf(...
                'Aero properties as a function of front ride height and velocity (rake=%.1f deg)'...
                , rake_const);
            clf;
            for i=1:length(v_z_vals)
                v_vec = (v_z_vals(i) / 3.6) * ones(1, length(RH_f_x_vals));
                rake_vec = rake_const * ones(1, length(RH_f_x_vals));

                DF_f = this.base_setup.DF_f(v_vec, RH_f_x_vals, rake_vec);
                DF_r = this.base_setup.DF_r(v_vec, RH_f_x_vals, rake_vec);
                drag = this.base_setup.drag(v_vec, RH_f_x_vals, rake_vec);

                fig = this.plotAeroProperties(fig, DF_f, DF_r, drag, v_z_vals(i) / 3.6, ...
                    RH_f_x_vals, 'Front Ride Height [mm]', sprintf('V=%.0f km/h', v_z_vals(i)));
            end

            % Properties vs. rake with multiple front ride height lines,
            % constant speed
            fig = figure(5);
            % fig.NumberTitle = 'off';
            fig.Name = sprintf(...
                'Aero properties as a function of rake and front ride height (V=%.0f km/h)', ...
                v_const);
            clf;
            for i=1:length(RH_f_z_vals)
                v_vec = (v_const / 3.6) * ones(1, length(rake_x_vals));
                RH_f_vec = RH_f_z_vals(i) * ones(1, length(rake_x_vals));

                DF_f = this.base_setup.DF_f(v_vec, RH_f_vec, rake_x_vals);
                DF_r = this.base_setup.DF_r(v_vec, RH_f_vec, rake_x_vals);
                drag = this.base_setup.drag(v_vec, RH_f_vec, rake_x_vals);

                fig = this.plotAeroProperties(fig, DF_f, DF_r, drag, v_const / 3.6, ...
                    rake_x_vals, 'Rake [\circ]', sprintf('RH_F=%.1fmm', RH_f_z_vals(i)));
            end
        end

        % plotAeroPropertiesElement
        %
        % Generates the following plots:
        %   -Downforce (total)
        %   -Drag
        %   -Aero Balance
        %   -ClA
        %   -CdA
        %   -L/D
        % As a function of aero element value and velocity
        %
        % INPUTS:
        %   test: Test with fitted aero element functions
        %   element_min: Minimum setting for the aero element
        %   element_max: Maximum setting for the aero element
        %   element_label: Name on the X axis
        %   v_z_vals: Constant velocity values to use for multi line plots [km/h]
        %   RH_f_const: Constant front ride height to use for reference forces
        %   rake_const: Constant rake angle to use for reference forces
        %   is_front: True when for the front aero elements, false for rear
        function plotAeroPropertiesElement(this, element_min, element_max, element_label, ...
                v_z_vals, RH_f_const, rake_const, is_front)
            fig = figure(7);
            clf;
            if is_front
                fig.Name = sprintf(...
                'Change in aero properties due to front elements (RH_f=%.1f mm, rake=%0.1f deg)', ...
                RH_f_const, rake_const);
            else
                fig.Name = sprintf(...
                'Change in aero properties due to rear elements (RH_f=%.1f mm, rake=%0.1f deg)', ...
                RH_f_const, rake_const);
            end

            element_vals = linspace(element_min, element_max);
            RH_f_const_vec = RH_f_const * ones(1, length(element_vals));
            rake_const_vec = rake_const * ones(1, length(element_vals));

            for i=1:length(v_z_vals)
                v_vec = (v_z_vals(i) / 3.6) * ones(1, length(element_vals));

                % Get the delta downforce levels from the aero element as well
                % as the reference levels without the aero element
                if is_front
                    DF_f_delta = this.base_setup.DF_f_front_element(v_vec, element_vals);
                    DF_r_delta = this.base_setup.DF_r_front_element(v_vec, element_vals);
                    drag_delta = this.base_setup.drag_front_element(v_vec, element_vals);
                else
                    DF_f_delta = this.base_setup.DF_f_rear_element(v_vec, element_vals);
                    DF_r_delta = this.base_setup.DF_r_rear_element(v_vec, element_vals);
                    drag_delta = this.base_setup.drag_rear_element(v_vec, element_vals);
                end


                DF_f_ref = this.base_setup.DF_f(v_vec, RH_f_const_vec, rake_const_vec);
                DF_r_ref = this.base_setup.DF_r(v_vec, RH_f_const_vec, rake_const_vec);
                drag_ref = this.base_setup.drag(v_vec, RH_f_const_vec, rake_const_vec);

                fig = this.plotDeltaAeroProperties(fig, DF_f_delta, DF_f_ref, DF_r_delta, ...
                    DF_r_ref, drag_delta, drag_ref, v_z_vals(i) / 3.6, element_vals, ...
                    element_label, sprintf('V=%.0f km/h', v_z_vals(i))) ;
            end
        end

        % extractCoastTestData
        %
        % Identifies a coast down periods for each run in a test and extracts
        % the following:
        %   -Timestamps
        %   -Velocity
        %   -Acceleration
        %   -Damper positions
        %   -Ride heights
        %   -Rake angle
        %   -Downforce (front and rear)
        %   -Drag
        %
        % Data for each individual run will be stored in a cell array.
        %
        % INPUTS:
        %   test: Test structure to add data to, must have motec and setup
        %         fields populated
        % OUTPUTS:
        %   test: test structure with all data fields added
        function test = extractCoastTestData(this, test)
            test.t_runs = cell(1, test.runs);
            test.v_runs = cell(1, test.runs);
            test.a_runs = cell(1, test.runs);
            test.xf_runs = cell(1, test.runs);
            test.xr_runs = cell(1, test.runs);
            test.RH_f_runs = cell(1, test.runs);
            test.RH_r_runs = cell(1, test.runs);
            test.rake_runs = cell(1, test.runs);
            test.DF_f_runs = cell(1, test.runs);
            test.DF_r_runs = cell(1, test.runs);
            test.drag_runs = cell(1, test.runs);

            % Load each setup
            for i=1:test.runs
                % Detect the coast period and get the speed, acceleration, and damper positions
                % during that time
                [t_start, t_end] = test.motec{i}.detectCoastPeriod(5, inf, 0.5);
                [test.v_runs{i}, test.t_runs{i}] = test.motec{i}.getSpeed(t_start, t_end);
                [test.a_runs{i}, ~] = test.motec{i}.getLongitudinalAccel(t_start, t_end);
                test.xf_runs{i} = test.motec{i}.getDamperPosAvgF(t_start, t_end);
                test.xr_runs{i} = test.motec{i}.getDamperPosAvgR(t_start, t_end);

                % Smooth the damper signal
                window_size = find((test.t_runs{i} - t_start) > this.dt_smooth_damper, 1);
                test.xf_runs{i} = smoothdata(test.xf_runs{i}, 'movmean', window_size);
                test.xr_runs{i} = smoothdata(test.xr_runs{i}, 'movmean', window_size);

                % Compute the ride height and rake
                [test.RH_f_runs{i}, test.RH_r_runs{i}] = ...
                    test.setups{i}.rideHeightFromDamperPos(test.xf_runs{i}, test.xr_runs{i});
                test.rake_runs{i} = ...
                    test.setups{i}.rakeFromRideHeights(test.RH_f_runs{i}, test.RH_r_runs{i});

                % Compute the downforce  and drag levels
                [DF, balance] = ...
                    test.setups{i}.downforceFromDamperPos(test.xf_runs{i}, test.xr_runs{i});
                test.DF_f_runs{i} = balance .* DF;
                test.DF_r_runs{i} = (1 - balance) .* DF;
                test.drag_runs{i} = test.setups{i}.totalMass() .* -test.a_runs{i};
            end

            % Concatenate all the results
            test.t = cell2mat(test.t_runs);
            test.v = cell2mat(test.v_runs);
            test.a = cell2mat(test.v_runs);
            test.RH_f = cell2mat(test.RH_f_runs);
            test.RH_r = cell2mat(test.RH_r_runs);
            test.rake = cell2mat(test.rake_runs);
            test.DF_f = cell2mat(test.DF_f_runs);
            test.DF_r = cell2mat(test.DF_r_runs);
            test.drag = cell2mat(test.drag_runs);
        end

    end

    methods(Static, Access = protected)
        % plotAeroProperties
        %
        % Utility to generate a set of 6 subplots for the following:
        %   -Downforce (total)
        %   -Drag
        %   -Aero Balance
        %   -ClA
        %   -CdA
        %   -L/D
        %
        % The independent variable is defined in the inputs.
        %
        % INPUTS:
        %   fig: Figure handle to plot on
        %   DF_f: Front downforce [N]
        %   DF_r: Rear downforce [N]
        %   drag: Drag force [N]
        %   velocity: Velocity for downforce/drag levels
        %   x_vals: Independent variable
        %   x_name: Name to print on x axis
        %   name: Name to print in the legend
        function fig = plotAeroProperties(fig, DF_f, DF_r, drag, velocity, x_vals, x_name, name)
            DF = DF_f + DF_r;
            balance = DF_f ./ DF;
            ClA = Setup().aeroForceCoefficient(velocity, DF);
            CdA = Setup().aeroForceCoefficient(velocity, drag);
            L_D = ClA ./ CdA;

            subplot(321);
            hold on;
            grid on;
            plot(x_vals, DF, 'DisplayName', name);
            xlabel(x_name);
            ylabel('Downforce [N]');
            legend(gca, 'show');

            subplot(323);
            hold on;
            grid on;
            plot(x_vals, drag, 'DisplayName', name);
            xlabel(x_name);
            ylabel('Drag [N]');
            legend(gca, 'show');

            subplot(322);
            hold on;
            grid on;
            plot(x_vals, ClA, 'DisplayName', name);
            xlabel(x_name);
            ylabel('C_LA');
            legend(gca, 'show');

            subplot(324);
            hold on;
            grid on;
            plot(x_vals, CdA, 'DisplayName', name);
            xlabel(x_name);
            ylabel('C_DA');
            legend(gca, 'show');

            subplot(325);
            hold on;
            grid on;
            plot(x_vals, 100 * balance, 'DisplayName', name);
            xlabel(x_name);
            ylabel('Balance (front) [%]');
            legend(gca, 'show');

            subplot(326);
            hold on;
            grid on;
            plot(x_vals, L_D, 'DisplayName', name);
            xlabel(x_name);
            ylabel('L/D');
            legend(gca, 'show');
        end

        % plotDeltaAeroProperties
        %
        % Utility to generate a set of 6 subplots for the following:
        %   -Delta Downforce (total)
        %   -Delta Drag
        %   -Delta Aero Balance
        %   -Delta ClA
        %   -Delta CdA
        %   -Delta L/D
        %
        % The independent variable is defined in the inputs. The delta force
        % values are computed with respect to the reference values provided.
        %
        % INPUTS:
        %   fig: Figure handle to plot on
        %   DF_f_delta: Change in front downforce [N]
        %   DF_f_ref: Reference front downforce [N]
        %   DF_r_delta: Change in rear downforce [N]
        %   DF_r_ref: Reference rear downforce [N]
        %   drag_delta: Change in drag force [N]
        %   drag_ref: Reference drag force [N]
        %   velocity: Velocity for downforce/drag levels
        %   x_vals: Independent variable
        %   x_name: Name to print on x axis
        %   name: Name to print in the legend
        function fig = plotDeltaAeroProperties(fig, DF_f_delta, DF_f_ref, DF_r_delta, DF_r_ref, ...
                drag_delta, drag_ref, velocity, x_vals, x_name, name)
            DF_delta = DF_f_delta + DF_r_delta;
            DF_ref = DF_f_ref + DF_r_ref;

            balance = (DF_f_ref + DF_f_delta) ./ (DF_delta + DF_ref);
            balance_ref = DF_f_ref ./ DF_ref;
            balance_delta = balance - balance_ref;

            ClA = Setup().aeroForceCoefficient(velocity, DF_delta + DF_ref);
            ClA_ref = Setup().aeroForceCoefficient(velocity, DF_ref);
            ClA_delta = ClA - ClA_ref;

            CdA = Setup().aeroForceCoefficient(velocity, drag_delta + drag_ref);
            CdA_ref = Setup().aeroForceCoefficient(velocity, drag_ref);
            CdA_delta = CdA - CdA_ref;

            L_D = ClA ./ CdA;
            L_D_ref = ClA_ref ./ CdA_ref;
            L_D_delta = L_D - L_D_ref;

            subplot(321);
            hold on;
            grid on;
            plot(x_vals, DF_delta, 'DisplayName', name);
            xlabel(x_name);
            ylabel('\Delta Downforce [N]');
            legend(gca, 'show');

            subplot(323);
            hold on;
            grid on;
            plot(x_vals, drag_delta, 'DisplayName', name);
            xlabel(x_name);
            ylabel('\Delta Drag [N]');
            legend(gca, 'show');

            subplot(322);
            hold on;
            grid on;
            plot(x_vals, ClA_delta, 'DisplayName', name);
            xlabel(x_name);
            ylabel('\Delta C_LA');
            legend(gca, 'show');

            subplot(324);
            hold on;
            grid on;
            plot(x_vals, CdA_delta, 'DisplayName', name);
            xlabel(x_name);
            ylabel('\Delta C_DA');
            legend(gca, 'show');

            subplot(325);
            hold on;
            grid on;
            plot(x_vals, 100 * balance_delta, 'DisplayName', name);
            xlabel(x_name);
            ylabel('\Delta Balance (front) [%]');
            legend(gca, 'show');

            subplot(326);
            hold on;
            grid on;
            plot(x_vals, L_D_delta, 'DisplayName', name);
            xlabel(x_name);
            ylabel('\Delta L/D');
            legend(gca, 'show');
        end
    end

    methods (Access = protected)

        % loadTestData
        %
        % Loads the setup and motec data defined in the test spreadsheet for
        % a specific test.
        %
        % INPUTS:
        %   test_name: Name of test to load (i.e. sheet and directory name)
        % OUTPUTS:
        %   test: Structure with loaded motec and setup data
        function test = loadTestData(this, test_name)
            addpath(strcat(this.data_dir, '/', test_name));

            warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
            test.table = ...
                readtable(this.spreadsheet, 'Sheet', test_name);
            test.runs = height(test.table);

            % Load all the setups and motec data
            test.setups = cell(1, test.runs);
            test.motec = cell(1, test.runs);
            for i=1:test.runs
                test.motec{i} = this.motec_handler.loadFromFile(test.table.MoTeCExportFilename{i});
                test.setups{i} = ...
                    this.base_setup.loadFromJSON(test.table.SetupFilename{i});

                % Make sure we compute the offset to get ride height from damper
                % positions, which we can do since we know the mass properties
                test.setups{i} = test.setups{i}.setRideHeightOffset();
            end
        end

    end
end
