classdef TestMomentModels < matlab.unittest.TestCase
    

    methods (TestClassSetup)
        function addSrcToPath(tc)
            % Projekt-Root ausgehend von dieser Testdatei bestimmen:
            testFile = mfilename('fullpath');           % .../tests/TestMomentModels.m
            rootDir  = fileparts(fileparts(testFile));  % .../ (eine Ebene über tests)
            srcDir   = fullfile(rootDir, 'src');

            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                srcDir, 'IncludeSubfolders', true));
        end
    end
    
    methods (Test)
        function testSmokeAllFamilies(testCase)
            cfg = mm_default_config();
            reg = mm_model_registry(cfg);
            psiIso = @(mu) 0.25 + 0.*mu;

            for i = 1:numel(reg)
                model = mm_build_model(reg(i).name, reg(i).order, cfg.models);
                quad = mm_build_quadrature(model, cfg.quad, 'moment');
                u = mm_project_density_to_moments(psiIso, model, quad);
                [f, ~] = mm_eval_flux_function(u, model, quad, cfg.optimizer, struct());
                testCase.verifyFalse(any(~isfinite(f)), sprintf('Non-finite flux for %s-%d', reg(i).name, reg(i).order));
            end
        end

        function testHatPartitionOfUnity(testCase)
            cfg = mm_default_config();
            model = mm_build_model('HFMn', 8, cfg.models);
            mu = linspace(-1, 1, 1001).';
            B = model.basis_eval(mu);
            s = sum(B, 2);
            testCase.verifyLessThan(max(abs(s - 1.0)), 5e-12);
        end

        function testLimiterThetaRange(testCase)
            cfg = mm_default_config();
            model = mm_build_model('HFMn', 4, cfg.models);
            quad = mm_build_quadrature(model, cfg.quad, 'lp');

            uCell = repmat([0.2; 0.3; 0.1; 0.15], 1, 3);
            uL = uCell - 0.4;
            uR = uCell + [0.5; -0.6; 0.4; -0.5];

            [uL2, uR2, st] = mm_apply_realizability_limiter(uCell, uL, uR, model, cfg.limiter, quad);
            testCase.verifyGreaterThanOrEqual(min(st.theta), 0.0);
            testCase.verifyLessThanOrEqual(max(st.theta), 1.0 + 1e-12);

            for i = 1:size(uCell, 2)
                [okL, ~] = mm_is_realizable(uL2(:, i), model, quad, cfg.limiter);
                [okR, ~] = mm_is_realizable(uR2(:, i), model, quad, cfg.limiter);
                testCase.verifyTrue(okL && okR);
            end
        end

        function testPaperCharacteristicLPSmoke(testCase)
            cfg = mm_default_config();
            model = mm_build_model('PN', 3, cfg.models);
            quad = mm_build_quadrature(model, cfg.quad, 'lp');

            nCells = 3;
            uCell = repmat([1.0; 0.0; 0.3; 0.0], 1, nCells);
            uL = uCell + repmat([0.0; 0.3; -0.6; 0.2], 1, nCells);
            uR = uCell + repmat([0.0; -0.25; 0.55; -0.15], 1, nCells);

            lim = cfg.limiter;
            lim.type = 'paper';
            lim.paper_lp_characteristic = true;
            lim.characteristic_basis = struct();
            lim.characteristic_basis.V = repmat({eye(model.nMom)}, 1, nCells);
            lim.characteristic_basis.Vinv = repmat({eye(model.nMom)}, 1, nCells);

            [uL2, uR2, st] = mm_apply_realizability_limiter(uCell, uL, uR, model, lim, quad);

            testCase.verifyTrue(isfield(st, 'theta_components'));
            testCase.verifyEqual(size(st.theta_components), [model.nMom, nCells]);
            testCase.verifyGreaterThanOrEqual(min(st.theta), 0.0);
            testCase.verifyLessThanOrEqual(max(st.theta), 1.0 + 1e-12);

            for i = 1:nCells
                [okL, ~] = mm_is_realizable(uL2(:, i), model, quad, lim);
                [okR, ~] = mm_is_realizable(uR2(:, i), model, quad, lim);
                testCase.verifyTrue(okL && okR);
            end
        end

        function testEntropySolverIsotropic(testCase)
            cfg = mm_default_config();
            model = mm_build_model('PMMn', 4, cfg.models);
            quad = mm_build_quadrature(model, cfg.quad, 'moment');

            psi = @(mu) 0.2 + 0.*mu;
            u = mm_project_density_to_moments(psi, model, quad);
            [alpha, info] = mm_entropy_dual_solve(u, model, quad, cfg.optimizer, struct());

            testCase.verifyTrue(all(isfinite(alpha)));
            testCase.verifyTrue(info.converged);
        end

        function testPaper1TrendAndExactness(testCase)
            cfg = mm_default_config();
            cfg.models.families = {'HFMn', 'PMMn', 'PMPn'};
            cfg.models.orders.HFMn = [4, 8];
            cfg.models.orders.PMMn = [4, 8];
            cfg.models.orders.PMPn = [4, 8];
            res = run_paper1_gauss_heaviside(cfg);
            T = res.table;

            tH = T(strcmp(T.case, 'Gauss') & strcmp(T.model, 'HFMn'), :);
            testCase.verifyLessThan(tH.L1(tH.nMom == 8), tH.L1(tH.nMom == 4));

            tP = T(strcmp(T.case, 'Gauss') & strcmp(T.model, 'PMMn'), :);
            testCase.verifyLessThan(tP.L1(tP.nMom == 8), tP.L1(tP.nMom == 4));

            tExact = T(strcmp(T.case, 'Heaviside') & strcmp(T.model, 'PMPn') & T.nMom == 4, :);
            testCase.verifyLessThan(tExact.L1, 3e-2);
        end

        function testPlaneSourceSymmetryNonnegAndTrend(testCase)
            cfg = mm_default_config();
            cfg.models.families = {'HFMn'};
            cfg.models.orders.HFMn = [4, 8];
            cfg.paper2.n_cells = 40;
            cfg.paper2.tf = 0.25;
            cfg.reference.n_cells = 180;
            cfg.reference.n_mu = 64;
            cfg.reference.tf = cfg.paper2.tf;
            cfg.reference.domain = cfg.paper2.domain;

            res = run_paper2_plane_source(cfg);
            T = res.table;

            testCase.verifyLessThan(max(T.symmetry_L1_rel), 0.3);
            testCase.verifyGreaterThan(min(T.min_rho), -1e-8);

            t4 = T.L1(T.nMom == 4);
            t8 = T.L1(T.nMom == 8);
            testCase.verifyLessThanOrEqual(t8, 1.2 * t4);
        end

        function testMCLLimiterSmoke(testCase)
            cfg = mm_default_config();
            cfg.models.families = {'HFMn'};
            cfg.models.orders.HFMn = [4];
            cfg.paper2.n_cells = 30;
            cfg.paper2.tf = 0.08;
            cfg.limiter.type = 'mcl';

            res = run_paper2_plane_source(cfg);
            T = res.table;

            testCase.verifyTrue(all(isfinite(T.L1)));
            testCase.verifyGreaterThan(min(T.min_rho), -1e-8);
            testCase.verifyLessThan(max(T.symmetry_L1_rel), 0.5);
        end

        function testM1StandaloneMCLSmoke(testCase)
            cfg = mm_default_config();
            cfg.paths.results = tempname;
            mkdir(cfg.paths.results);
            cfg.models.families = {'MN'};
            cfg.models.orders.MN = [1];
            cfg.paper2.n_cells = 24;
            cfg.paper2.tf = 0.08;
            cfg.reference.n_cells = 120;
            cfg.reference.n_mu = 48;
            cfg.reference.tf = cfg.paper2.tf;
            cfg.reference.domain = cfg.paper2.domain;
            cfg.limiter.type = 'mcl';
            cfg.io.close_figures = true;

            res = run_paper2_plane_source(cfg);
            T = res.table;
            testCase.verifyTrue(all(isfinite(T.L1)));
            testCase.verifyGreaterThan(min(T.min_rho), -1e-8);
        end

        function testHoloSyncConstraint(testCase)
            cfg = mm_default_config();
            cfg.reconstruction.use_characteristic = false;
            cfg.holo.sync_each_stage = true;
            cfg.holo.use_mcl_pn = true;
            cfg.holo.use_mcl_m1 = true;
            cfg.quad.base_legendre_order = 20;
            cfg.limiter.mcl_bisect_iter = 16;

            model_pn = mm_build_model('PN', 3, cfg.models);
            model_m1 = mm_build_model('MN', 1, cfg.models);

            nCells = 8;
            z = linspace(-0.35, 0.35, nCells);
            dz = z(2) - z(1);
            grid = struct('z', z(:), 'dz', dz, 'nCells', nCells, 'ghost', 1);

            uIso = model_pn.b_iso * 0.08;
            state = struct();
            state.u_pn = repmat(uIso, 1, nCells);
            state.u_pn(2, :) = linspace(-0.01, 0.01, nCells);
            state.u_m1 = state.u_pn(1:2, :);

            phys = struct();
            phys.sigma_s = 0.0;
            phys.sigma_a = 0.0;
            phys.Q = 0.0;
            phys.psi_vac_density = cfg.physics.psi_vac_density;
            phys.boundary = struct('left', cfg.physics.psi_vac_density, 'right', cfg.physics.psi_vac_density);

            step_cfg = struct();
            step_cfg.optimizer = cfg.optimizer;
            step_cfg.reconstruction = cfg.reconstruction;
            step_cfg.limiter = cfg.limiter;
            step_cfg.quad_cfg = cfg.quad;
            step_cfg.holo = cfg.holo;
            step_cfg.sync = struct('weight', 'identity', 'each_stage', true);

            [state_next, st] = mm_step_holo_m1pn_flux_rk2(state, struct('pn', model_pn, 'm1', model_m1), phys, grid, 0.25 * dz, step_cfg, struct());
            testCase.verifyLessThan(st.sync_stage1.max_constraint_violation, 1e-10);
            testCase.verifyLessThan(st.sync_final.max_constraint_violation, 1e-10);
            testCase.verifyLessThan(max(abs(state_next.u_pn(1:2, :) - state_next.u_m1), [], 'all'), 1e-10);
        end

        function testHoloLimiterThetaRange(testCase)
            cfg = mm_default_config();
            cfg.reconstruction.use_characteristic = false;
            cfg.holo.sync_each_stage = true;
            cfg.holo.use_mcl_pn = true;
            cfg.holo.use_mcl_m1 = true;
            cfg.quad.base_legendre_order = 20;
            cfg.limiter.mcl_bisect_iter = 12;

            model_pn = mm_build_model('PN', 3, cfg.models);
            model_m1 = mm_build_model('MN', 1, cfg.models);

            nCells = 7;
            z = linspace(-0.3, 0.3, nCells);
            dz = z(2) - z(1);
            grid = struct('z', z(:), 'dz', dz, 'nCells', nCells, 'ghost', 1);

            uIso = model_pn.b_iso * 0.07;
            state = struct('u_pn', repmat(uIso, 1, nCells));
            state.u_pn(2, :) = 0.015 * sin(linspace(0, pi, nCells));
            state.u_m1 = state.u_pn(1:2, :);

            phys = struct();
            phys.sigma_s = 0.0;
            phys.sigma_a = 0.0;
            phys.Q = 0.0;
            phys.psi_vac_density = cfg.physics.psi_vac_density;
            phys.boundary = struct('left', cfg.physics.psi_vac_density, 'right', cfg.physics.psi_vac_density);

            step_cfg = struct();
            step_cfg.optimizer = cfg.optimizer;
            step_cfg.reconstruction = cfg.reconstruction;
            step_cfg.limiter = cfg.limiter;
            step_cfg.quad_cfg = cfg.quad;
            step_cfg.holo = cfg.holo;
            step_cfg.sync = struct('weight', 'identity', 'each_stage', true);

            [~, st] = mm_step_holo_m1pn_flux_rk2(state, struct('pn', model_pn, 'm1', model_m1), phys, grid, 0.2 * dz, step_cfg, struct());

            tp = [st.stage1.pn_limiter.theta(:); st.stage2.pn_limiter.theta(:)];
            tm = [st.stage1.m1_limiter.theta(:); st.stage2.m1_limiter.theta(:)];
            testCase.verifyGreaterThanOrEqual(min(tp), 0.0);
            testCase.verifyLessThanOrEqual(max(tp), 1.0 + 1e-12);
            testCase.verifyGreaterThanOrEqual(min(tm), 0.0);
            testCase.verifyLessThanOrEqual(max(tm), 1.0 + 1e-12);
        end

        function testHoloPlaneSourceSymmetry(testCase)
            cfg = mm_default_config();
            cfg.paths.results = tempname;
            mkdir(cfg.paths.results);
            cfg.holo.pn_orders = [3];
            cfg.paper2.n_cells = 24;
            cfg.paper2.tf = 0.12;
            cfg.reference.n_cells = 120;
            cfg.reference.n_mu = 48;
            cfg.reference.tf = cfg.paper2.tf;
            cfg.reference.domain = cfg.paper2.domain;
            cfg.io.close_figures = true;

            res = run_holo_m1pn_plane_source(cfg);
            T = res.table;
            testCase.verifyLessThan(max(T.symmetry_L1_rel), 0.35);
            testCase.verifyGreaterThan(min(T.min_rho), -1e-8);
        end

        function testHoloOutputArtifacts(testCase)
            cfg = mm_default_config();
            cfg.paths.results = tempname;
            mkdir(cfg.paths.results);
            cfg.holo.pn_orders = [3, 7];
            cfg.paper2.n_cells = 16;
            cfg.paper2.tf = 0.05;
            cfg.reference.n_cells = 80;
            cfg.reference.n_mu = 32;
            cfg.reference.tf = cfg.paper2.tf;
            cfg.reference.domain = cfg.paper2.domain;
            cfg.io.close_figures = true;

            res = run_holo_m1pn_plane_source(cfg);
            testCase.verifyTrue(isfile(res.csv));

            testCase.verifyTrue(isfile(fullfile(cfg.paths.results, 'paper2_plane_source_HOLOM1PN_3_profile.csv')));
            testCase.verifyTrue(isfile(fullfile(cfg.paths.results, 'paper2_plane_source_HOLOM1PN_7_profile.csv')));
            testCase.verifyTrue(isfile(fullfile(cfg.paths.results, 'paper2_plane_source_HOLOM1PN_3.png')));
            testCase.verifyTrue(isfile(fullfile(cfg.paths.results, 'paper2_plane_source_HOLOM1PN_7.png')));
        end

        function testCompareHoloVsAll(testCase)
            cfg = mm_default_config();
            cfg.paths.results = tempname;
            mkdir(cfg.paths.results);
            cfg.io.close_figures = true;
            cfg.models.families = {'PN', 'MN'};
            cfg.models.orders.PN = [3, 7];
            cfg.models.orders.MN = [1];
            cfg.holo.pn_orders = [3, 7];
            cfg.holo.compare_include_baselines = true;
            cfg.paper2.n_cells = 18;
            cfg.paper2.tf = 0.06;
            cfg.reference.n_cells = 90;
            cfg.reference.n_mu = 40;
            cfg.reference.tf = cfg.paper2.tf;
            cfg.reference.domain = cfg.paper2.domain;

            cmp = compare_holo_m1pn_vs_all(cfg);
            labels = strcat(string(cmp.table.model), "-", string(cmp.table.order));
            mustHave = ["PN-3", "PN-7", "MN-1", "HOLOM1PN-3", "HOLOM1PN-7"];
            testCase.verifyTrue(all(ismember(mustHave, labels)));
        end

        function testFigureToggleRespected(testCase)
            close all force;

            cfg = mm_default_config();
            cfg.paths.results = tempname;
            mkdir(cfg.paths.results);
            cfg.io.close_figures = true;
            cfg.holo.compare_include_baselines = false;
            cfg.holo.pn_orders = [3];
            cfg.paper2.n_cells = 14;
            cfg.paper2.tf = 0.04;
            cfg.reference.n_cells = 70;
            cfg.reference.n_mu = 32;
            cfg.reference.tf = cfg.paper2.tf;
            cfg.reference.domain = cfg.paper2.domain;

            compare_holo_m1pn_vs_all(cfg);
            figs = findall(0, 'Type', 'figure');
            testCase.verifyEqual(numel(figs), 0);
        end
    end
end
