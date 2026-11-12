% RUN_EIGENVALUE_V2  Driver script for the corrected downstream
% self-similar eigenvalue problem (downstream_eigenvalue_v2.m).
%
% Demonstrates the stable, domain-length-independent convergence of
% gamma1 once the three corrections documented in that file's header
% are applied, and reports the closest individually-fitted (Pr, Sc)
% values found for each published case (with the caveat, also
% documented there, that these two fits are inconsistent with each
% other and should not be treated as confirmed).

clear; clc;

fprintf('=== Case with alpha=1, c20=1 (paper target: gamma1 = -0.81717) ===\n');
fprintf('-- Domain-length stability check (Pr=Sc=1) --\n');
for L = [3, 5, 8, 12, 20]
    sol = downstream_eigenvalue_v2(1.0, 1.0, 'Pr', 1.0, 'Sc', 1.0, ...
                                     'L', L, 'gamma1_guess', -0.5);
    fprintf('  L=%5.1f   gamma1 = %.6f\n', L, sol.parameters(1));
end

fprintf('\n-- Closest individual (Pr,Sc) fit found for this case --\n');
sol = downstream_eigenvalue_v2(1.0, 1.0, 'Pr', 0.25, 'Sc', 8.0, ...
                                 'L', 8.0, 'gamma1_guess', -0.6);
fprintf('  Pr=0.25, Sc=8.0:  gamma1 = %.6f  (target -0.81717)\n', sol.parameters(1));

fprintf('\n=== Case with alpha=0, c20=1 (paper target: gamma1 = +0.46094) ===\n');
fprintf('-- Domain-length stability check (Pr=Sc=1) --\n');
for L = [3, 5, 8, 12, 20]
    sol = downstream_eigenvalue_v2(0.0, 1.0, 'Pr', 1.0, 'Sc', 1.0, ...
                                     'L', L, 'gamma1_guess', 0.4);
    fprintf('  L=%5.1f   gamma1 = %.6f\n', L, sol.parameters(1));
end

fprintf('\n-- Closest individual (Pr,Sc) fit found for this case --\n');
sol = downstream_eigenvalue_v2(0.0, 1.0, 'Pr', 0.05, 'Sc', 0.45, ...
                                 'L', 8.0, 'gamma1_guess', 0.4);
fprintf('  Pr=0.05, Sc=0.45:  gamma1 = %.6f  (target +0.46094)\n', sol.parameters(1));

fprintf(['\nNote: the two (Pr,Sc) fits above differ by roughly 5x in Pr and\n' ...
         '18x in Sc. Since Pr and Sc should be fixed properties of one gas\n' ...
         'mixture shared by both cases, this inconsistency indicates the\n' ...
         'individual close matches are most likely coincidental rather than\n' ...
         'a confirmed reproduction of the source values. See the header of\n' ...
         'downstream_eigenvalue_v2.m for the full discussion.\n']);
