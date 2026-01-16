% =============================================================================
%            NON-MAGNETIC KIRIGAMI SHEET LAMMPS DATA FILE GENERATOR
% =============================================================================
% Generates LAMMPS data file for a non-magnetic kirigami sheet.
%
% Units: L = 10^-3 m,  T = 10^-3 s,  F = 10^-3 N
%        E = 10^3 N/m^2,  rho = 10^3 kg/m^3,  a = 10^3 m/s^2
% =============================================================================

clear;


% =============================================================================
%                           CONFIGURATION
% =============================================================================

% File paths
inpfile = 'kirigami_sheet_d11_t015_m03.inp';
lamfile = 'kirigami_sheet_d11_t030_m03.lam';

% Unit scaling factors
stress_scale  = 1e3;
density_scale = 1e3;

% Material properties
G   = (1.5/3) * 1e6 / stress_scale;    % Shear modulus
rho = 4e3 / density_scale;             % Density

% Geometric parameters
t = 0.30;                              % Film thickness
a = 0.35;                              % Lattice size (for diameter)

% Derived elastic properties
E  = 3 * G;                            % Young's modulus
nu = 0.4;                              % Poisson's ratio
S  = E * t;                            % Stretching stiffness
B  = (E * t^3) / (12 * (1 - nu^2));    % Bending stiffness


% =============================================================================
%                           READ INP FILE
% =============================================================================

fid = fopen(inpfile);

% Skip header lines
for i = 1:9
    fgetl(fid);
end

% Read node data (only first half of nodes)
data = (fscanf(fid, '%i, %f, %f, %f \n', [4 inf]))';
NP   = length(data);
P    = data(1:NP/2, 2:4);

% Read triangle connectivity
fgetl(fid);
data = (fscanf(fid, '%i, %i, %i, %i, %i, %i, %i\n', [7 inf]))';
T    = data(:, 2:4);

fclose(fid);

% Update counts
NP = size(P, 1);
NT = size(T, 1);


% =============================================================================
%                           CENTER THE MESH
% =============================================================================

Px = (max(P(:,1)) + min(P(:,1))) / 2;
Py = (max(P(:,2)) + min(P(:,2))) / 2;

P(:,1) = P(:,1) - Px;
P(:,2) = P(:,2) - Py;

% Store initial position
p = P;


% =============================================================================
%                       CORRECT TRIANGLE ORIENTATION
% =============================================================================

P1 = P(T(:,1), :);
P2 = P(T(:,2), :);
P3 = P(T(:,3), :);

% Compute signed area (positive = CCW orientation)
Area = P1(:,1).*P2(:,2) - P2(:,1).*P1(:,2) + ...
       P2(:,1).*P3(:,2) - P3(:,1).*P2(:,2) + ...
       P3(:,1).*P1(:,2) - P1(:,1).*P3(:,2);
Area = Area / 2;

% Flip triangles with negative area (CW orientation)
Ineg = find(Area < 0);

if ~isempty(Ineg)
    tmp          = T(Ineg, 1);
    T(Ineg, 1)   = T(Ineg, 2);
    T(Ineg, 2)   = tmp;
    Area(Ineg)   = -Area(Ineg);
end


% =============================================================================
%                           VISUALIZE MESH
% =============================================================================

figure; hold on;
triplot(T, P(:,1), P(:,2), 'r');
triplot(T, P(:,1), P(:,2));
title('Kirigami Sheet Mesh');
xlabel('X'); ylabel('Y');
axis equal;


% =============================================================================
%                       RECOMPUTE AREAS & NODAL MASS
% =============================================================================

P1 = P(T(:,1), :);
P2 = P(T(:,2), :);
P3 = P(T(:,3), :);

Area = P1(:,1).*P2(:,2) - P2(:,1).*P1(:,2) + ...
       P2(:,1).*P3(:,2) - P3(:,1).*P2(:,2) + ...
       P3(:,1).*P1(:,2) - P1(:,1).*P3(:,2);
Area = Area / 2;

% Distribute triangle mass to vertices (1/3 each)
pmass = zeros(NP, 1);
pmass(T(:,1)) = pmass(T(:,1)) + (1/3) * rho * Area * t;
pmass(T(:,2)) = pmass(T(:,2)) + (1/3) * rho * Area * t;
pmass(T(:,3)) = pmass(T(:,3)) + (1/3) * rho * Area * t;

% Clamp small masses to improve timestep stability
apmass  = sum(pmass) / length(pmass);
small_m = pmass < apmass / 3.0;
pmass(small_m) = apmass / 3.0;


% =============================================================================
%                           BUILD TOPOLOGY
% =============================================================================

TR = triangulation(T, P);
E  = edges(TR);
NE = size(E, 1);


% =============================================================================
%                           BOND LIST & COEFFICIENTS
% =============================================================================

Bond     = E;
BondCoef = zeros(NE, 2);
TID      = edgeAttachments(TR, Bond);

for ii = 1:NE
    tid = TID{ii};
    et  = length(tid);
    for jj = 1:et
        BondCoef(ii, 1) = BondCoef(ii, 1) + sqrt(3)/4 * S / et;
    end
end

% Compute rest lengths
pdiff         = P(Bond(:,1), :) - P(Bond(:,2), :);
BondCoef(:,2) = sqrt(sum(pdiff.^2, 2));


% =============================================================================
%                       DIHEDRAL LIST & COEFFICIENTS
% =============================================================================

TID      = edgeAttachments(TR, Bond);
Dihedral = zeros(NE, 4);
kd       = zeros(NE, 1);
cnt      = 0;

for i = 1:NE
    tid = TID{i};
    
    % Only interior edges have two attached triangles
    if length(tid) > 1
        cnt = cnt + 1;
        
        pid1 = T(tid(1), :);
        pid2 = T(tid(2), :);
        pid3 = E(i, :);
        
        [~, loc1] = ismember(pid1, pid3);
        [~, loc2] = ismember(pid2, pid3);
        
        Dihedral(cnt, :) = [pid1(~loc1), pid3(1), pid3(2), pid2(~loc2)];
        
        EdgeLengthSq = (P(pid3(1),1) - P(pid3(2),1))^2 + ...
                       (P(pid3(1),2) - P(pid3(2),2))^2 + ...
                       (P(pid3(1),3) - P(pid3(2),3))^2;
        AreaCombined = Area(tid(1)) + Area(tid(2));
        
        kd(cnt) = B * EdgeLengthSq / AreaCombined;
    end
end

% Trim unused entries
Dihedral(cnt+1:end, :) = [];
kd(cnt+1:end, :)       = [];
DihedralCoef           = [kd, ones(cnt, 1), ones(cnt, 1)];


% =============================================================================
%                           IMPROPER ANGLES
% =============================================================================

Improper         = zeros(NT, 4);
Improper(:, 1)   = NP + 1 : NP + NT;
Improper(:, 2:4) = T;


% =============================================================================
%                       TRIANGLE CENTER ATOMS
% =============================================================================

NT = size(T, 1);

n1 = T(:, 1);
n2 = T(:, 2);
n3 = T(:, 3);

atom_c = zeros(NT, 3);
atom_c(:, 1) = P(n1, 1) + P(n2, 1) + P(n3, 1);
atom_c(:, 2) = P(n1, 2) + P(n2, 2) + P(n3, 2);
atom_c(:, 3) = P(n1, 3) + P(n2, 3) + P(n3, 3);
atom_c = atom_c / 3;

atom_c_t = ones(NT, 1) * 2;

% Center atom masses
cmass   = rho * Area * t;
acmass  = sum(cmass) / NT;
small_c = cmass < acmass / 2.0;
cmass(small_c) = acmass / 2.0;


mu_c = zeros(NT, 3);    % Dipole moments needed to define angles
mu_c(:, 3) = 1.0;       % initially along +z

% =============================================================================
%                           DIPOLE ANGLES
% =============================================================================

Angle      = [];
Angle_type = [];
Angle_g0   = [];
ina        = 0;

for ii = 1:NT
    mu_r = sqrt(mu_c(ii, :) * mu_c(ii, :)');
    
    if mu_r > 1e-10
        r1 = P(T(ii, 1), :) - atom_c(ii, :);
        r2 = P(T(ii, 2), :) - atom_c(ii, :);
        r3 = P(T(ii, 3), :) - atom_c(ii, :);
        
        g1 = acos(mu_c(ii, :) * r1' / sqrt(r1 * r1') / mu_r) / pi * 180;
        g2 = acos(mu_c(ii, :) * r2' / sqrt(r2 * r2') / mu_r) / pi * 180;
        g3 = acos(mu_c(ii, :) * r3' / sqrt(r3 * r3') / mu_r) / pi * 180;
        
        Angle_1 = [NP + ii, T(ii, 1), T(ii, 2)];
        Angle_2 = [NP + ii, T(ii, 2), T(ii, 3)];
        Angle_3 = [NP + ii, T(ii, 3), T(ii, 1)];
        
        Angle      = [Angle; Angle_1; Angle_2; Angle_3];
        Angle_g0   = [Angle_g0; g1; g2; g3];
        Angle_type = [Angle_type; ina + 1; ina + 2; ina + 3];
        
        ina = ina + 3;
    end
end

NB = size(Bond, 1);


% =============================================================================
%                       DIPOLE-DIPOLE BONDS
% =============================================================================

tri_neigh   = neighbors(TR);
dipole_bond = [];

for ii = 1:length(tri_neigh)
    for jj = 1:3
        if ii < tri_neigh(ii, jj)
            new_bond    = [ii + NP, tri_neigh(ii, jj) + NP];
            dipole_bond = [dipole_bond; new_bond];
        end
    end
end

% Dipole bond properties
dipole_bond_id  = (1:length(dipole_bond)) + NB;
dipole_bond_t   = ones(1, length(dipole_bond)) * (NB + 1);
dipole_bond_coe = [0, 1];

% Combine all bonds
bond_t   = [1:NB, dipole_bond_t];
bond_id  = [1:NB, dipole_bond_id];
Bond     = [Bond; dipole_bond];
BondCoef = [BondCoef; dipole_bond_coe];
BondType = 1:NB + 1;


% =============================================================================
%                       COMBINE ALL ATOMS
% =============================================================================

% Molecular IDs: nodes = 1, centers = 2
molID = [ones(NP, 1); 2 * ones(NT, 1)];

% Angle coefficients
AngleCoef = [ones(ina, 1) * B * 20, Angle_g0];

% Atom types and properties
AtomType = [ones(NP, 1); atom_c_t];
atom_xyz = [P; atom_c];
mu       = [zeros(NP, 3); mu_c];

NP    = size(atom_xyz, 1);
tmass = [pmass; cmass];

% Diameter and density
d    = a * ones(NP, 1);
trho = tmass ./ (4/3 * pi * (d/2).^3);

% Charge placeholder
q = zeros(NP, 1);


% =============================================================================
%                       WRITE LAMMPS DATA FILE
% =============================================================================

% Final counts
ND        = size(Dihedral, 1);
NA        = size(Angle, 1);
NI        = size(Improper, 1);
NB        = size(Bond, 1);
NBType    = size(BondCoef, 1);
NAtomType = max(AtomType);

% Box margins
margin = [15, 15; 10, 10; 10, 10];

% --- Open file and write header ---
fid = fopen(lamfile, 'w');

fprintf(fid, 'LAMMPS DATA FILE \n\n');

% Counts
fprintf(fid, '%i atoms \n',     NP);
fprintf(fid, '%i bonds \n',     NB);
fprintf(fid, '%i angles \n',    NA);
fprintf(fid, '%i dihedrals \n', ND);
fprintf(fid, '%i impropers \n', NI);
fprintf(fid, '\n');

% Type counts
fprintf(fid, '%i atom types \n',     NAtomType);
fprintf(fid, '%i bond types \n',     NBType);
fprintf(fid, '%i angle types \n',    ina);
fprintf(fid, '%i dihedral types \n', ND);
fprintf(fid, '%i improper types \n', 1);
fprintf(fid, '\n');

% Box dimensions
fprintf(fid, '%g %g xlo xhi \n', min(p(:,1)) - margin(1,1), max(p(:,1)) + margin(1,2));
fprintf(fid, '%g %g ylo yhi \n', min(p(:,2)) - margin(2,1), max(p(:,2)) + margin(2,2));
fprintf(fid, '%g %g zlo zhi \n', min(p(:,3)) - margin(3,1), max(p(:,3)) + margin(3,2));
fprintf(fid, '\n');

% --- Masses ---
fprintf(fid, 'Masses\n\n');
for i = 1:NAtomType
    fprintf(fid, '%i %g \n', i, 1);
end
fprintf(fid, '\n');

% --- Coefficients ---
fprintf(fid, 'Bond Coeffs \n\n');
fprintf(fid, '%i %.10g %.10g \n', [BondType; BondCoef']);
fprintf(fid, '\n');

fprintf(fid, 'Angle Coeffs \n\n');
fprintf(fid, '%i %.10g %.10g \n', [1:NA; AngleCoef']);
fprintf(fid, '\n');

fprintf(fid, 'Dihedral Coeffs \n\n');
fprintf(fid, '%i %.10g %.10g %.10g \n', [1:ND; DihedralCoef']);
fprintf(fid, '\n');

% --- Topology ---
fprintf(fid, 'Atoms \n\n');
fprintf(fid, '%i %i %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %i\n', ...
        [1:NP; AtomType'; atom_xyz'; d'; trho'; q'; mu'; molID']);
fprintf(fid, '\n');

fprintf(fid, 'Bonds \n\n');
fprintf(fid, '%i %i %i %i \n', [bond_id; bond_t; Bond']);
fprintf(fid, '\n');

fprintf(fid, 'Angles \n\n');
fprintf(fid, '%i %i %i %i %i \n', [1:NA; Angle_type'; Angle']);
fprintf(fid, '\n');

fprintf(fid, 'Dihedrals \n\n');
fprintf(fid, '%i %i %i %i %i %i \n', [1:ND; 1:ND; Dihedral']);
fprintf(fid, '\n');

fprintf(fid, 'Impropers \n\n');
fprintf(fid, '%i %i %i %i %i %i \n', [1:NI; ones(1, NI); Improper']);

fclose(fid);

fprintf('Successfully wrote %s\n', lamfile);