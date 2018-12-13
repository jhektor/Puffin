%
% Used for generating *.i moose file 
%
%
%numSn =input('number of sn grain   :');
%numImc=input('number of imc grains :');
numSn=6;      % must be >1
numCu6Sn5=4;  % must be >1
%
ElementSize5='80e-9';  % 5*element size
%
% Generate header
%
fileId=fopen('mooseout.i','w');
fprintf(fileId,'# Auto generated moose input file\n');
fprintf(fileId,'#\n');
fprintf(fileId,['# Number of Cu  grains ',num2str(1),'\n']);
fprintf(fileId,'#   Phase field eta0\n');
fprintf(fileId,'#   concentration field c0\n');
fprintf(fileId,['# Number of Cu6Sn5 grains ',num2str(numCu6Sn5),'\n']);
fprintf(fileId,'#   Phase field(s): eta1');
for i=2:numCu6Sn5
  fprintf(fileId,[', eta',num2str(i),'']);
end
fprintf(fileId,'\n');
fprintf(fileId,'#   Concentration field(s): c1');
for i=2:numCu6Sn5
  fprintf(fileId,[', c',num2str(i),'']);
end
fprintf(fileId,'\n');
fprintf(fileId,['# Number of Sn  grains ',num2str(numSn),'\n']);
fprintf(fileId,['#   Phase field(s): eta',num2str(numCu6Sn5+1)]);
for i=(2:numSn)+numCu6Sn5
  fprintf(fileId,[', eta',num2str(i),'']);
end
fprintf(fileId,'\n');
fprintf(fileId,['#   Concentration field(s): c',num2str(numCu6Sn5+1)]);
for i=(2:numSn)+numCu6Sn5
  fprintf(fileId,[', c',num2str(i),'']);
end
fprintf(fileId,'\n');
fprintf(fileId,'#\n\n');
%
% [Mesh] 
%
fprintf(fileId,'[Mesh]\n');
fprintf(fileId,'  type = GeneratedMesh\n');
fprintf(fileId,'  dim = 3\n');
fprintf(fileId,'  elem_type = HEX8\n');
fprintf(fileId,'#\n');
fprintf(fileId,'# test mesh\n');
fprintf(fileId,'  nx = 15\n');
fprintf(fileId,'  ny = 27\n');
fprintf(fileId,'  nz = 15\n');
fprintf(fileId,'# dimensions in nm\n');
fprintf(fileId,'  xmin = -396\n');
fprintf(fileId,'  xmax =  396\n');
fprintf(fileId,'  ymin =  0\n');
fprintf(fileId,'  ymax =  20\n');
fprintf(fileId,'  zmin = -396\n');
fprintf(fileId,'  zmax =  396\n');
fprintf(fileId,'  displacements = ''disp_x disp_y disp_z''\n');
fprintf(fileId,'[]\n\n');
%
% [Problem]
%
fprintf(fileId,'#[Problem]\n');
fprintf(fileId,'#  solve=false\n');
fprintf(fileId,'#[]\n\n');
%
% [GlobalParams]
%
fprintf(fileId,'[GlobalParams]\n');
fprintf(fileId,'  # CahnHilliard needs the third derivatives\n');
fprintf(fileId,'  derivative_order = 3\n');
fprintf(fileId,'  displacements = ''disp_x disp_y disp_z''\n');
fprintf(fileId,'[]\n\n');
%
% [Variables]
%
fprintf(fileId,'[Variables]\n');
fprintf(fileId,'  [./disp_x]\n');
fprintf(fileId,'    scaling = 1e-2\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./disp_y]\n');
fprintf(fileId,'    scaling = 1e-2\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./disp_z]\n');
fprintf(fileId,'    scaling = 1e-2\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./c]\n');
fprintf(fileId,'    #scaling = 1e2\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  # Chemical potential\n');
fprintf(fileId,'  [./w]\n');
fprintf(fileId,'    #scaling = 1e1\n');
fprintf(fileId,'  [../]\n');
%
% -Phase fields
%
fprintf(fileId,'  # Phase fields\n');
fprintf(fileId,'  # Cu\n');
fprintf(fileId,'  [./eta0]\n');
fprintf(fileId,'    #scaling = 1e1\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  # Cu6Sn5 grains\n');
for i=1:numCu6Sn5
  fprintf(fileId,['  [./eta',num2str(i),']\n']);
  fprintf(fileId,'    #scaling = 1e1\n');
  fprintf(fileId,'  [../]\n');
end
fprintf(fileId,'  # Sn grains\n');
for i=(1:numSn)+numCu6Sn5
  fprintf(fileId,['  [./eta',num2str(i),']\n']);
  fprintf(fileId,'    #scaling = 1e1\n');
  fprintf(fileId,'  [../]\n');
end
%
% -Concentrations
%
fprintf(fileId,'  # Concentrations\n');
fprintf(fileId,'  # Cu\n');
fprintf(fileId,'  [./c0]\n');
fprintf(fileId,'    initial_condition = 0.0\n');
fprintf(fileId,'  [../]\n');
%
fprintf(fileId,'  # Cu6Sn5\n');
for i=1:numCu6Sn5
  fprintf(fileId,['  [./c',num2str(i),']\n']);
  fprintf(fileId,'    initial_condition = 0.4350\n');
  fprintf(fileId,'  [../]\n');
end
fprintf(fileId,'  # Sn\n');
for i=(1:numSn)+numCu6Sn5    
  fprintf(fileId,['  [./c',num2str(i),']\n']);
  fprintf(fileId,'    initial_condition = 1.0\n');
  fprintf(fileId,'  [../]\n');
end
fprintf(fileId,'[]\n\n');
%
% [ICs]
%
fprintf(fileId,'[ICs]\n');
fprintf(fileId,'  [./eta0] # Cu grain\n');
fprintf(fileId,'  [../]\n');
for i=1:numCu6Sn5
  fprintf(fileId,['  [./eta',num2str(i),'] # Cu6Sn5 grain\n']);
  fprintf(fileId,'  [../]\n');
end
for i=(1:numSn)+numCu6Sn5
  fprintf(fileId,['  [./eta',num2str(i),'] # Sn grain\n']);
  fprintf(fileId,'  [../]\n');
end
fprintf(fileId,'  [./c] # Concentration of Sn\n');
fprintf(fileId,'    type = VarDepIC\n');
fprintf(fileId,'    variable = c\n');
fprintf(fileId,'    cis = ''c0');
for i=1:(numCu6Sn5+numSn)
  fprintf(fileId,[' c',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    etas = ''eta0');
for i=1:(numCu6Sn5+numSn)
  fprintf(fileId,[' eta',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./w]\n');
fprintf(fileId,'    type = FunctionIC\n');
fprintf(fileId,'    variable = w\n');
fprintf(fileId,'    function = ''if(y<0,-16.4486,-0.74034)'' # generated\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'[]\n\n');
%
% [BCs]
%
fprintf(fileId,'[BCs]\n');
fprintf(fileId,'  [./Periodic]\n');
fprintf(fileId,'    [./xy]\n');
fprintf(fileId,'      auto_direction = ''x z''\n');
fprintf(fileId,'      variable = ''eta0');
for i=1:(numCu6Sn5+numSn)
  fprintf(fileId,[' eta',num2str(i)]);
end
fprintf(fileId,' w c c0');
for i=1:(numCu6Sn5+numSn)
  fprintf(fileId,[' c',num2str(i)]);
end
fprintf(fileId,' disp_x disp_y disp_z''\n');
fprintf(fileId,'    [../]\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./symmy]\n');
fprintf(fileId,'    type = PresetBC\n');
fprintf(fileId,'    variable = disp_y\n');
fprintf(fileId,'    boundary = ''bottom''\n');
fprintf(fileId,'    value = 0\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./symmx]\n');
fprintf(fileId,'    type = PresetBC\n');
fprintf(fileId,'    variable = disp_x\n');
fprintf(fileId,'    boundary = ''bottom''\n');
fprintf(fileId,'    value = 0\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./symmz]\n');
fprintf(fileId,'    type = PresetBC\n');
fprintf(fileId,'    variable = disp_z\n');
fprintf(fileId,'    boundary = ''bottom''\n');
fprintf(fileId,'    value = 0\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'[]\n\n');
%
% [Functions]
%
fprintf(fileId,'[Functions]\n');
fprintf(fileId,'  [./tanh_arg]\n');
fprintf(fileId,'    type = ParsedFunction\n');
fprintf(fileId,'    vars = ''t_start_cool t_end_cool k_tan''\n');
fprintf(fileId,'    vals = ''1 21 0.333333''\n');
fprintf(fileId,'    #value = ''k_tan*(t - (t_start_cool + t_end_cool)/2)''\n');
fprintf(fileId,'    value = ''1000000''\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'[]\n\n');
%
% [UserObjects]
%
fprintf(fileId,'[UserObjects]\n');
fprintf(fileId,'  #Crystal plasticity BCT\n');
for i=(1:numSn)+numCu6Sn5
  fprintf(fileId,['  [./slip_rate_gss',num2str(i),']  # Sn-eta',num2str(i),'\n']);
  fprintf(fileId,'    type = CrystalPlasticitySlipRateGSSBaseName\n');
  fprintf(fileId,'    variable_size = 32\n');
  fprintf(fileId,'    slip_sys_file_name = slip_systems_bct.txt\n');
  fprintf(fileId,'    num_slip_sys_flowrate_props = 2\n');
  fprintf(fileId,'    flowprops = ''1 32 0.001 0.166666667'' #start_ss end_ss gamma0 1/m\n');
  fprintf(fileId,['    uo_state_var_name = state_var_gss',num2str(i),'\n']);
  fprintf(fileId,['    base_name = ''eta',num2str(i),'''\n']);
  fprintf(fileId,'  [../]\n');
  fprintf(fileId,['  [./slip_resistance_gss',num2str(i),']\n']);
  fprintf(fileId,'    type = CrystalPlasticitySlipResistanceGSS\n');
  fprintf(fileId,'    variable_size = 32\n');
  fprintf(fileId,['    uo_state_var_name = state_var_gss',num2str(i),'\n']);
  fprintf(fileId,'  [../]\n');
  fprintf(fileId,['  [./state_var_gss',num2str(i),']\n']);
  fprintf(fileId,'    type = CrystalPlasticityStateVariable\n');
  fprintf(fileId,'    variable_size = 32\n');
  %    #groups = '0 32'
  %    #group_values = '0.144' # 23 MPa in eV/nm^3 initial values of slip resistance
  fprintf(fileId,'    groups = ''0 2 4 6 10 12 16 18 20 24 32''\n');
  fprintf(fileId,'    group_values = ''0.05306 0.02684 0.06492 0.02809 0.03496 0.03184 0.04619 0.09363 0.04120 0.07491''\n');
  fprintf(fileId,['    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss',num2str(i),'\n']);
  fprintf(fileId,'    scale_factor = 1.0\n');
  fprintf(fileId,'  [../]\n');
  fprintf(fileId,['  [./state_var_evol_rate_comp_gss',num2str(i),']\n']);
  fprintf(fileId,'    type = CrystalPlasticityStateVarRateComponentVoceP\n');
  fprintf(fileId,'    variable_size = 32\n');
  fprintf(fileId,'    groups = ''0 2 4 6 10 20 24 32''\n');
  fprintf(fileId,'    h0_group_values = ''0.12484 0.12484 0.12484 0.12484 0.12484 0.12484 0.12484''\n');
  fprintf(fileId,'    tau0_group_values = ''0.0 0.0 0.0 0.0 0.0 0.0 0.0''\n');
  fprintf(fileId,'    tauSat_group_values = ''0.06866 0.05618 0.06866 0.05618 0.06242 0.05618 0.08115''\n');
  fprintf(fileId,'    hardeningExponent_group_values = ''2.0 2.0 2.0 2.0 2.0 2.0 2.0''\n');
  fprintf(fileId,'    coplanarHardening_group_values = ''1.0 1.0 1.0 1.0 1.0 1.0 1.0'' #q_aa = 1\n');
  fprintf(fileId,'    selfHardening_group_values = ''1.4 1.4 1.4 1.4 1.4 1.4 1.4''\n');
  fprintf(fileId,'    crystal_lattice_type = BCT # default is BCT \n');
  fprintf(fileId,['    uo_slip_rate_name = slip_rate_gss',num2str(i),'\n']);
  fprintf(fileId,['    uo_state_var_name = state_var_gss',num2str(i),'\n']);
  fprintf(fileId,'  [../]\n');
end
fprintf(fileId,'[]\n\n');
%
% [Materials]
%
fprintf(fileId,'[Materials]\n');
fprintf(fileId,'  [./time]\n');
fprintf(fileId,'    type = TimeStepMaterial\n');
fprintf(fileId,'    prop_time = time\n');
fprintf(fileId,'    prop_dt = dt\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
%    # output_properties = 'dt'
%    # outputs = exodus
fprintf(fileId,'  [../]\n\n');
% -Cu
fprintf(fileId,'  # Cu\n');
fprintf(fileId,'  [./elasticity_tensor_cu]\n');
fprintf(fileId,'    type = ComputeIsotropicElasticityTensor\n');
fprintf(fileId,'    youngs_modulus = 936. #150 GPa in eV/nm^3\n');
fprintf(fileId,'    poissons_ratio = 0.35\n');
fprintf(fileId,'    base_name = eta0\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./strain_cu]\n');
fprintf(fileId,'    type = ComputeFiniteStrain\n');
fprintf(fileId,'    base_name = eta0\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./stress_cu]\n');
fprintf(fileId,'    type = ComputeFiniteStrainElasticStress\n');
fprintf(fileId,'    base_name = eta0\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./fel_cu]\n');
fprintf(fileId,'    type = ElasticEnergyMaterialGreenPK2\n');
fprintf(fileId,'    args = '' ''\n');
fprintf(fileId,'    base_name = eta0\n');
fprintf(fileId,'    f_name = fe0\n');
fprintf(fileId,'    outputs = exodus\n');
fprintf(fileId,'    output_properties = fe0\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n\n');
% -Cu6Sn5
fprintf(fileId,'  #  Cu6Sn5-material\n');
for i=1:numCu6Sn5
  fprintf(fileId,['  [./elasticity_tensor_eta',num2str(i),']\n']);
  fprintf(fileId,'    type = ComputeIsotropicElasticityTensor\n');
  fprintf(fileId,'    youngs_modulus = 701. #112.3 GPa in eV/nm^3\n');
  fprintf(fileId,'    poissons_ratio = 0.31\n');
  fprintf(fileId,['    base_name = eta',num2str(i),'\n']);
  fprintf(fileId,'  [../]\n');
  fprintf(fileId,['  [./strain_eta',num2str(i),']\n']);
  fprintf(fileId,'    type = ComputeFiniteStrain\n');
  fprintf(fileId,['    base_name = eta',num2str(i),'\n']);
  fprintf(fileId,'    eigenstrain_names = eT_eta\n');
  fprintf(fileId,'  [../]\n');
  fprintf(fileId,['  [./stress_eta',num2str(i),']\n']);
  fprintf(fileId,'    type = ComputeFiniteStrainElasticStress\n');
  fprintf(fileId,['    base_name = eta',num2str(i),'\n']);
  fprintf(fileId,'  [../]\n');
  fprintf(fileId,['  [./eigenstrain_eta',num2str(i),']\n']);
  fprintf(fileId,'    type = ComputeVariableEigenstrain\n');
  fprintf(fileId,'    args = ''eta0');
  for j=1:(numCu6Sn5+numSn)
    fprintf(fileId,[' eta',num2str(j)]);
  end
  fprintf(fileId,'''\n');
  fprintf(fileId,['    base_name = eta',num2str(i),'\n']);
  fprintf(fileId,'    eigen_base = ''1 1 1 0 0 0''\n');
  fprintf(fileId,'    eigenstrain_name = eT_eta\n');
  fprintf(fileId,'    prefactor = pre\n');
  fprintf(fileId,'  [../]\n');
  fprintf(fileId,['  [./fe',num2str(i),']\n']);
  fprintf(fileId,'    type = ElasticEnergyMaterialGreenPK2\n');
  fprintf(fileId,'    args = ''eta0');
  for j=1:(numCu6Sn5+numSn)
    fprintf(fileId,[' eta',num2str(j)]);
  end
  fprintf(fileId,'''\n');
  fprintf(fileId,['    f_name = fe',num2str(i),'\n']);
  fprintf(fileId,['    base_name = eta',num2str(i),'\n']);
  fprintf(fileId,'    eigenstrain = true\n');
  fprintf(fileId,'    eigenstrain_name = eT_eta\n');
  fprintf(fileId,'    use_displaced_mesh = true\n');
  fprintf(fileId,'    outputs = exodus\n');
  fprintf(fileId,['    output_properties = fe',num2str(i),'\n']);
  fprintf(fileId,'  [../]\n');
end
fprintf(fileId,'  [./pre]\n');
fprintf(fileId,'    type = DerivativeParsedMaterial\n');
fprintf(fileId,'    args = ''');    
for i=0:(numSn+numCu6Sn5)
for j=(i+1):(numSn+numCu6Sn5)
    if (i+j)>1 fprintf(fileId,' '); end
    fprintf(fileId,['h',num2str(i),'h',num2str(j)]); 
end
end
fprintf(fileId,'''\n');
fprintf(fileId,'    # first part is Cu-Cu6Sn5 part second is Cu6Sn5-Sn part\n');
fprintf(fileId,'    function = '' 0 ''\n');
fprintf(fileId,'    #function = ''-0.003*(');
i=0;
for j=(i+1):numCu6Sn5
    if j>1 fprintf(fileId,'+'); end
    fprintf(fileId,['h',num2str(i),'h',num2str(j)]);
end
fprintf(fileId,')\n');
fprintf(fileId,'    #            -0.003*(');
for i=1:numCu6Sn5
for j=(numCu6Sn5+1):(numSn+numCu6Sn5)
    if (i+j)>numCu6Sn5+2 fprintf(fileId,'+'); end
    fprintf(fileId,['h',num2str(i),'h',num2str(j)]);
end
end
fprintf(fileId,')''\n');
fprintf(fileId,'    f_name = pre\n');
fprintf(fileId,'    outputs = exodus\n');
fprintf(fileId,'  [../]\n\n');
%
% -Sn
%
angle=0;
fprintf(fileId,'  #  Sn-material\n');
for i=(1:numSn)+numCu6Sn5
  fprintf(fileId,['  [./crysp',num2str(i),']\n']);
  fprintf(fileId,'    type = FiniteStrainUObasedCPBaseName\n');
  fprintf(fileId,'    rtol = 1e-6\n');
  fprintf(fileId,'    abs_tol = 1e-6\n');
  fprintf(fileId,'    stol = 1e-2\n');
  fprintf(fileId,['    uo_slip_rates = ''slip_rate_gss',num2str(i),'''\n']);
  fprintf(fileId,['    uo_slip_resistances = ''slip_resistance_gss',num2str(i),'''\n']);
  fprintf(fileId,['    uo_state_vars = ''state_var_gss',num2str(i),'''\n']);
  fprintf(fileId,['    uo_state_var_evol_rate_comps = ''state_var_evol_rate_comp_gss',num2str(i),'''\n']);
  fprintf(fileId,['    base_name = ''eta',num2str(i),'''\n']);
  fprintf(fileId,'    maximum_substep_iteration = 10\n');
  fprintf(fileId,'    tan_mod_type = exact\n');
  fprintf(fileId,'    # output_properties = ''slip_rate_gss2''\n');
  fprintf(fileId,'    # outputs = exodus\n');
  fprintf(fileId,'  [../]\n');
  fprintf(fileId,['  [./strain',num2str(i),']\n']);
  fprintf(fileId,'    type = ComputeFiniteStrain\n');
  fprintf(fileId,'    #displacements = ''disp_x disp_y disp_z''\n');
  fprintf(fileId,['    base_name = ''eta',num2str(i),'''\n']);
  fprintf(fileId,'  [../]\n');
  fprintf(fileId,['  [./elasticity_tensor',num2str(i),']\n']);
  fprintf(fileId,'    type = ComputeElasticityTensorCPBaseName #Allows for changes due to crystal re-orientation\n');
  fprintf(fileId,'    #C_ijkl = ''72.3e3 59.4e3 35.8e3 72.3e3 35.8e3 88.4e3 22e3 22e3 24e3'' #MPa #From Darbandi 2014 table I\n');
  fprintf(fileId,'    C_ijkl = ''451.26 370.75 223.45 451.26 223.45 551.75 137.31 137.31 149.8022'' #eV/nm^3 #From Darbandi 2013 table I\n');
  fprintf(fileId,'    fill_method = symmetric9\n');
  fprintf(fileId,['    euler_angle_1 = ',num2str(angle),'\n']);
  fprintf(fileId,['    euler_angle_2 = ',num2str(angle+10),'\n']);
  fprintf(fileId,['    euler_angle_3 = ',num2str(angle+20),'\n']);
  fprintf(fileId,['    base_name = ''eta',num2str(i),'''\n']);
    angle=angle+50;
  fprintf(fileId,'  [../]\n');
  fprintf(fileId,['  [./fe',num2str(i),']\n']);
  fprintf(fileId,'    type = ElasticEnergyMaterialGreenPK2\n');
  fprintf(fileId,'    args = '' ''\n');
  fprintf(fileId,['    f_name = fe',num2str(i),'\n']);
  fprintf(fileId,['    base_name = eta',num2str(i),'\n']);
  fprintf(fileId,'    use_displaced_mesh = true\n');
  fprintf(fileId,'    plasticity = true\n');
  fprintf(fileId,'    outputs = exodus\n');
  fprintf(fileId,['    output_properties = fe',num2str(i),'\n']);
  fprintf(fileId,'  [../]\n');
  fprintf(fileId,['  [./fp',num2str(i),']\n']);
  fprintf(fileId,'    type = CPPlasticEnergyMaterial\n');
  fprintf(fileId,'    # Q = 1 in class\n');
  fprintf(fileId,'    variable_size = 32\n');
  fprintf(fileId,['    uo_state_var_name = state_var_gss',num2str(i),'\n']);
  fprintf(fileId,['    f_name = fp',num2str(i),'\n']);
  fprintf(fileId,'    use_displaced_mesh = true\n');
  fprintf(fileId,'    args = '' ''\n');
  fprintf(fileId,'    #s0 = 0.144\n');
  fprintf(fileId,'    groups = ''0 2 4 6 10 12 16 18 20 24 32'' #values calibrated on 110\n');
  fprintf(fileId,'    group_values = ''0.05122774 0.03452174 0.03696807 0.00912421 0.02046358 0.01612225 0.04525029 0.08612754 0.29181706 0.02277457''\n');
  fprintf(fileId,'    # outputs = exodus\n');
  fprintf(fileId,'    # output_properties = fp2\n');
  fprintf(fileId,'  [../]\n\n');
end
%
%
fprintf(fileId,'  [./scale]\n');
fprintf(fileId,'    type = GenericConstantMaterial\n');
fprintf(fileId,'    prop_names = ''length_scale energy_scale time_scale''\n');
fprintf(fileId,'    prop_values = ''1e9 6.24150943e18 1.'' #m to nm J to eV s to h\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./model_constants]\n');
fprintf(fileId,'    type = GenericConstantMaterial # delta is diff reg. width\n');
fprintf(fileId,'    prop_names = ''sigma delta delta_real gamma tgrad_corr_mult''\n');
fprintf(fileId,'    prop_values = ''0.5 ');
fprintf(fileId,ElementSize5);
fprintf(fileId,' 5e-10 1.5 0'' #J/m^2 m - ? \n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./kappa]\n');
fprintf(fileId,'    type = ParsedMaterial\n');
fprintf(fileId,'    material_property_names = ''sigma delta length_scale energy_scale''\n');
fprintf(fileId,'    f_name = kappa\n');
fprintf(fileId,'    function = ''0.75*sigma*delta*energy_scale/length_scale'' #eV/nm\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./mu]\n');
fprintf(fileId,'    type = ParsedMaterial\n');
fprintf(fileId,'    material_property_names = ''sigma delta length_scale energy_scale''\n');
fprintf(fileId,'    f_name = mu\n');
fprintf(fileId,'    function = ''6*(sigma/delta)*energy_scale/length_scale^3'' #eV/nm^3\n');
fprintf(fileId,'    output_properties = mu\n');
fprintf(fileId,'  [../]\n\n');
%
% Switch and Barrier functions
%
for i=0:(numSn+numCu6Sn5)
  fprintf(fileId,['  [./h',num2str(i),']\n']);
  fprintf(fileId,'      type = SwitchingFunctionMultiPhaseMaterial\n');
  fprintf(fileId,['      h_name = h',num2str(i),'\n']);
  fprintf(fileId,'      all_etas = ''eta0');
  for j=1:(numCu6Sn5+numSn)
    fprintf(fileId,[' eta',num2str(j)]);
  end
  fprintf(fileId,'''\n');
  fprintf(fileId,['      phase_etas = eta',num2str(i),'\n']);
  fprintf(fileId,'      use_displaced_mesh = true\n');
  fprintf(fileId,'      outputs = exodus\n');
  fprintf(fileId,['      output_properties = h',num2str(i),'\n']);
  fprintf(fileId,'  [../]\n');
end
for i=0:(numSn+numCu6Sn5)
  fprintf(fileId,['  [./g',num2str(i),']\n']);
  fprintf(fileId,'      type = BarrierFunctionMaterial\n');
  fprintf(fileId,['      eta = eta',num2str(i),'\n']);
  fprintf(fileId,'      well_only = true\n');
  fprintf(fileId,['      function_name = g',num2str(i),'\n']);
  fprintf(fileId,'      g_order = SIMPLE\n');
  fprintf(fileId,'      use_displaced_mesh = true\n');
  fprintf(fileId,'  [../]\n');
end
fprintf(fileId,'\n');
%
%
%
fprintf(fileId,'  [./ACMobility]\n'); 
fprintf(fileId,'      type = GenericConstantMaterial\n');
fprintf(fileId,'      prop_names = L\n');
fprintf(fileId,'      prop_values = 2.7 #2.7\n');
fprintf(fileId,'  [../]\n\n');
%
% -noise 
%
fprintf(fileId,'  [./noise_constants]\n');
fprintf(fileId,'    type = GenericConstantMaterial\n');
fprintf(fileId,'    prop_names = ''T kb lambda dim'' #eta2erature Boltzmann gridsize dimensionality\n');
fprintf(fileId,'    #prop_values = ''493 8.6173303e-5  3''\n');
fprintf(fileId,'    #prop_values = ''493 8.6173303e-5 8 3''\n');
fprintf(fileId,'    prop_values = ''493 8.6173303e-5 16 3''\n');
fprintf(fileId,'  [../]\n');
%
% - nuc
%
fprintf(fileId,'  [./nuc]\n');
fprintf(fileId,'    type =  DerivativeParsedMaterial\n');
fprintf(fileId,'    args = ''eta0');
  for j=1:(numCu6Sn5+numSn)
    fprintf(fileId,[' eta',num2str(j)]);
  end
fprintf(fileId,'''\n');
fprintf(fileId,'    f_name = nuc\n');
fprintf(fileId,'    material_property_names = ''time dt T kb lambda dim L h0 h2 h3 h4''\n');
fprintf(fileId,'    #function = ''if(time<1&(h0*h2>0.09 | h0*h3>0.09),sqrt(2*kb*T*L/(lambda^dim*dt)),0)'' \n');
fprintf(fileId,'    function = '' 1 '' \n');
fprintf(fileId,'    #expression from Shen (2007) not sure about diff reg coeff hi*hj>X?\n');
fprintf(fileId,'    outputs = exodus\n');
fprintf(fileId,'    output_properties = nuc\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n\n');
%
% -energy paramaters
%
fprintf(fileId,'# Constants The energy parameters are for 220 C and for 25 C \n');
fprintf(fileId,'  [./energy_constants_A]\n');
fprintf(fileId,'    type =  GenericFunctionMaterial # GenericConstantMaterial\n');
fprintf(fileId,'    prop_names = ''A_cu_220  A_eps_220   A_eta_220   A_sn_220\n');
fprintf(fileId,'                  A_cu_25   A_eps_25    A_eta_25    A_sn_25''\n');
fprintf(fileId,'    #prop_values = ''1.0133e5/Vm 4e5/Vm 4.2059e6/Vm'' #J/m^3\n');
fprintf(fileId,'    prop_values = ''1.7756e10  2.4555e11   2.4555e11   2.3033e10\n');
fprintf(fileId,'                   6.2204e9   2.4555e10   2.4555e10   2.5819e11'' #J/m^3 Aeps=Aeta=2e6 prev. used; Vm := molar vol = 16.29 [cm^3 mol^-1]\n');
fprintf(fileId,'    #prop_values = ''1.5929e10 2.4555e12 2.4555e12 2.3020e10'' #J/m^3 Aeps = 2e7 Aeta = 2e7\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_constants_B]\n');
fprintf(fileId,'    type = GenericConstantMaterial\n');
fprintf(fileId,'    prop_names = ''B_cu_220  B_eps_220 B_eta_220 B_sn_220\n');
fprintf(fileId,'                  B_cu_25   B_eps_25  B_eta_25  B_sn_25''\n');
fprintf(fileId,'    #prop_values = ''-2.1146e4/Vm -6.9892e3/Vm 7.168e3/Vm'' #J/m^3\n');
fprintf(fileId,'    prop_values = ''-2.6351e9  -1.4014e9   2.3251e7    2.14216e8\n');
fprintf(fileId,'                   -1.2981e9  -4.2905e8   -4.2905e8   4.4002e8'' #J/m^3\n');
fprintf(fileId,'    #prop_values = ''-2.5789e9 -1.3733e9 2.3175e7 2.1406e8'' #J/m^3\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_C]\n');
fprintf(fileId,'    type = GenericConstantMaterial\n');
fprintf(fileId,'    prop_names = ''C_cu_220  C_eps_220  C_eta_220  C_sn_220\n');
fprintf(fileId,'                  C_cu_25   C_eps_25   C_eta_25   C_sn_25''\n');
fprintf(fileId,'    #prop_values = ''-1.2842e4/Vm -1.9185e4/Vm -1.5265e4/Vm'' #J/m^3\n');
fprintf(fileId,'    prop_values = ''-1.1441e9  -1.7294e9   -1.7646e9   -1.646e9\n');
fprintf(fileId,'                   -7.8834e8  -1.1778e9   -1.1778e9   -9.3708e8'' #J/m^3\n');
fprintf(fileId,'    #prop_values = ''-1.1529e9 -1.7330e9 -1.7646e9 -1.646e9'' #J/m^3\n');
fprintf(fileId,'  [../]\n\n');
%
fprintf(fileId,'  [./tanh_arg]\n');
fprintf(fileId,'    type = GenericFunctionMaterial\n');
fprintf(fileId,'    #times for the ramping of mat. constants\n');
fprintf(fileId,'    prop_names = ''phi_t''\n');
fprintf(fileId,'    prop_values = ''tanh_arg'' # [s]\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  # GIBBS ENERGY CONSTANTS ABC\n');
fprintf(fileId,'  [./energy_A_cu]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''A_cu_220 A_cu_25 phi_t''\n');
fprintf(fileId,'      function = ''A_cu_25 - (A_cu_25 - A_cu_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = A_cu\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_A_eps]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''A_eps_220 A_eps_25 phi_t''\n');
fprintf(fileId,'      function = ''A_eps_25 - (A_eps_25 - A_eps_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = A_eps\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_A_eta]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''A_eta_220 A_eta_25 phi_t''\n');
fprintf(fileId,'      function = ''A_eta_25 - (A_eta_25 - A_eta_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = A_eta\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_A_sn]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''A_sn_220 A_sn_25 phi_t''\n');
fprintf(fileId,'      function = ''A_sn_25 - (A_sn_25 - A_sn_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = A_sn\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_B_cu]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''B_cu_220 B_cu_25 phi_t''\n');
fprintf(fileId,'      function = ''B_cu_25 - (B_cu_25 - B_cu_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = B_cu\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_B_eps]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''B_eps_220 B_eps_25 phi_t''\n');
fprintf(fileId,'      function = ''B_eps_25 - (B_eps_25 - B_eps_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = B_eps\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_B_eta]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''B_eta_220 B_eta_25 phi_t''\n');
fprintf(fileId,'      function = ''B_eta_25 - (B_eta_25 - B_eta_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = B_eta\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_B_sn]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''B_sn_220 B_sn_25 phi_t''\n');
fprintf(fileId,'      function = ''B_sn_25 - (B_sn_25 - B_sn_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = B_sn\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_C_cu]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''C_cu_220 C_cu_25 phi_t''\n');
fprintf(fileId,'      function = ''C_cu_25 - (C_cu_25 - C_cu_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = C_cu\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_C_eps]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''C_eps_220 C_eps_25 phi_t''\n');
fprintf(fileId,'      function = ''C_eps_25 - (C_eps_25 - C_eps_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = C_eps\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_C_eta]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''C_eta_220 C_eta_25 phi_t''\n');
fprintf(fileId,'      function = ''C_eta_25 - (C_eta_25 - C_eta_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = C_eta\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_C_sn]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''C_sn_220 C_sn_25 phi_t''\n');
fprintf(fileId,'      function = ''C_sn_25 - (C_sn_25 - C_sn_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = C_sn\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  #END OF GIBBS ENERGY PARAMETER CONSTANTS ABC\n\n');
fprintf(fileId,'  # [./energy_c_ab]\n');
fprintf(fileId,'  #   type = GenericConstantMaterial\n');
fprintf(fileId,'  #   prop_names = ''c_cu_eps c_cu_eta c_cu_sn c_eps_cu c_eps_eta c_eps_sn c_eta_cu c_eta_eps c_eta_sn c_sn_cu c_sn_eps c_sn_eta''\n');
fprintf(fileId,'  #   prop_values = ''0.02 0.1957 0.6088 0.2383 0.2483 0.2495 0.4299 0.4343 0.4359 0.9789 0.9839 0.9889'' #-\n');
fprintf(fileId,'  #   #prop_values = ''0.0234 0.198 0.6088 0.2479 0.2489 0.000 0.4345 0.4349 0.4351 0.9789 0.000 0.9889'' #- Aeps = 2e7 Aeta = 2e7\n');
fprintf(fileId,'  # [../]\n');
fprintf(fileId,'  #GIBBS ENERGY EQ. MOLAR CONC., eps phase is redundant\n');
fprintf(fileId,'  [./energy_constants_chat]\n');
fprintf(fileId,'    type = GenericConstantMaterial\n');
fprintf(fileId,'    prop_names = ''chat_cu_220   chat_eps_220  chat_eta_220  chat_sn_220\n');
fprintf(fileId,'                  chat_cu_25    chat_eps_25   chat_eta_25   chat_sn_25''\n');
fprintf(fileId,'    prop_values = ''0.02     0.2433    0.4351    0.9889\n');
fprintf(fileId,'                   0.10569  0.41753   0.41753   0.99941'' #-\n');
fprintf(fileId,'    #prop_values = ''0.0234 0.2484 0.4350 0.9889'' #-\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./diffusion_constants]\n');
fprintf(fileId,'    type = GenericConstantMaterial\n');
fprintf(fileId,'    prop_names = ''D_cu_220  D_eps_220   D_eta_220   D_sn_220\n');
fprintf(fileId,'                  D_cu_25   D_eps_25    D_eta_25    D_sn_25''\n');
fprintf(fileId,'    #prop_values = ''1e-20 6e-16 1.5e-14 1e-13'' # m^2/s #D12 best slightly slow\n');
fprintf(fileId,'    #prop_values = ''1e-20 9.5e-16 3e-14 1e-13'' # m^2/s #D15\n');
fprintf(fileId,'    prop_values = ''1e-20      1.25e-15  3.1e-14   1e-13\n');
fprintf(fileId,'                   2.877e-36  6.575e-19 6.575e-19 2.452e-17'' # m^2/s #D16 BEST\n');
fprintf(fileId,'    #prop_values = ''1e-16 1.25e-15 3.1e-14 1e-13'' # m^2/s #D17\n');
fprintf(fileId,'    #outputs = exodus\n');
fprintf(fileId,'  [../]\n\n');
fprintf(fileId,'  [./energy_chat_cu]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''chat_cu_220 chat_cu_25 phi_t''\n');
fprintf(fileId,'      function = ''chat_cu_25 - (chat_cu_25 - chat_cu_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = chat_cu\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_chat_eps]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''chat_eps_220 chat_eps_25 phi_t''\n');
fprintf(fileId,'      function = ''chat_eps_25 - (chat_eps_25 - chat_eps_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = chat_eps\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_chat_eta]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''chat_eta_220 chat_eta_25 phi_t''\n');
fprintf(fileId,'      function = ''chat_eta_25 - (chat_eta_25 - chat_eta_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = chat_eta\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./energy_chat_sn]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''chat_sn_220 chat_sn_25 phi_t''\n');
fprintf(fileId,'      function = ''chat_sn_25 - (chat_sn_25 - chat_sn_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = chat_sn\n');
fprintf(fileId,'  [../]\n\n');
fprintf(fileId,'  [./bulk_diff_D_cu]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''D_cu_220 D_cu_25 phi_t''\n');
fprintf(fileId,'      function = ''D_cu_25 - (D_cu_25 - D_cu_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = D_cu\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./bulk_diff_D_eps]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''D_eps_220 D_eps_25 phi_t''\n');
fprintf(fileId,'      function = ''D_eps_25 - (D_eps_25 - D_eps_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = D_eps\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./bulk_diff_D_eta]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''D_eta_220 D_eta_25 phi_t''\n');
fprintf(fileId,'      function = ''D_eta_25 - (D_eta_25 - D_eta_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = D_eta\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./bulk_diff_D_sn]\n');
fprintf(fileId,'      type =  ParsedMaterial\n');
fprintf(fileId,'      material_property_names = ''D_sn_220 D_sn_25 phi_t''\n');
fprintf(fileId,'      function = ''D_sn_25 - (D_sn_25 - D_sn_220)*(1-tanh(phi_t))/2''\n');
fprintf(fileId,'      f_name = D_sn\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./D_gb]\n');
fprintf(fileId,'    type = ParsedMaterial\n');
fprintf(fileId,'    material_property_names = ''D_eta D_sn ');
for i=1:(numCu6Sn5+numSn)
  if i>1 fprintf(fileId,' '); end
  fprintf(fileId,['h',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    f_name = D_gb\n');
fprintf(fileId,'    #function = ''200*D_eta''\n');
fprintf(fileId,'    #check if Cu6Sn5 grain boundary\n');
fprintf(fileId,'    function = ''if(');
for i=1:numCu6Sn5
   for j=(i+1):numCu6Sn5
   if i+j>3 fprintf(fileId,' | '); end
   fprintf(fileId,['h',num2str(i),'*h',num2str(j),'>0.02']);
   end
end   
fprintf(fileId,', 200*D_eta, 200*D_sn)''\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./fch_cu] #Chemical energy cu phase\n');
fprintf(fileId,'      type = DerivativeParsedMaterial\n');
fprintf(fileId,'      f_name = fch0\n');
fprintf(fileId,'      args = ''c0''\n');
fprintf(fileId,'      material_property_names = ''A_cu B_cu C_cu chat_cu length_scale energy_scale''\n');
fprintf(fileId,'      function = ''(energy_scale/length_scale^3)*(0.5*A_cu*(c0-chat_cu)^2+B_cu*(c0-chat_cu)+C_cu)'' #eV/nm^3\n');
fprintf(fileId,'      derivative_order = 2\n');
fprintf(fileId,'      use_displaced_mesh = true\n');
fprintf(fileId,'      outputs = exodus\n');
fprintf(fileId,'      output_properties = fch0\n');
fprintf(fileId,'  [../]\n');
for i=1:numCu6Sn5
  fprintf(fileId,['  [./fch_imc',num2str(i),'] #Chemical energy Cu6Sn5\n']);
  fprintf(fileId,'      type = DerivativeParsedMaterial\n');
  fprintf(fileId,['      f_name = fch',num2str(i),'\n']);
  fprintf(fileId,['      args = ''c',num2str(i),'''\n']);
  fprintf(fileId,'      material_property_names = ''A_eta B_eta C_eta chat_eta length_scale energy_scale''\n');
  fprintf(fileId,['      function = ''(energy_scale/length_scale^3)*(0.5*A_eta*(c',num2str(i),'-chat_eta)^2+B_eta*(c',num2str(i),'-chat_eta)+C_eta)'' #eV/nm^3\n']);
  fprintf(fileId,'      derivative_order = 2\n');
  fprintf(fileId,'      use_displaced_mesh = true\n');
  fprintf(fileId,'      outputs = exodus\n');
  fprintf(fileId,['      output_properties = fch',num2str(i),'\n']);
  fprintf(fileId,'  [../]\n');
end
for i=(1:numSn)+numCu6Sn5
fprintf(fileId,['  [./fch_sn',num2str(i),'] #Chemical energy Sn central grain\n']);
fprintf(fileId,'      type = DerivativeParsedMaterial\n');
fprintf(fileId,['      f_name = fch',num2str(i),'\n']);
fprintf(fileId,['      args = ''c',num2str(i),'''\n']);
fprintf(fileId,'      material_property_names = ''A_sn B_sn C_sn chat_sn length_scale energy_scale''\n');
fprintf(fileId,['      function = ''(energy_scale/length_scale^3)*(0.5*A_sn*(c',num2str(i),'-chat_sn)^2+B_sn*(c',num2str(i),'-chat_sn)+C_sn)'' #eV/nm^3\n']);
fprintf(fileId,'      derivative_order = 2\n');
fprintf(fileId,'      use_displaced_mesh = true\n');
fprintf(fileId,'      outputs = exodus\n');
fprintf(fileId,['      output_properties = fch',num2str(i),'\n']);
fprintf(fileId,'  [../]\n');
end
fprintf(fileId,'  [./Mgb]\n');
fprintf(fileId,'    type=DerivativeParsedMaterial\n');
fprintf(fileId,'    args = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0 fprintf(fileId,' '); end
  fprintf(fileId,['eta',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    material_property_names = ''D_gb delta delta_real A_cu A_eta A_sn length_scale energy_scale time_scale ');
for i=0:(numCu6Sn5+numSn);
    fprintf(fileId,[' h',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    f_name = Mgb\n');
fprintf(fileId,'    # Previuos lower boundary 0.09 testing 0.02 h2*h3>0.02\n');
%WRONG
fprintf(fileId,'    function = ''if(');
for i=1:numCu6Sn5
  for j=(i+1):numCu6Sn5
    if i+j>3 fprintf(fileId,' | '); end
    fprintf(fileId,['h',num2str(i),'*h',num2str(j),'>0.02']);
  end
end    
for i=(1:numSn)+numCu6Sn5
  for j=(i+1):(numSn+numCu6Sn5)
    if i+j>2*numCu6Sn5 fprintf(fileId,' | '); end
    fprintf(fileId,['h',num2str(i),'*h',num2str(j),'>0.02']);
  end
end    
fprintf(fileId,',(length_scale^5/(energy_scale*time_scale))*3.*D_gb*delta_real/((h0*A_cu+(');
for i=1:numCu6Sn5
  if i>1 fprintf(fileId,'+');  end
  fprintf(fileId,['h',num2str(i)]);
end
fprintf(fileId,')*A_eta+(');
for i=(1:numSn)+numCu6Sn5
  if i>numCu6Sn5+1 fprintf(fileId,'+');  end
  fprintf(fileId,['h',num2str(i)]);
end
fprintf(fileId,')*A_sn)*delta),0)''\n');
fprintf(fileId,'    outputs = exodus\n');
fprintf(fileId,'    output_properties = Mgb\n');
fprintf(fileId,'  [../]\n\n');
%
fprintf(fileId,'  [./Mbulk]\n');
fprintf(fileId,'      type = DerivativeParsedMaterial\n');
fprintf(fileId,'      f_name = Mbulk\n');
fprintf(fileId,'    args = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0 fprintf(fileId,' '); end
  fprintf(fileId,['eta',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,['      material_property_names = ''']);
for i=0:(numCu6Sn5+numSn)
  if i>0 fprintf(fileId,' '); end
  fprintf(fileId,['h',num2str(i)]);
end
fprintf(fileId,' D_cu D_eta D_sn A_cu A_eta A_sn Mgb length_scale energy_scale time_scale''\n');
fprintf(fileId,'      function = ''(length_scale^5/(energy_scale*time_scale))*(h0*D_cu/A_cu+(');
for i=1:numCu6Sn5
    if i>1 fprintf(fileId,'+'); end
    fprintf(fileId,['h',num2str(i)]);
end        
fprintf(fileId,')*D_eta/A_eta+(');
for i=(1:numSn)+numCu6Sn5
   if i>numCu6Sn5+1 fprintf(fileId,'+'); end
   fprintf(fileId,['h',num2str(i)]);
end
fprintf(fileId,')*D_sn/A_sn)'' #nm^5/eVs\n');
fprintf(fileId,'      derivative_order = 2\n');
fprintf(fileId,'      #outputs = exodus\n');
fprintf(fileId,'      #output_properties = M\n');
fprintf(fileId,'      use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./CHMobility]\n');
fprintf(fileId,'    type = DerivativeSumMaterial\n');
fprintf(fileId,'    f_name = M\n');
fprintf(fileId,'    sum_materials = ''Mgb Mbulk''\n');
fprintf(fileId,'    args = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0 fprintf(fileId,' '); end
  fprintf(fileId,['eta',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    outputs = exodus\n');
fprintf(fileId,'    output_properties = M\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n\n');
fprintf(fileId,'  [./F_cu]\n');
fprintf(fileId,'    type = DerivativeSumMaterial\n');
fprintf(fileId,'    f_name = F0\n');
fprintf(fileId,'    args = ''c0''\n');
fprintf(fileId,'    sum_materials = ''fch0 fe0''\n');
fprintf(fileId,'    #sum_materials = ''fch2''\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n');
for i=1:numCu6Sn5
fprintf(fileId,['  [./F_eta',num2str(i),']\n']);
fprintf(fileId,'    type = DerivativeSumMaterial\n');
fprintf(fileId,['    f_name = F',num2str(i),'\n']);
fprintf(fileId,['    args = ''c',num2str(i),' ']);
for j=0:(numCu6Sn5+numSn)
  fprintf(fileId,[' eta',num2str(j)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,['    sum_materials = ''fch',num2str(i),' fe',num2str(i),'''\n']);
fprintf(fileId,'    #sum_materials = ''fch1''\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n');
end
for i=(1:numSn)+numCu6Sn5
fprintf(fileId,['  [./F_sn',num2str(i),']\n']);
fprintf(fileId,'    type = DerivativeSumMaterial\n');
fprintf(fileId,['    f_name = F',num2str(i),'\n']);
fprintf(fileId,['    args = ''c',num2str(i),' ']);
for j=0:(numCu6Sn5+numSn)
  fprintf(fileId,[' eta',num2str(j)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,['    sum_materials = ''fch',num2str(i),' fe',num2str(i),' fp',num2str(i),'''\n']);
fprintf(fileId,'    #sum_materials = ''fch1''\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n');
end
fprintf(fileId,'\n');
%
fprintf(fileId,'  # Generate the global stress from the phase stresses\n');
fprintf(fileId,'  [./global_stress] #homogeniserar bara Cauchy stress\n');
fprintf(fileId,'    type = MultiPhaseStressMaterial\n');
fprintf(fileId,['    phase_base =''']);
for i=0:(numCu6Sn5+numSn)
  fprintf(fileId,[' eta',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,['    h       = ''']);
for i=0:(numCu6Sn5+numSn)
  fprintf(fileId,[' h',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    base_name = global\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./global_strain]\n');
fprintf(fileId,'    type = ComputeFiniteStrain\n');
fprintf(fileId,'    #displacements = ''disp_x disp_y disp_z''\n');
fprintf(fileId,'   base_name = global\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'[]\n\n');
%
% Kernels
%
fprintf(fileId,'[Kernels]\n');
fprintf(fileId,'  #Stress divergence\n');
fprintf(fileId,' [./TensorMechanics]\n');
fprintf(fileId,'    #displacements = ''disp_x disp_y disp_z''\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'    base_name = global\n');
fprintf(fileId,'    strain = FINITE #this also sets incremental strain =true\n');
fprintf(fileId,'  [../]\n\n');

for i=1:numCu6Sn5
fprintf(fileId,'  #Nucleation of Cu6Sn5, not activated\n');
fprintf(fileId,['  [./nuceta',num2str(i),']\n']);
fprintf(fileId,'    type = LangevinNoisePositive\n');
fprintf(fileId,['    variable = eta',num2str(i),'\n']);
fprintf(fileId,'    #amplitude = 1\n');
fprintf(fileId,'    amplitude = 0\n');
fprintf(fileId,'    seed = 1e9 #123456789\n');
fprintf(fileId,'    multiplier = nuc\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n');
end
fprintf(fileId,'\n');

fprintf(fileId,'  # Cahn-Hilliard Equation\n');
fprintf(fileId,'  [./CHBulk] # Gives the residual for the concentration, dF/dc-mu\n');
fprintf(fileId,'      type = KKSSplitCHCRes\n');
fprintf(fileId,'      variable = c\n');
fprintf(fileId,'      ca       = c1\n');
fprintf(fileId,'      cb       = c2\n');
fprintf(fileId,'      fa_name  = F1 #only fa is used\n');
fprintf(fileId,'      fb_name  = F2\n');
fprintf(fileId,'      w        = w\n');
fprintf(fileId,'      h_name   = h1\n');
fprintf(fileId,'      args_a = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0 fprintf(fileId,' '); end
  fprintf(fileId,['eta',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'      use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n\n');

fprintf(fileId,'  [./dcdt] # Gives dc/dt\n');
fprintf(fileId,'      type = CoupledTimeDerivative\n');
fprintf(fileId,'      variable = w\n');
fprintf(fileId,'      v = c\n');
fprintf(fileId,'      use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n\n');

fprintf(fileId,'  [./ckernel] # Gives residual for chemical potential dc/dt+M\\grad(mu)\n');
fprintf(fileId,'      type = SplitCHWRes\n');
fprintf(fileId,'      mob_name = M\n');
fprintf(fileId,'      variable = w\n');
fprintf(fileId,'    args = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0 fprintf(fileId,' '); end
  fprintf(fileId,['eta',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'      use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n\n');

fprintf(fileId,'  #KKS conditions\n');
fprintf(fileId,'  # enforce pointwise equality of chemical potentials, n-1 kernels needed (mu_1=mu_2, mu_2=mu_3, ..., mu_n-1=mu_n\n');
for i=0:(numCu6Sn5+numSn-1)
fprintf(fileId,['  [./chempot_cu_eta',num2str(i),']\n']);
fprintf(fileId,'    type = KKSPhaseChemicalPotential\n');
fprintf(fileId,['    variable = c',num2str(i),'\n']);
fprintf(fileId,['    cb       = c',num2str(i+1),'\n']);
fprintf(fileId,['    fa_name  = F',num2str(i),'\n']);
fprintf(fileId,['    fb_name  = F',num2str(i+1),'\n']);
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'    args_a = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0 fprintf(fileId,' '); end
  fprintf(fileId,['eta',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    args_b = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0 fprintf(fileId,' '); end
  fprintf(fileId,['eta',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'  [../]\n');
end
%
fprintf(fileId,'  [./phaseconcentration] # enforce c = sum h_i*c_i\n');
fprintf(fileId,'    type = KKSMultiPhaseConcentration\n');
fprintf(fileId,['    variable = c',num2str(numCu6Sn5+numSn),'\n']);
fprintf(fileId,'    cj = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0 fprintf(fileId,' '); end
  fprintf(fileId,['c',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    hj_names = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0 fprintf(fileId,' '); end
  fprintf(fileId,['h',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    etas = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0 fprintf(fileId,' '); end
  fprintf(fileId,['eta',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    c = c\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  # AC Kernels\n');
fprintf(fileId,'  [./KKSMultiACKernel]\n');
fprintf(fileId,['    op_num = ',num2str(numCu6Sn5+numSn+1),'\n']);
fprintf(fileId,'    op_name_base = ''eta''\n');
fprintf(fileId,'    ci_name_base = ''c''\n');
fprintf(fileId,'    f_name_base = F\n');
fprintf(fileId,'    #wi = 0.0624\n');
fprintf(fileId,'    #wi = 10.\n');
fprintf(fileId,'    #wi = 4.\n');
fprintf(fileId,'    wi = 0.\n');
fprintf(fileId,'    g_name_base = g\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'[]\n\n');
%
% [AuxVariables]
%
fprintf(fileId,'[AuxVariables]\n');
for i=0:(numSn+numCu6Sn5)
for j=(i+1):(numSn+numCu6Sn5)
fprintf(fileId,['  [./h',num2str(i),'h',num2str(j),']\n']); 
fprintf(fileId,'    order = CONSTANT\n');
fprintf(fileId,'    family = MONOMIAL\n');
fprintf(fileId,'    outputs = exodus\n');
fprintf(fileId,'    initial_condition = 0\n');
fprintf(fileId,'  [../]\n');
end
end
fprintf(fileId,'  [./f_int]\n');
fprintf(fileId,'    order = CONSTANT\n');
fprintf(fileId,'    family = MONOMIAL\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./f_density]\n');
fprintf(fileId,'    order = CONSTANT\n');
fprintf(fileId,'    family = MONOMIAL\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./step_dt]\n');
fprintf(fileId,'    order = CONSTANT\n');
fprintf(fileId,'    family = MONOMIAL\n');
fprintf(fileId,'  [../]\n');
for i=(1:numSn)+numCu6Sn5
  fprintf(fileId,['  [./sxx',num2str(i),']\n']);
fprintf(fileId,'    order = CONSTANT\n');
fprintf(fileId,'    family = MONOMIAL\n');
  fprintf(fileId,'  [../]\n');
end    
for i=(1:numSn)+numCu6Sn5
  fprintf(fileId,['  [./szz',num2str(i),']\n']);
fprintf(fileId,'    order = CONSTANT\n');
fprintf(fileId,'    family = MONOMIAL\n');
  fprintf(fileId,'  [../]\n');
end    
fprintf(fileId,'[]\n\n');
%
% AuxKernels
%
fprintf(fileId,'[AuxKernels]\n');
for i=0:(numSn+numCu6Sn5)
for j=(i+1):(numSn+numCu6Sn5)
fprintf(fileId,['  [./h',num2str(i),'h',num2str(j),']\n']); 
fprintf(fileId,'    type = ParsedAux\n');
fprintf(fileId,['    variable = h',num2str(i),'h',num2str(j),'\n']);
fprintf(fileId,['    args = ''h',num2str(i),' h',num2str(j),' h',num2str(i),'h',num2str(j),'''\n']);
fprintf(fileId,['    function = ''if(h',num2str(i),'*h',num2str(j),'>h',num2str(i),'h',num2str(j),',h',num2str(i),'*h',num2str(j),',h',num2str(i),'h',num2str(j),')''\n']);
fprintf(fileId,'    execute_on = ''LINEAR TIMESTEP_END''\n');
fprintf(fileId,'  [../]\n');
end
end
%
fprintf(fileId,'  [./step_dt]\n');
fprintf(fileId,'    type = MaterialRealAux\n');
fprintf(fileId,'    property = dt\n');
fprintf(fileId,'    variable = step_dt\n');
fprintf(fileId,'    execute_on = TIMESTEP_END\n');
fprintf(fileId,'  [../]\n');
% -Stresses in Sn
for i=(1:numSn)+numCu6Sn5
  fprintf(fileId,['  [./sxx',num2str(i),']\n']);
  fprintf(fileId,'    type = RankTwoAux\n');
  fprintf(fileId,['    variable = sxx',num2str(i),'\n']);
  fprintf(fileId,['    rank_two_tensor = eta',num2str(i),'_stress\n']);
  fprintf(fileId,'    index_i = 0\n');
  fprintf(fileId,'    index_j = 0\n');
  fprintf(fileId,'  [../]\n');
  fprintf(fileId,['  [./szz',num2str(i),']\n']);
  fprintf(fileId,'    type = RankTwoAux\n');
  fprintf(fileId,['    variable = szz',num2str(i),'\n']);
  fprintf(fileId,['    rank_two_tensor = eta',num2str(i),'_stress\n']);
  fprintf(fileId,'    index_i = 2\n');
  fprintf(fileId,'    index_j = 2\n');
  fprintf(fileId,'  [../]\n');
end
%
fprintf(fileId,'  [./f_density]\n');
fprintf(fileId,'    type = KKSMultiFreeEnergy\n');
fprintf(fileId,'    variable = f_density\n');
fprintf(fileId,'    hj_names = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0; fprintf(fileId,' '); end
  fprintf(fileId,['h',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    Fj_names = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0; fprintf(fileId,' '); end
  fprintf(fileId,['F',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    gj_names = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0; fprintf(fileId,' '); end
  fprintf(fileId,['g',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    additional_free_energy = f_int\n');
fprintf(fileId,'    interfacial_vars = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0; fprintf(fileId,' '); end
  fprintf(fileId,['eta',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    kappa_names = ''');
for i=0:(numCu6Sn5+numSn)
    fprintf(fileId,' kappa');
end
fprintf(fileId,'''\n');
fprintf(fileId,'    w = 0.\n');
fprintf(fileId,'    execute_on = ''initial timestep_end''\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./f_int]\n');
fprintf(fileId,'    type = ParsedAux\n');
fprintf(fileId,'    variable = f_int\n');
fprintf(fileId,'    args = ''');
for i=0:(numCu6Sn5+numSn)
  if i>0; fprintf(fileId,' '); end
  fprintf(fileId,['eta',num2str(i)]);
end
fprintf(fileId,'''\n');
fprintf(fileId,'    constant_names = ''sigma delta gamma length_scale energy_scale''\n');
fprintf(fileId,'    constant_expressions = ''0.5 ');
fprintf(fileId,ElementSize5);
fprintf(fileId,' 1.5 1e9 6.24150943e18''\n');
fprintf(fileId,'    function =''mu:=(6*sigma/delta)*(energy_scale/length_scale^3); mu*(');
for i=0:(numCu6Sn5+numSn)
    if i>0 fprintf(fileId,'+'); end;
    fprintf(fileId,['0.25*eta',num2str(i),'^4-0.5*eta',num2str(i),'^2']);
end
fprintf(fileId,'+gamma*(');
for i=0:(numCu6Sn5+numSn)
  for j=(i+1):(numCu6Sn5+numSn)
    if i+j>1 fprintf(fileId,'+'); end;
    fprintf(fileId,['eta',num2str(i),'^2*eta',num2str(j),'^2']);
  end
end
fprintf(fileId,')+0.25)''\n');
fprintf(fileId,'    execute_on = ''initial timestep_end''\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n'); 
fprintf(fileId,'[]\n\n');
%
fprintf(fileId,'[Postprocessors]\n');
fprintf(fileId,'  [./total_energy]\n');
fprintf(fileId,'    type = ElementIntegralVariablePostprocessor\n');
fprintf(fileId,'    variable = f_density\n');
fprintf(fileId,'    execute_on = ''Initial TIMESTEP_END''\n');
fprintf(fileId,'    use_displaced_mesh = true\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./step_size]\n');
fprintf(fileId,'    type = TimestepSize\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'[]\n\n');
%
fprintf(fileId,'[Debug]\n');
fprintf(fileId,'  show_var_residual_norms = true\n');
fprintf(fileId,'  show_material_props = true\n');
fprintf(fileId,'[]\n\n');
%
fprintf(fileId,'[Preconditioning]\n');
fprintf(fileId,'  [./smp]\n');
fprintf(fileId,'    type = SMP\n');
fprintf(fileId,'    full = true\n');
fprintf(fileId,'    solve_type = PJFNK\n');
fprintf(fileId,'    petsc_options_iname = ''-pc_asm_overlap -pc_type -sub_pc_type  -pc_factor_shift_type  -sub_pc_factor_shift_type -sub_ksp_type''\n');
fprintf(fileId,'    petsc_options_value = ''2                  asm       lu             nonzero             nonzero                    preonly''\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'[]\n\n');
%
fprintf(fileId,'[Executioner]\n');
fprintf(fileId,'   type = Transient\n');
fprintf(fileId,'   solve_type = ''PJFNK''\n');
fprintf(fileId,'   line_search = basic\n');
fprintf(fileId,'   petsc_options = ''-snes_converged_reason -ksp_converged_reason -snes_ksp_ew''\n');
fprintf(fileId,'   l_max_its = 10\n');
fprintf(fileId,'   nl_max_its = 25\n');
fprintf(fileId,'   l_tol = 1.0e-3\n');
fprintf(fileId,'   nl_rel_tol = 1.0e-6\n');
fprintf(fileId,'   nl_abs_tol = 1.0e-9\n');
fprintf(fileId,'   #num_steps = 0\n');
fprintf(fileId,'   end_time = 1.2672e6 #2x176h\n');
fprintf(fileId,'   dtmax= 1e5\n\n');
%
fprintf(fileId,'   [./TimeIntegrator]\n');
fprintf(fileId,'       type = LStableDirk2\n');
fprintf(fileId,'   [../]\n\n');
%
fprintf(fileId,'   [./TimeStepper]\n');
fprintf(fileId,'       # Turn on time stepping\n');
fprintf(fileId,'       type = IterationAdaptiveDT\n');
fprintf(fileId,'       dt = 5e-3\n');
fprintf(fileId,'       cutback_factor = 0.5\n');
fprintf(fileId,'       growth_factor = 1.3\n');
fprintf(fileId,'       optimal_iterations = 15 # 20\n');
fprintf(fileId,'       linear_iteration_ratio = 25\n');
fprintf(fileId,'   [../]\n');
fprintf(fileId,'[]\n\n');
%
fprintf(fileId,'[Outputs]\n');
fprintf(fileId,'  [./exodus]\n');
fprintf(fileId,'    type = Exodus\n');
fprintf(fileId,'    file_base = CuSn6Sn5Sn\n');
fprintf(fileId,'    append_date = true\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./perf_log]\n');
fprintf(fileId,'    type =  CSV\n');
fprintf(fileId,'    file_base = perf_log\n');
fprintf(fileId,'    append_date = true\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./sim_log]\n');
fprintf(fileId,'    type = Console\n');
fprintf(fileId,'    append_date = true\n');
fprintf(fileId,'    output_file = true\n');
fprintf(fileId,'    file_base = console_sim_log\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'  [./checkpoint]\n');
fprintf(fileId,'    type = Checkpoint\n');
fprintf(fileId,'    num_files = 2\n');
fprintf(fileId,'    suffix = cp\n');
fprintf(fileId,'    use_displaced = true\n');
fprintf(fileId,'  [../]\n');
fprintf(fileId,'[]\n\n');
%
%
fclose(fileId);
