
[DBG][ACT] Action Dependency Sets:
[DBG][ACT] ([33mno_action, setup_oversampling, deprecated_block, finish_input_file_output, meta_action[39m)
[DBG][ACT]	[32mPeriodic[39m
[DBG][ACT]	EmptyAction
[DBG][ACT] ([33mdynamic_object_registration[39m)
[DBG][ACT]	[32mProblem[39m
[DBG][ACT]	DynamicObjectRegistrationAction
[DBG][ACT] ([33mcommon_output[39m)
[DBG][ACT]	[32mOutputs[39m
[DBG][ACT]	CommonOutputAction
[DBG][ACT] ([33mset_global_params[39m)
[DBG][ACT] ([33msetup_recover_file_base[39m)
[DBG][ACT]	SetupRecoverFileBaseAction
[DBG][ACT] ([33mcheck_copy_nodal_vars[39m)
[DBG][ACT]	[32mc[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mw[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mc1[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mc2[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mc3[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mc0[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mf_density[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mf_int[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32ms[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT] ([33msetup_mesh[39m)
[DBG][ACT]	[32mMesh[39m
[DBG][ACT]	SetupMeshAction ([36msetup_mesh[35m, init_mesh[39m)
[DBG][ACT] ([33madd_partitioner[39m)
[DBG][ACT] ([33minit_mesh[39m)
[DBG][ACT]	[32mMesh[39m
[DBG][ACT]	SetupMeshAction ([36minit_mesh[35m, setup_mesh[39m)
[DBG][ACT] ([33mprepare_mesh[39m)
[DBG][ACT]	[32mMesh[39m
[DBG][ACT]	SetupMeshCompleteAction ([36mprepare_mesh[35m, execute_mesh_modifiers, setup_mesh_complete, uniform_refine_mesh[39m)
[DBG][ACT] ([33madd_mesh_modifier[39m)
[DBG][ACT] ([33mexecute_mesh_modifiers[39m)
[DBG][ACT]	[32mMesh[39m
[DBG][ACT]	SetupMeshCompleteAction ([36mexecute_mesh_modifiers[35m, prepare_mesh, setup_mesh_complete, uniform_refine_mesh[39m)
[DBG][ACT] ([33madd_mortar_interface[39m)
[DBG][ACT] ([33muniform_refine_mesh[39m)
[DBG][ACT]	[32mMesh[39m
[DBG][ACT]	SetupMeshCompleteAction ([36muniform_refine_mesh[35m, execute_mesh_modifiers, prepare_mesh, setup_mesh_complete[39m)
[DBG][ACT] ([33msetup_mesh_complete[39m)
[DBG][ACT]	[32mMesh[39m
[DBG][ACT]	SetupMeshCompleteAction ([36msetup_mesh_complete[35m, execute_mesh_modifiers, prepare_mesh, uniform_refine_mesh[39m)
[DBG][ACT] ([33mdetermine_system_type[39m)
[DBG][ACT]	[32mExecutioner[39m
[DBG][ACT]	DetermineSystemType
[DBG][ACT] ([33mcreate_problem[39m)
[DBG][ACT]	[32mProblem[39m
[DBG][ACT]	CreateProblemAction
[DBG][ACT] ([33mvalidate_coordinate_systems, setup_time_integrator[39m)
[DBG][ACT] ([33msetup_executioner[39m)
[DBG][ACT]	[32mExecutioner[39m
[DBG][ACT]	CreateExecutionerAction
[DBG][ACT] ([33msetup_time_stepper[39m)
[DBG][ACT]	[32mTimeStepper[39m
[DBG][ACT]	SetupTimeStepperAction
[DBG][ACT] ([33msetup_predictor[39m)
[DBG][ACT] ([33msetup_postprocessor_data[39m)
[DBG][ACT] ([33minit_displaced_problem[39m)
[DBG][ACT]	[32mMesh[39m
[DBG][ACT]	CreateDisplacedProblemAction
[DBG][ACT] ([33madd_elemental_field_variable, add_aux_variable, add_variable[39m)
[DBG][ACT]	[32mf_density[39m
[DBG][ACT]	AddAuxVariableAction
[DBG][ACT]	[32mf_int[39m
[DBG][ACT]	AddAuxVariableAction
[DBG][ACT]	[32ms[39m
[DBG][ACT]	AddAuxVariableAction
[DBG][ACT]	[32mc[39m
[DBG][ACT]	AddVariableAction
[DBG][ACT]	[32mw[39m
[DBG][ACT]	AddVariableAction
[DBG][ACT]	[32mc1[39m
[DBG][ACT]	AddVariableAction
[DBG][ACT]	[32mc2[39m
[DBG][ACT]	AddVariableAction
[DBG][ACT]	[32mc3[39m
[DBG][ACT]	AddVariableAction
[DBG][ACT]	[32mc0[39m
[DBG][ACT]	AddVariableAction
[DBG][ACT]	[32mPolycrystalVariables[39m
[DBG][ACT]	PolycrystalVariablesAction
[DBG][ACT] ([33msetup_variable_complete[39m)
[DBG][ACT] ([33msetup_quadrature[39m)
[DBG][ACT]	SetupQuadratureAction
[DBG][ACT] ([33madd_function[39m)
[DBG][ACT] ([33madd_distribution[39m)
[DBG][ACT] ([33madd_periodic_bc[39m)
[DBG][ACT]	[32mx[39m
[DBG][ACT]	AddPeriodicBCAction
[DBG][ACT] ([33madd_user_object[39m)
[DBG][ACT] ([33msetup_function_complete[39m)
[DBG][ACT] ([33msetup_adaptivity[39m)
[DBG][ACT] ([33mset_adaptivity_options[39m)
[DBG][ACT] ([33madd_ic[39m)
[DBG][ACT]	[32meta1[39m
[DBG][ACT]	AddInitialConditionAction
[DBG][ACT]	[32meta2[39m
[DBG][ACT]	AddInitialConditionAction
[DBG][ACT]	[32meta3[39m
[DBG][ACT]	AddInitialConditionAction
[DBG][ACT]	[32meta0[39m
[DBG][ACT]	AddInitialConditionAction
[DBG][ACT]	[32mc[39m
[DBG][ACT]	AddInitialConditionAction
[DBG][ACT] ([33madd_constraint, add_field_split[39m)
[DBG][ACT] ([33madd_preconditioning[39m)
[DBG][ACT]	[32mfull[39m
[DBG][ACT]	SetupPreconditionerAction
[DBG][ACT] ([33mready_to_init[39m)
[DBG][ACT]	EmptyAction
[DBG][ACT] ([33msetup_dampers[39m)
[DBG][ACT]	SetupDampersAction
[DBG][ACT] ([33msetup_residual_debug[39m)
[DBG][ACT]	[32mDebug[39m
[DBG][ACT]	SetupResidualDebugAction
[DBG][ACT] ([33madd_bounds_vectors[39m)
[DBG][ACT] ([33madd_multi_app[39m)
[DBG][ACT] ([33madd_transfer[39m)
[DBG][ACT] ([33mcopy_nodal_vars, copy_nodal_aux_vars[39m)
[DBG][ACT]	[32mc[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mw[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mc1[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mc2[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mc3[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mc0[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mf_density[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32mf_int[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT]	[32ms[39m
[DBG][ACT]	CopyNodalVarsAction
[DBG][ACT] ([33madd_material[39m)
[DBG][ACT]	[32mscale[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mmodel_constants[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32menergy_A[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32menergy_B[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32menergy_C[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32menergy_c_ab[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32menergy_chat[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mdiffusion_constants[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mD_gb[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mkappa[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mmu[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mL_imc0[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mfch1[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mfch2[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mfch3[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mfch0[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mh_cu[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mh_imc1[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mh_imc2[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mh_sn[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mg_cu[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mg_imc1[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mg_imc2[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mg_sn[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mCHMobility[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT]	[32mACMobility[39m
[DBG][ACT]	AddMaterialAction
[DBG][ACT] ([33msetup_material_output[39m)
[DBG][ACT]	MaterialOutputAction
[DBG][ACT] ([33minit_problem[39m)
[DBG][ACT]	InitProblemAction
[DBG][ACT] ([33msetup_debug[39m)
[DBG][ACT]	[32mDebug[39m
[DBG][ACT]	SetupDebugAction
[DBG][ACT] ([33madd_output[39m)
[DBG][ACT]	[32mexodus_out[39m
[DBG][ACT]	AddOutputAction
[DBG][ACT] ([33madd_postprocessor[39m)
[DBG][ACT] ([33madd_vector_postprocessor[39m)
[DBG][ACT] ([33madd_kernel, add_nodal_kernel, add_bc, add_aux_kernel, add_scalar_kernel, add_aux_scalar_kernel, add_dirac_kernel, add_dg_kernel, add_interface_kernel, add_damper, add_indicator, add_marker[39m)
[DBG][ACT]	[32mCHBulk[39m
[DBG][ACT]	AddKernelAction
[DBG][ACT]	[32mdcdt[39m
[DBG][ACT]	AddKernelAction
[DBG][ACT]	[32mckernel[39m
[DBG][ACT]	AddKernelAction
[DBG][ACT]	[32mchempot_cu_imc[39m
[DBG][ACT]	AddKernelAction
[DBG][ACT]	[32mchempot_imc_imc[39m
[DBG][ACT]	AddKernelAction
[DBG][ACT]	[32mchempot_sn_cu[39m
[DBG][ACT]	AddKernelAction
[DBG][ACT]	[32mphaseconcentration[39m
[DBG][ACT]	AddKernelAction
[DBG][ACT]	[32mKKSMultiACKernel[39m
[DBG][ACT]	KKSMultiACKernelAction
[DBG][ACT]	[32mneumann[39m
[DBG][ACT]	AddBCAction
[DBG][ACT] ([33madd_control[39m)
[DBG][ACT] ([33mcheck_output[39m)
[DBG][ACT]	CheckOutputAction
[DBG][ACT] ([33mcheck_integrity[39m)
[DBG][ACT]	CheckIntegrityAction


[DBG][ACT] Executing actions:
[DBG][ACT] TASK ([33m               no_action[39m) TYPE ([33m                     EmptyAction[39m) NAME ([33m        Periodic[39m) 
[DBG][ACT] TASK ([33mdynamic_object_registration[39m) TYPE ([33m DynamicObjectRegistrationAction[39m) NAME ([33m         Problem[39m) 
[DBG][ACT] TASK ([33m           common_output[39m) TYPE ([33m              CommonOutputAction[39m) NAME ([33m         Outputs[39m) 
[DBG][ACT] TASK ([33m setup_recover_file_base[39m) TYPE ([33m      SetupRecoverFileBaseAction[39m) NAME ([33m                [39m) 
[DBG][ACT] TASK ([33m   check_copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m               c[39m) 
[DBG][ACT] TASK ([33m   check_copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m               w[39m) 
[DBG][ACT] TASK ([33m   check_copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m              c1[39m) 
[DBG][ACT] TASK ([33m   check_copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m              c2[39m) 
[DBG][ACT] TASK ([33m   check_copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m              c3[39m) 
[DBG][ACT] TASK ([33m   check_copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m              c0[39m) 
[DBG][ACT] TASK ([33m   check_copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m       f_density[39m) 
[DBG][ACT] TASK ([33m   check_copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m           f_int[39m) 
[DBG][ACT] TASK ([33m   check_copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m               s[39m) 
[DBG][ACT] TASK ([33m              setup_mesh[39m) TYPE ([33m                 SetupMeshAction[39m) NAME ([33m            Mesh[39m) 
[DBG][ACT] TASK ([33m               init_mesh[39m) TYPE ([33m                 SetupMeshAction[39m) NAME ([33m            Mesh[39m) 
[DBG][ACT] TASK ([33m            prepare_mesh[39m) TYPE ([33m         SetupMeshCompleteAction[39m) NAME ([33m            Mesh[39m) 
[DBG][ACT] TASK ([33m  execute_mesh_modifiers[39m) TYPE ([33m         SetupMeshCompleteAction[39m) NAME ([33m            Mesh[39m) 
[DBG][ACT] TASK ([33m     uniform_refine_mesh[39m) TYPE ([33m         SetupMeshCompleteAction[39m) NAME ([33m            Mesh[39m) 
[DBG][ACT] TASK ([33m     setup_mesh_complete[39m) TYPE ([33m         SetupMeshCompleteAction[39m) NAME ([33m            Mesh[39m) 
[DBG][ACT] TASK ([33m   determine_system_type[39m) TYPE ([33m             DetermineSystemType[39m) NAME ([33m     Executioner[39m) 
[DBG][ACT] TASK ([33m          create_problem[39m) TYPE ([33m             CreateProblemAction[39m) NAME ([33m         Problem[39m) 
[DBG][ACT] TASK ([33m       setup_executioner[39m) TYPE ([33m         CreateExecutionerAction[39m) NAME ([33m     Executioner[39m) 
[DBG][ACT] TASK ([33m      setup_time_stepper[39m) TYPE ([33m          SetupTimeStepperAction[39m) NAME ([33m     TimeStepper[39m) 
[DBG][ACT] TASK ([33m  init_displaced_problem[39m) TYPE ([33m    CreateDisplacedProblemAction[39m) NAME ([33m            Mesh[39m) 
[DBG][ACT] TASK ([33m        add_aux_variable[39m) TYPE ([33m            AddAuxVariableAction[39m) NAME ([33m       f_density[39m) 
[DBG][ACT] TASK ([33m        add_aux_variable[39m) TYPE ([33m            AddAuxVariableAction[39m) NAME ([33m           f_int[39m) 
[DBG][ACT] TASK ([33m        add_aux_variable[39m) TYPE ([33m            AddAuxVariableAction[39m) NAME ([33m               s[39m) 
[DBG][ACT] TASK ([33m            add_variable[39m) TYPE ([33m               AddVariableAction[39m) NAME ([33m               c[39m) 
[DBG][ACT] TASK ([33m            add_variable[39m) TYPE ([33m               AddVariableAction[39m) NAME ([33m               w[39m) 
[DBG][ACT] TASK ([33m            add_variable[39m) TYPE ([33m               AddVariableAction[39m) NAME ([33m              c1[39m) 
[DBG][ACT] TASK ([33m            add_variable[39m) TYPE ([33m               AddVariableAction[39m) NAME ([33m              c2[39m) 
[DBG][ACT] TASK ([33m            add_variable[39m) TYPE ([33m               AddVariableAction[39m) NAME ([33m              c3[39m) 
[DBG][ACT] TASK ([33m            add_variable[39m) TYPE ([33m               AddVariableAction[39m) NAME ([33m              c0[39m) 
[DBG][ACT] TASK ([33m            add_variable[39m) TYPE ([33m      PolycrystalVariablesAction[39m) NAME ([33mPolycrystalVariables[39m) 
[DBG][ACT] TASK ([33m        setup_quadrature[39m) TYPE ([33m           SetupQuadratureAction[39m) NAME ([33m                [39m) 
[DBG][ACT] TASK ([33m         add_periodic_bc[39m) TYPE ([33m             AddPeriodicBCAction[39m) NAME ([33m               x[39m) 
[DBG][ACT] TASK ([33m                  add_ic[39m) TYPE ([33m       AddInitialConditionAction[39m) NAME ([33m            eta1[39m) 
[DBG][ACT] TASK ([33m                  add_ic[39m) TYPE ([33m       AddInitialConditionAction[39m) NAME ([33m            eta2[39m) 
[DBG][ACT] TASK ([33m                  add_ic[39m) TYPE ([33m       AddInitialConditionAction[39m) NAME ([33m            eta3[39m) 
[DBG][ACT] TASK ([33m                  add_ic[39m) TYPE ([33m       AddInitialConditionAction[39m) NAME ([33m            eta0[39m) 
[DBG][ACT] TASK ([33m                  add_ic[39m) TYPE ([33m       AddInitialConditionAction[39m) NAME ([33m               c[39m) 
[DBG][ACT] TASK ([33m                  add_ic[39m) TYPE ([33m       AddInitialConditionAction[39m) NAME ([33m        c1_moose[39m) 
[DBG][ACT] TASK ([33m                  add_ic[39m) TYPE ([33m       AddInitialConditionAction[39m) NAME ([33m        c2_moose[39m) 
[DBG][ACT] TASK ([33m                  add_ic[39m) TYPE ([33m       AddInitialConditionAction[39m) NAME ([33m        c3_moose[39m) 
[DBG][ACT] TASK ([33m                  add_ic[39m) TYPE ([33m       AddInitialConditionAction[39m) NAME ([33m        c0_moose[39m) 
[DBG][ACT] TASK ([33m     add_preconditioning[39m) TYPE ([33m       SetupPreconditionerAction[39m) NAME ([33m            full[39m) 
[DBG][ACT] TASK ([33m           ready_to_init[39m) TYPE ([33m                     EmptyAction[39m) NAME ([33m                [39m) 
[DBG][ACT] TASK ([33m           setup_dampers[39m) TYPE ([33m              SetupDampersAction[39m) NAME ([33m                [39m) 
[DBG][ACT] TASK ([33m    setup_residual_debug[39m) TYPE ([33m        SetupResidualDebugAction[39m) NAME ([33m           Debug[39m) 
[DBG][ACT] TASK ([33m         copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m               c[39m) 
[DBG][ACT] TASK ([33m         copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m               w[39m) 
[DBG][ACT] TASK ([33m         copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m              c1[39m) 
[DBG][ACT] TASK ([33m         copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m              c2[39m) 
[DBG][ACT] TASK ([33m         copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m              c3[39m) 
[DBG][ACT] TASK ([33m         copy_nodal_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m              c0[39m) 
[DBG][ACT] TASK ([33m     copy_nodal_aux_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m       f_density[39m) 
[DBG][ACT] TASK ([33m     copy_nodal_aux_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m           f_int[39m) 
[DBG][ACT] TASK ([33m     copy_nodal_aux_vars[39m) TYPE ([33m             CopyNodalVarsAction[39m) NAME ([33m               s[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m           scale[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m model_constants[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m        energy_A[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m        energy_B[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m        energy_C[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m     energy_c_ab[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m     energy_chat[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33mdiffusion_constants[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m            D_gb[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m           kappa[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m              mu[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m          L_imc0[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m            fch1[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m            fch2[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m            fch3[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m            fch0[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m            h_cu[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m          h_imc1[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m          h_imc2[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m            h_sn[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m            g_cu[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m          g_imc1[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m          g_imc2[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m            g_sn[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m      CHMobility[39m) 
[DBG][ACT] TASK ([33m            add_material[39m) TYPE ([33m               AddMaterialAction[39m) NAME ([33m      ACMobility[39m) 
[DBG][ACT] TASK ([33m   setup_material_output[39m) TYPE ([33m            MaterialOutputAction[39m) NAME ([33m                [39m) 
[DBG][ACT] TASK ([33m            init_problem[39m) TYPE ([33m               InitProblemAction[39m) NAME ([33m                [39m) 
[DBG][ACT] TASK ([33m             setup_debug[39m) TYPE ([33m                SetupDebugAction[39m) NAME ([33m           Debug[39m) 
[DBG][ACT] TASK ([33m              add_output[39m) TYPE ([33m                 AddOutputAction[39m) NAME ([33m      exodus_out[39m) 
[DBG][ACT] TASK ([33m              add_output[39m) TYPE ([33m                 AddOutputAction[39m) NAME ([33m         console[39m) 
[DBG][ACT] TASK ([33m              add_output[39m) TYPE ([33m                 AddOutputAction[39m) NAME ([33m             csv[39m) 
[DBG][ACT] TASK ([33m              add_output[39m) TYPE ([33m                 AddOutputAction[39m) NAME ([33m_moose_variable_residual_norms_debug_output[39m) 
[DBG][ACT] TASK ([33m              add_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m          CHBulk[39m) 
[DBG][ACT] TASK ([33m              add_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m            dcdt[39m) 
[DBG][ACT] TASK ([33m              add_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m         ckernel[39m) 
[DBG][ACT] TASK ([33m              add_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m  chempot_cu_imc[39m) 
[DBG][ACT] TASK ([33m              add_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m chempot_imc_imc[39m) 
[DBG][ACT] TASK ([33m              add_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m   chempot_sn_cu[39m) 
[DBG][ACT] TASK ([33m              add_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mphaseconcentration[39m) 
[DBG][ACT] TASK ([33m              add_kernel[39m) TYPE ([33m          KKSMultiACKernelAction[39m) NAME ([33mKKSMultiACKernel[39m) 
[DBG][ACT] TASK ([33m                  add_bc[39m) TYPE ([33m                     AddBCAction[39m) NAME ([33m         neumann[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mfch1_d^2fch1/dc1^2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m  fch1_dfch1/dc1[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m       fch1_fch1[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mfch2_d^2fch2/dc2^2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m  fch2_dfch2/dc2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m       fch2_fch2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mfch3_d^2fch3/dc3^2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m  fch3_dfch3/dc3[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m       fch3_fch3[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mfch0_d^2fch0/dc0^2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m  fch0_dfch0/dc0[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m       fch0_fch0[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m    CHMobility_M[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_dM/deta0[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_dM/deta1[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_dM/deta2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_dM/deta3[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_d^2M/deta0^2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_d^2M/deta0deta1[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_d^2M/deta0deta2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_d^2M/deta0deta3[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_d^2M/deta1^2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_d^2M/deta1deta2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_d^2M/deta1deta3[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_d^2M/deta2^2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_d^2M/deta2deta3[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33mCHMobility_d^2M/deta3^2[39m) 
[DBG][ACT] TASK ([33m          add_aux_kernel[39m) TYPE ([33m                 AddKernelAction[39m) NAME ([33m    ACMobility_L[39m) 
[DBG][ACT] TASK ([33m            check_output[39m) TYPE ([33m               CheckOutputAction[39m) NAME ([33m                [39m) 
[DBG][ACT] TASK ([33m         check_integrity[39m) TYPE ([33m            CheckIntegrityAction[39m) NAME ([33m                [39m) 
Framework Information:
MOOSE version:           git commit 1456dfb on 2017-06-15
PETSc Version:           3.7.5
Current Time:            Mon Nov  6 10:21:51 2017
Executable Timestamp:    Mon Nov  6 09:40:51 2017

Parallelism:
  Num Processors:          16
  Num Threads:             1

Mesh: 
  Parallel Type:           replicated
  Mesh Dimension:          2
  Spatial Dimension:       2
  Nodes:                   
    Total:                 2501
    Local:                 180
  Elems:                   
    Total:                 2400
    Local:                 154
  Num Subdomains:          1
  Num Partitions:          16
  Partitioner:             metis

Nonlinear System:
  Num DOFs:                25010
  Num Local DOFs:          1800
  Variables:               { "c" "w" "c1" "c2" "c3" "c0" "eta0" "eta1" "eta2" "eta3" } 
  Finite Element Types:    "LAGRANGE" 
  Approximation Orders:    "FIRST" 

Auxiliary System:
  Num DOFs:                74501
  Num Local DOFs:          4800
  Variables:               { "f_density" "f_int" } "s" { "L" "M" "dM/deta0" "dM/deta1" "dM/deta2" "dM/deta3" 
                             "d^2M/deta0^2" "d^2M/deta0deta1" "d^2M/deta0deta2" "d^2M/deta0deta3" "d^2M/deta1^2" 
                             "d^2M/deta1deta2" "d^2M/deta1deta3" "d^2M/deta2^2" "d^2M/deta2deta3" "d^2M/deta3^2" 
                             "d^2fch0/dc0^2" "d^2fch1/dc1^2" "d^2fch2/dc2^2" "d^2fch3/dc3^2" "dfch0/dc0" 
                             "dfch1/dc1" "dfch2/dc2" "dfch3/dc3" "fch0" "fch1" "fch2" "fch3" } 
  Finite Element Types:    "MONOMIAL" "LAGRANGE" "MONOMIAL" 
  Approximation Orders:    "CONSTANT" "FIRST" "CONSTANT" 

Execution Information:
  Executioner:             Transient
  TimeStepper:             IterationAdaptiveDT
  Solver Mode:             Preconditioned JFNK



--- action-test.i ------------------------------------------------------
[]
  element_order                  = AUTO
  order                          = AUTO
  side_order                     = AUTO
  type                           = GAUSS
[]

[AuxVariables]

  [./f_density]
    block                        = INVALID
    family                       = MONOMIAL
    initial_condition            = INVALID
    order                        = CONSTANT
    outputs                      = INVALID
    initial_from_file_timestep   = LATEST
    initial_from_file_var        = INVALID
  [../]

  [./f_int]
    block                        = INVALID
    family                       = MONOMIAL
    initial_condition            = INVALID
    order                        = CONSTANT
    outputs                      = INVALID
    initial_from_file_timestep   = LATEST
    initial_from_file_var        = INVALID
  [../]

  [./s]
    block                        = INVALID
    family                       = LAGRANGE
    initial_condition            = INVALID
    order                        = FIRST
    outputs                      = INVALID
    initial_from_file_timestep   = LATEST
    initial_from_file_var        = INVALID
  [../]
[]

[BCs]

  [./Periodic]

    [./x]
      auto_direction             = x
      inv_transform_func         = INVALID
      primary                    = INVALID
      secondary                  = INVALID
      transform_func             = INVALID
      translation                = INVALID
      variable                   = 'c w c1 c2 c3 c0 eta1 eta2 eta3 eta0'
    [../]
  [../]

  [./neumann]
    boundary                     = 'top bottom'
    control_tags                 = INVALID
    enable                       = 1
    implicit                     = 1
    type                         = NeumannBC
    use_displaced_mesh           = 0
    variable                     = c
    diag_save_in                 = INVALID
    save_in                      = INVALID
    seed                         = 0
    value                        = 0
  [../]
[]

[Debug]
  show_actions                   = 1
  show_material_props            = 0
  show_parser                    = 0
  show_top_residuals             = 0
  show_var_residual_norms        = 1
  show_var_residual              = INVALID
[]

[Executioner]
  type                           = Transient
  abort_on_solve_fail            = 0
  compute_initial_residual_before_preset_bcs = 0
  control_tags                   = 
  dt                             = 1
  dtmax                          = 1e+30
  dtmin                          = 2e-14
  enable                         = 1
  end_time                       = 2.16e+06
  l_abs_step_tol                 = -1
  l_max_its                      = 30
  l_tol                          = 0.0001
  line_search                    = default
  max_xfem_update                = 4294967295
  n_startup_steps                = 0
  nl_abs_step_tol                = 1e-50
  nl_abs_tol                     = 1e-07
  nl_max_funcs                   = 10000
  nl_max_its                     = 10
  nl_rel_step_tol                = 1e-50
  nl_rel_tol                     = 1e-08
  no_fe_reinit                   = 0
  num_steps                      = 4294967295
  petsc_options                  = '-KSP_CONVERGED_REASON -SNES_CONVERGED_REASON'
  petsc_options_iname            = '-PC_TYPE -SUB_PC_FACTOR_SHIFT_TYPE -SUB_PC_TYPE'
  petsc_options_value            = 'asm ilu nonzero'
  picard_abs_tol                 = 1e-50
  picard_max_its                 = 1
  picard_rel_tol                 = 1e-08
  reset_dt                       = 0
  restart_file_base              = 
  scheme                         = INVALID
  solve_type                     = PJFNK
  splitting                      = INVALID
  ss_check_tol                   = 1e-08
  ss_tmin                        = 0
  start_time                     = 0
  time_period_ends               = INVALID
  time_period_starts             = INVALID
  time_periods                   = INVALID
  timestep_tolerance             = 2e-14
  trans_ss_check                 = 0
  use_multiapp_dt                = 0
  verbose                        = 0

  [./TimeStepper]
    type                         = IterationAdaptiveDT
    _executioner                 = 0x1db8f10
    _fe_problem_base             = 0x1d544d0
    control_tags                 = Executioner
    cutback_factor               = 0.2
    dt                           = 1
    enable                       = 1
    force_step_every_function_point = 0
    growth_factor                = 2
    iteration_window             = INVALID
    linear_iteration_ratio       = INVALID
    max_function_change          = INVALID
    optimal_iterations           = 5
    postprocessor_dtlim          = INVALID
    reset_dt                     = 0
    time_dt                      = INVALID
    time_t                       = INVALID
    timestep_limiting_function   = INVALID
  [../]
[]

[Executioner]
  _fe_problem                    = 0x1d544d0
  _fe_problem_base               = 0x1d544d0

  [./TimeStepper]
  [../]
[]

[ICs]

  [./c]
    type                         = MultiBoundingBoxIC
    block                        = INVALID
    boundary                     = INVALID
    control_tags                 = ICs
    corners                      = '(x,y,z)=(       0,        0,        0) (x,y,z)=(       0,      650,        0) (x,y,z)=(       0,     1050,        0)'
    enable                       = 1
    ignore_uo_dependency         = 0
    inside                       = '0.002 0.417 0.999'
    opposite_corners             = '(x,y,z)=(    2000,      600,        0) (x,y,z)=(    2000,     1000,        0) (x,y,z)=(    2000,     3000,        0)'
    outside                      = 0
    variable                     = c
  [../]

  [./eta0]
    type                         = MultiBoundingBoxIC
    block                        = INVALID
    boundary                     = INVALID
    control_tags                 = ICs
    corners                      = '(x,y,z)=(       0,     1050,        0)'
    enable                       = 1
    ignore_uo_dependency         = 0
    inside                       = 1
    opposite_corners             = '(x,y,z)=(    2000,     3000,        0)'
    outside                      = 0
    variable                     = eta0
  [../]

  [./eta1]
    type                         = MultiBoundingBoxIC
    block                        = INVALID
    boundary                     = INVALID
    control_tags                 = ICs
    corners                      = '(x,y,z)=(       0,        0,        0)'
    enable                       = 1
    ignore_uo_dependency         = 0
    inside                       = 1
    opposite_corners             = '(x,y,z)=(    2000,      600,        0)'
    outside                      = 0
    variable                     = eta1
  [../]

  [./eta2]
    type                         = MultiBoundingBoxIC
    block                        = INVALID
    boundary                     = INVALID
    control_tags                 = ICs
    corners                      = '(x,y,z)=(       0,      650,        0) (x,y,z)=(    1550,      650,        0)'
    enable                       = 1
    ignore_uo_dependency         = 0
    inside                       = 1
    opposite_corners             = '(x,y,z)=(     500,     1000,        0) (x,y,z)=(    2000,     1000,        0)'
    outside                      = 0
    variable                     = eta2
  [../]

  [./eta3]
    type                         = MultiBoundingBoxIC
    block                        = INVALID
    boundary                     = INVALID
    control_tags                 = ICs
    corners                      = '(x,y,z)=(     550,      650,        0)'
    enable                       = 1
    ignore_uo_dependency         = 0
    inside                       = 1
    opposite_corners             = '(x,y,z)=(    1500,     1000,        0)'
    outside                      = 0
    variable                     = eta3
  [../]
[]

[Kernels]

  [./CHBulk]
    type                         = KKSSplitCHCRes
    args_a                       = INVALID
    block                        = INVALID
    ca                           = c2
    cb                           = c0
    control_tags                 = Kernels
    diag_save_in                 = INVALID
    eigen_kernel                 = 0
    enable                       = 1
    fa_name                      = fch2
    fb_name                      = fch0
    h_name                       = h2
    implicit                     = 1
    save_in                      = INVALID
    seed                         = 0
    use_displaced_mesh           = 0
    variable                     = c
    w                            = w
  [../]

  [./KKSMultiACKernel]
    ci_name_base                 = c
    fch_name_base                = fch
    g_name_base                  = g
    gamma                        = gamma
    h_name_base                  = h
    implicit                     = 1
    kappa                        = kappa
    mob_name                     = L
    op_name_base                 = eta
    op_num                       = 4
    use_displaced_mesh           = 0
    wi                           = 0
  [../]

  [./chempot_cu_imc]
    type                         = KKSPhaseChemicalPotential
    args_a                       = INVALID
    args_b                       = INVALID
    block                        = INVALID
    cb                           = c2
    control_tags                 = Kernels
    diag_save_in                 = INVALID
    eigen_kernel                 = 0
    enable                       = 1
    fa_name                      = fch1
    fb_name                      = fch2
    implicit                     = 1
    save_in                      = INVALID
    seed                         = 0
    use_displaced_mesh           = 0
    variable                     = c1
  [../]

  [./chempot_imc_imc]
    type                         = KKSPhaseChemicalPotential
    args_a                       = INVALID
    args_b                       = INVALID
    block                        = INVALID
    cb                           = c3
    control_tags                 = Kernels
    diag_save_in                 = INVALID
    eigen_kernel                 = 0
    enable                       = 1
    fa_name                      = fch2
    fb_name                      = fch3
    implicit                     = 1
    save_in                      = INVALID
    seed                         = 0
    use_displaced_mesh           = 0
    variable                     = c2
  [../]

  [./chempot_sn_cu]
    type                         = KKSPhaseChemicalPotential
    args_a                       = INVALID
    args_b                       = INVALID
    block                        = INVALID
    cb                           = c0
    control_tags                 = Kernels
    diag_save_in                 = INVALID
    eigen_kernel                 = 0
    enable                       = 1
    fa_name                      = fch3
    fb_name                      = fch0
    implicit                     = 1
    save_in                      = INVALID
    seed                         = 0
    use_displaced_mesh           = 0
    variable                     = c3
  [../]

  [./ckernel]
    type                         = SplitCHWRes
    args                         = 'eta1 eta2 eta3 eta0'
    block                        = INVALID
    control_tags                 = Kernels
    diag_save_in                 = INVALID
    eigen_kernel                 = 0
    enable                       = 1
    implicit                     = 1
    mob_name                     = M
    save_in                      = INVALID
    seed                         = 0
    use_displaced_mesh           = 0
    variable                     = w
  [../]

  [./dcdt]
    type                         = CoupledTimeDerivative
    block                        = INVALID
    control_tags                 = Kernels
    diag_save_in                 = INVALID
    eigen_kernel                 = 0
    enable                       = 1
    implicit                     = 1
    save_in                      = INVALID
    seed                         = 0
    use_displaced_mesh           = 0
    v                            = c
    variable                     = w
  [../]

  [./phaseconcentration]
    type                         = KKSMultiPhaseConcentration
    block                        = INVALID
    c                            = c
    cj                           = 'c1 c2 c3 c0'
    control_tags                 = Kernels
    diag_save_in                 = INVALID
    eigen_kernel                 = 0
    enable                       = 1
    etas                         = 'eta1 eta2 eta3 eta0'
    hj_names                     = 'h1 h2 h3 h0'
    implicit                     = 1
    save_in                      = INVALID
    seed                         = 0
    use_displaced_mesh           = 0
    variable                     = c0
  [../]
[]

[Materials]

  [./ACMobility]
    type                         = DerivativeParsedMaterial
    args                         = 'eta1 eta2 eta3 eta0'
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_expressions         = 
    constant_names               = 
    constant_on_elem             = 0
    control_tags                 = Materials
    derivative_order             = 2
    disable_fpoptimizer          = 0
    enable                       = 1
    enable_ad_cache              = 1
    enable_auto_optimize         = 1
    enable_jit                   = 1
    f_name                       = L
    fail_on_evalerror            = 0
    function                     = L_imc0
    implicit                     = 1
    material_property_names      = 'L_cu_imc L_imc0 L_cu_sn L_imc_imc'
    output_properties            = INVALID
    outputs                      = exodus_out
    seed                         = 0
    third_derivatives            = INVALID
    tol_names                    = 
    tol_values                   = 
    use_displaced_mesh           = 0
  [../]

  [./CHMobility]
    type                         = DerivativeParsedMaterial
    args                         = 'eta1 eta2 eta3 eta0'
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_expressions         = 
    constant_names               = 
    constant_on_elem             = 0
    control_tags                 = Materials
    derivative_order             = 2
    disable_fpoptimizer          = 0
    enable                       = 1
    enable_ad_cache              = 1
    enable_auto_optimize         = 1
    enable_jit                   = 1
    f_name                       = M
    fail_on_evalerror            = 0
    function                     = (length_scale^5/(energy_scale*time_scale))*(h1*D_cu/A_cu+h2*D_imc/A_imc+h3*D_imc/A_imc+h0*D_sn/A_sn)
    implicit                     = 1
    material_property_names      = 'h1(eta1,eta2,eta3,eta0) h2(eta1,eta2,eta3,eta0) h3(eta1,eta2,eta3,eta0) h0(eta1,eta2,eta3,eta0) D_cu D_imc D_sn A_cu A_imc A_sn length_scale energy_scale time_scale Mgb'
    output_properties            = INVALID
    outputs                      = exodus_out
    seed                         = 0
    third_derivatives            = INVALID
    tol_names                    = 
    tol_values                   = 
    use_displaced_mesh           = 0
  [../]

  [./D_gb]
    type                         = ParsedMaterial
    args                         = INVALID
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_expressions         = 
    constant_names               = 
    constant_on_elem             = 0
    control_tags                 = Materials
    disable_fpoptimizer          = 0
    enable                       = 1
    enable_ad_cache              = 1
    enable_auto_optimize         = 1
    enable_jit                   = 1
    f_name                       = D_gb
    fail_on_evalerror            = 0
    function                     = 200*D_imc
    implicit                     = 1
    material_property_names      = D_imc
    output_properties            = INVALID
    outputs                      = none
    seed                         = 0
    tol_names                    = 
    tol_values                   = 
    use_displaced_mesh           = 0
  [../]

  [./L_imc0]
    type                         = ParsedMaterial
    args                         = INVALID
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_expressions         = 
    constant_names               = 
    constant_on_elem             = 0
    control_tags                 = Materials
    disable_fpoptimizer          = 0
    enable                       = 1
    enable_ad_cache              = 1
    enable_auto_optimize         = 1
    enable_jit                   = 1
    f_name                       = L_imc0
    fail_on_evalerror            = 0
    function                     = (length_scale^5/(energy_scale*time_scale))*2*mu*(D_sn/A_sn+D_imc/A_imc)/(3*kappa*(c_imc_sn-c_sn_imc)^2)
    implicit                     = 1
    material_property_names      = 'mu kappa D_sn D_imc A_sn A_imc c_imc_sn c_sn_imc length_scale energy_scale time_scale'
    output_properties            = INVALID
    outputs                      = none
    seed                         = 0
    tol_names                    = 
    tol_values                   = 
    use_displaced_mesh           = 0
  [../]

  [./diffusion_constants]
    type                         = GenericConstantMaterial
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    prop_names                   = 'D_cu D_imc D_sn'
    prop_values                  = '2.877e-36 6.575e-19 2.452e-17'
    seed                         = 0
    use_displaced_mesh           = 0
  [../]

  [./energy_A]
    type                         = GenericConstantMaterial
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    prop_names                   = 'A_cu A_imc A_sn'
    prop_values                  = '6.2204e+09 2.4555e+10 2.5819e+11'
    seed                         = 0
    use_displaced_mesh           = 0
  [../]

  [./energy_B]
    type                         = GenericConstantMaterial
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    prop_names                   = 'B_cu B_imc B_sn'
    prop_values                  = '-1.2981e+09 -4.29e+08 4.4e+08'
    seed                         = 0
    use_displaced_mesh           = 0
  [../]

  [./energy_C]
    type                         = GenericConstantMaterial
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    prop_names                   = 'Cf_cu Cf_imc Cf_sn'
    prop_values                  = '-7.883e+08 -1.178e+09 -9.4e+08'
    seed                         = 0
    use_displaced_mesh           = 0
  [../]

  [./energy_c_ab]
    type                         = GenericConstantMaterial
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    prop_names                   = 'c_imc_sn c_sn_imc'
    prop_values                  = '0.4529 0.9994'
    seed                         = 0
    use_displaced_mesh           = 0
  [../]

  [./energy_chat]
    type                         = GenericConstantMaterial
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    prop_names                   = 'chat_cu chat_imc chat_sn'
    prop_values                  = '0.1057 0.41753 0.9994'
    seed                         = 0
    use_displaced_mesh           = 0
  [../]

  [./fch0]
    type                         = DerivativeParsedMaterial
    args                         = c0
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_expressions         = 
    constant_names               = 
    constant_on_elem             = 0
    control_tags                 = Materials
    derivative_order             = 2
    disable_fpoptimizer          = 0
    enable                       = 1
    enable_ad_cache              = 1
    enable_auto_optimize         = 1
    enable_jit                   = 1
    f_name                       = fch0
    fail_on_evalerror            = 0
    function                     = (energy_scale/length_scale^3)*(0.5*A_sn*(c0-chat_sn)^2+B_sn*(c0-chat_sn)+Cf_sn)
    implicit                     = 1
    material_property_names      = 'A_sn B_sn Cf_sn chat_sn length_scale energy_scale'
    output_properties            = INVALID
    outputs                      = exodus_out
    seed                         = 0
    third_derivatives            = INVALID
    tol_names                    = 
    tol_values                   = 
    use_displaced_mesh           = 0
  [../]

  [./fch1]
    type                         = DerivativeParsedMaterial
    args                         = c1
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_expressions         = 
    constant_names               = 
    constant_on_elem             = 0
    control_tags                 = Materials
    derivative_order             = 2
    disable_fpoptimizer          = 0
    enable                       = 1
    enable_ad_cache              = 1
    enable_auto_optimize         = 1
    enable_jit                   = 1
    f_name                       = fch1
    fail_on_evalerror            = 0
    function                     = (energy_scale/length_scale^3)*(0.5*A_cu*(c1-chat_cu)^2+B_cu*(c1-chat_cu)+Cf_cu)
    implicit                     = 1
    material_property_names      = 'A_cu B_cu Cf_cu chat_cu length_scale energy_scale'
    output_properties            = INVALID
    outputs                      = exodus_out
    seed                         = 0
    third_derivatives            = INVALID
    tol_names                    = 
    tol_values                   = 
    use_displaced_mesh           = 0
  [../]

  [./fch2]
    type                         = DerivativeParsedMaterial
    args                         = c2
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_expressions         = 
    constant_names               = 
    constant_on_elem             = 0
    control_tags                 = Materials
    derivative_order             = 2
    disable_fpoptimizer          = 0
    enable                       = 1
    enable_ad_cache              = 1
    enable_auto_optimize         = 1
    enable_jit                   = 1
    f_name                       = fch2
    fail_on_evalerror            = 0
    function                     = (energy_scale/length_scale^3)*(0.5*A_imc*(c2-chat_imc)^2+B_imc*(c2-chat_imc)+Cf_imc)
    implicit                     = 1
    material_property_names      = 'A_imc B_imc Cf_imc chat_imc length_scale energy_scale'
    output_properties            = INVALID
    outputs                      = exodus_out
    seed                         = 0
    third_derivatives            = INVALID
    tol_names                    = 
    tol_values                   = 
    use_displaced_mesh           = 0
  [../]

  [./fch3]
    type                         = DerivativeParsedMaterial
    args                         = c3
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_expressions         = 
    constant_names               = 
    constant_on_elem             = 0
    control_tags                 = Materials
    derivative_order             = 2
    disable_fpoptimizer          = 0
    enable                       = 1
    enable_ad_cache              = 1
    enable_auto_optimize         = 1
    enable_jit                   = 1
    f_name                       = fch3
    fail_on_evalerror            = 0
    function                     = (energy_scale/length_scale^3)*(0.5*A_imc*(c3-chat_imc)^2+B_imc*(c3-chat_imc)+Cf_imc)
    implicit                     = 1
    material_property_names      = 'A_imc B_imc Cf_imc chat_imc length_scale energy_scale'
    output_properties            = INVALID
    outputs                      = exodus_out
    seed                         = 0
    third_derivatives            = INVALID
    tol_names                    = 
    tol_values                   = 
    use_displaced_mesh           = 0
  [../]

  [./g_cu]
    type                         = BarrierFunctionMaterial
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    eta                          = eta1
    function_name                = g1
    g_order                      = SIMPLE
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    seed                         = 0
    use_displaced_mesh           = 0
    well_only                    = 1
  [../]

  [./g_imc1]
    type                         = BarrierFunctionMaterial
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    eta                          = eta2
    function_name                = g2
    g_order                      = SIMPLE
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    seed                         = 0
    use_displaced_mesh           = 0
    well_only                    = 1
  [../]

  [./g_imc2]
    type                         = BarrierFunctionMaterial
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    eta                          = eta3
    function_name                = g3
    g_order                      = SIMPLE
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    seed                         = 0
    use_displaced_mesh           = 0
    well_only                    = 1
  [../]

  [./g_sn]
    type                         = BarrierFunctionMaterial
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    eta                          = eta0
    function_name                = g0
    g_order                      = SIMPLE
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    seed                         = 0
    use_displaced_mesh           = 0
    well_only                    = 1
  [../]

  [./h_cu]
    type                         = SwitchingFunctionMultiPhaseMaterial
    all_etas                     = 'eta1 eta2 eta3 eta0'
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    h_name                       = h1
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    phase_etas                   = eta1
    seed                         = 0
    use_displaced_mesh           = 0
  [../]

  [./h_imc1]
    type                         = SwitchingFunctionMultiPhaseMaterial
    all_etas                     = 'eta1 eta2 eta3 eta0'
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    h_name                       = h2
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    phase_etas                   = eta2
    seed                         = 0
    use_displaced_mesh           = 0
  [../]

  [./h_imc2]
    type                         = SwitchingFunctionMultiPhaseMaterial
    all_etas                     = 'eta1 eta2 eta3 eta0'
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    h_name                       = h3
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    phase_etas                   = eta3
    seed                         = 0
    use_displaced_mesh           = 0
  [../]

  [./h_sn]
    type                         = SwitchingFunctionMultiPhaseMaterial
    all_etas                     = 'eta1 eta2 eta3 eta0'
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    h_name                       = h0
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    phase_etas                   = eta0
    seed                         = 0
    use_displaced_mesh           = 0
  [../]

  [./kappa]
    type                         = ParsedMaterial
    args                         = INVALID
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_expressions         = 
    constant_names               = 
    constant_on_elem             = 0
    control_tags                 = Materials
    disable_fpoptimizer          = 0
    enable                       = 1
    enable_ad_cache              = 1
    enable_auto_optimize         = 1
    enable_jit                   = 1
    f_name                       = kappa
    fail_on_evalerror            = 0
    function                     = 0.75*sigma*delta*energy_scale/length_scale
    implicit                     = 1
    material_property_names      = 'sigma delta length_scale energy_scale'
    output_properties            = INVALID
    outputs                      = none
    seed                         = 0
    tol_names                    = 
    tol_values                   = 
    use_displaced_mesh           = 0
  [../]

  [./model_constants]
    type                         = GenericConstantMaterial
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    prop_names                   = 'sigma delta delta_real gamma Vm tgrad_corr_mult'
    prop_values                  = '0.5 3e-07 5e-10 1.5 1.629e-05 0'
    seed                         = 0
    use_displaced_mesh           = 0
  [../]

  [./mu]
    type                         = ParsedMaterial
    args                         = INVALID
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_expressions         = 
    constant_names               = 
    constant_on_elem             = 0
    control_tags                 = Materials
    disable_fpoptimizer          = 0
    enable                       = 1
    enable_ad_cache              = 1
    enable_auto_optimize         = 1
    enable_jit                   = 1
    f_name                       = mu
    fail_on_evalerror            = 0
    function                     = 6*(sigma/delta)*energy_scale/length_scale^3
    implicit                     = 1
    material_property_names      = 'sigma delta length_scale energy_scale'
    output_properties            = INVALID
    outputs                      = none
    seed                         = 0
    tol_names                    = 
    tol_values                   = 
    use_displaced_mesh           = 0
  [../]

  [./scale]
    type                         = GenericConstantMaterial
    block                        = INVALID
    boundary                     = INVALID
    compute                      = 1
    constant_on_elem             = 0
    control_tags                 = Materials
    enable                       = 1
    implicit                     = 1
    output_properties            = INVALID
    outputs                      = none
    prop_names                   = 'length_scale energy_scale time_scale'
    prop_values                  = '1e+09 6.24151e+18 1'
    seed                         = 0
    use_displaced_mesh           = 0
  [../]
[]

[Mesh]
  displacements                  = INVALID
  block_id                       = INVALID
  block_name                     = INVALID
  boundary_id                    = INVALID
  boundary_name                  = INVALID
  construct_side_list_from_node_list = 0
  ghosted_boundaries             = INVALID
  ghosted_boundaries_inflation   = INVALID
  patch_size                     = 40
  second_order                   = 0
  skip_partitioning              = 0
  type                           = GeneratedMesh
  uniform_refine                 = 0
  allow_renumbering              = 1
  bias_x                         = 1
  bias_y                         = 1
  bias_z                         = 1
  centroid_partitioner_direction = INVALID
  construct_node_list_from_side_list = 1
  control_tags                   = 
  dim                            = 2
  distribution                   = DEFAULT
  elem_type                      = QUAD4
  enable                         = 1
  gauss_lobatto_grid             = 0
  ghost_point_neighbors          = 0
  nemesis                        = 0
  num_ghosted_layers             = 1
  nx                             = 40
  ny                             = 60
  nz                             = 1
  parallel_type                  = DEFAULT
  partitioner                    = default
  patch_update_strategy          = never
  xmax                           = 2000
  xmin                           = 0
  ymax                           = 3000
  ymin                           = 0
  zmax                           = 1
  zmin                           = 0
[]

[Mesh]
[]

[Outputs]
  append_date                    = 0
  append_date_format             = INVALID
  checkpoint                     = 0
  color                          = 1
  console                        = 1
  controls                       = 0
  csv                            = 1
  dofmap                         = 0
  execute_on                     = 'INITIAL TIMESTEP_END'
  exodus                         = 0
  file_base                      = kernelaction
  gmv                            = 0
  gnuplot                        = 0
  hide                           = INVALID
  interval                       = 1
  nemesis                        = 0
  output_if_base_contains        = INVALID
  print_linear_residuals         = 1
  print_mesh_changed_info        = 0
  print_perf_log                 = 1
  show                           = INVALID
  solution_history               = 0
  sync_times                     = 
  tecplot                        = 0
  vtk                            = 0
  xda                            = 0
  xdr                            = 0

  [./exodus_out]
    type                         = Exodus
    additional_execute_on        = INVALID
    append_date                  = 0
    append_date_format           = INVALID
    append_oversample            = 0
    control_tags                 = Outputs
    elemental_as_nodal           = 0
    enable                       = 1
    end_time                     = INVALID
    execute_elemental_on         = INVALID
    execute_elemental_variables  = 1
    execute_input                = 1
    execute_input_on             = INITIAL
    execute_nodal_on             = INVALID
    execute_nodal_variables      = 1
    execute_on                   = 'INITIAL TIMESTEP_END'
    execute_postprocessors_on    = INVALID
    execute_scalar_variables     = 1
    execute_scalars_on           = INVALID
    execute_system_information   = 1
    execute_vector_postprocessors = 1
    file                         = INVALID
    file_base                    = kernelaction
    hide                         = INVALID
    interval                     = 1
    linear_residual_dt_divisor   = 1000
    linear_residual_end_time     = INVALID
    linear_residual_start_time   = INVALID
    linear_residuals             = 0
    nonlinear_residual_dt_divisor = 1000
    nonlinear_residual_end_time  = INVALID
    nonlinear_residual_start_time = INVALID
    nonlinear_residuals          = 0
    output_if_base_contains      = 
    output_linear                = 0
    output_material_properties   = 0
    output_nonlinear             = 0
    output_postprocessors        = 1
    oversample                   = 0
    overwrite                    = 0
    padding                      = 3
    position                     = INVALID
    refinements                  = 0
    scalar_as_nodal              = 0
    sequence                     = INVALID
    show                         = INVALID
    show_material_properties     = INVALID
    start_time                   = INVALID
    sync_only                    = 0
    sync_times                   = 
    time_tolerance               = 1e-14
    use_displaced                = 0
    use_problem_dimension        = INVALID
  [../]
[]

[Preconditioning]

  [./full]
    type                         = SMP
    control_tags                 = Preconditioning
    coupled_groups               = INVALID
    enable                       = 1
    full                         = 1
    ksp_norm                     = unpreconditioned
    line_search                  = default
    off_diag_column              = INVALID
    off_diag_row                 = INVALID
    pc_side                      = default
    petsc_options                = INVALID
    petsc_options_iname          = INVALID
    petsc_options_value          = INVALID
    solve_type                   = INVALID
  [../]
[]

[Problem]
  block                          = INVALID
  coord_type                     = XYZ
  fe_cache                       = 0
  kernel_coverage_check          = 1
  material_coverage_check        = 1
  name                           = 'MOOSE Problem'
  restart_file_base              = INVALID
  rz_coord_axis                  = Y
  type                           = FEProblem
  library_path                   = 
  object_names                   = INVALID
  register_objects_from          = INVALID
  control_tags                   = 
  enable                         = 1
  error_on_jacobian_nonzero_reallocation = 0
  force_restart                  = 0
  near_null_space_dimension      = 0
  null_space_dimension           = 0
  petsc_inames                   = 
  petsc_options                  = INVALID
  petsc_values                   = 
  solve                          = 0
  transpose_null_space_dimension = 0
  use_nonlinear                  = 1
[]

[Variables]

  [./PolycrystalVariables]
    family                       = LAGRANGE
    op_num                       = 4
    order                        = FIRST
    scaling                      = 1
    var_name_base                = eta
  [../]

  [./c]
    block                        = INVALID
    eigen                        = 0
    family                       = LAGRANGE
    initial_condition            = INVALID
    order                        = FIRST
    outputs                      = INVALID
    scaling                      = 1
    initial_from_file_timestep   = LATEST
    initial_from_file_var        = INVALID
  [../]

  [./c0]
    block                        = INVALID
    eigen                        = 0
    family                       = LAGRANGE
    initial_condition            = 0.999
    order                        = FIRST
    outputs                      = INVALID
    scaling                      = 1
    initial_from_file_timestep   = LATEST
    initial_from_file_var        = INVALID
  [../]

  [./c1]
    block                        = INVALID
    eigen                        = 0
    family                       = LAGRANGE
    initial_condition            = 0.002
    order                        = FIRST
    outputs                      = INVALID
    scaling                      = 1
    initial_from_file_timestep   = LATEST
    initial_from_file_var        = INVALID
  [../]

  [./c2]
    block                        = INVALID
    eigen                        = 0
    family                       = LAGRANGE
    initial_condition            = 0.417
    order                        = FIRST
    outputs                      = INVALID
    scaling                      = 1
    initial_from_file_timestep   = LATEST
    initial_from_file_var        = INVALID
  [../]

  [./c3]
    block                        = INVALID
    eigen                        = 0
    family                       = LAGRANGE
    initial_condition            = 0.417
    order                        = FIRST
    outputs                      = INVALID
    scaling                      = 1
    initial_from_file_timestep   = LATEST
    initial_from_file_var        = INVALID
  [../]

  [./w]
    block                        = INVALID
    eigen                        = 0
    family                       = LAGRANGE
    initial_condition            = INVALID
    order                        = FIRST
    outputs                      = INVALID
    scaling                      = 1
    initial_from_file_timestep   = LATEST
    initial_from_file_var        = INVALID
  [../]
[]


Time Step  0, time = 0
                dt = 0

Time Step  1, time = 1
                dt = 1
[32m Solve Converged![39m

Time Step  2, time = 3
                dt = 2
[32m Solve Converged![39m

Time Step  3, time = 7
                dt = 4
[32m Solve Converged![39m

Time Step  4, time = 15
                dt = 8
[32m Solve Converged![39m

Time Step  5, time = 31
                dt = 16
[32m Solve Converged![39m

Time Step  6, time = 63
                dt = 32
[32m Solve Converged![39m

Time Step  7, time = 127
                dt = 64
[32m Solve Converged![39m

Time Step  8, time = 255
                dt = 128
[32m Solve Converged![39m

Time Step  9, time = 511
                dt = 256
[32m Solve Converged![39m

Time Step 10, time = 1023
                dt = 512
[32m Solve Converged![39m

Time Step 11, time = 2047
                dt = 1024
[32m Solve Converged![39m

Time Step 12, time = 4095
                dt = 2048
[32m Solve Converged![39m

Time Step 13, time = 8191
                dt = 4096
[32m Solve Converged![39m

Time Step 14, time = 16383
                dt = 8192
[32m Solve Converged![39m

Time Step 15, time = 32767
                dt = 16384
[32m Solve Converged![39m

Time Step 16, time = 65535
                dt = 32768
[32m Solve Converged![39m

Time Step 17, time = 131071
                dt = 65536
[32m Solve Converged![39m

Time Step 18, time = 262143
                dt = 131072
[32m Solve Converged![39m

Time Step 19, time = 524287
                dt = 262144
[32m Solve Converged![39m

Time Step 20, time = 1.04858e+06
                dt = 524288
[32m Solve Converged![39m

Time Step 21, time = 2.09715e+06
                dt = 1.04858e+06
[32m Solve Converged![39m

Time Step 22, time = 2.16e+06
                dt = 62849
[32m Solve Converged![39m

 -------------------------------------------------------------------------------------------------------------------------
| Puffin Performance: Alive time=18.9903, Active time=18.9197                                                             |
 -------------------------------------------------------------------------------------------------------------------------
| Event                                      nCalls     Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                       w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------------------|
|                                                                                                                         |
|                                                                                                                         |
| Application                                                                                                             |
|   Full Runtime                             1          0.0654      0.065405    18.9198     18.919753   0.35     100.00   |
|                                                                                                                         |
| Execution                                                                                                               |
|   computeElemAux(TIMESTEP_END)             22         0.3633      0.016513    0.3633      0.016513    1.92     1.92     |
|   solve()                                  22         0.0018      0.000082    0.0018      0.000082    0.01     0.01     |
|                                                                                                                         |
| Output                                                                                                                  |
|   CSV::output()                            46         0.0002      0.000005    0.0002      0.000005    0.00     0.00     |
|   Exodus::output()                         23         0.3132      0.013618    0.3132      0.013618    1.66     1.66     |
|                                                                                                                         |
| Setup                                                                                                                   |
|   Application Setup                        1          17.2866     17.286590   18.1140     18.113987   91.37    95.74    |
|   FEProblemBase::init::meshChanged()       1          0.0015      0.001492    0.0015      0.001492    0.01     0.01     |
|   Initial updateActiveSemiLocalNodeRange() 1          0.0001      0.000059    0.0001      0.000059    0.00     0.00     |
|   Initial updateGeomSearch()               2          0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   NonlinearSystem::update()                1          0.0010      0.000986    0.0010      0.000986    0.01     0.01     |
|   eq.init()                                1          0.8249      0.824918    0.8249      0.824918    4.36     4.36     |
|   execMultiApps()                          1          0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   execTransfers()                          1          0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   initial adaptivity                       1          0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   initialSetup()                           1          0.0555      0.055483    0.0618      0.061792    0.29     0.33     |
|   reinit() after updateGeomSearch()        1          0.0007      0.000651    0.0007      0.000651    0.00     0.00     |
|                                                                                                                         |
| Utility                                                                                                                 |
|   projectSolution()                        1          0.0056      0.005591    0.0056      0.005591    0.03     0.03     |
 -------------------------------------------------------------------------------------------------------------------------
| Totals:                                    127        18.9197                                         100.00            |
 -------------------------------------------------------------------------------------------------------------------------
