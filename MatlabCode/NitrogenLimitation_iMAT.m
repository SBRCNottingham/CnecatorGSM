

reac_expression_data = xlsread('CnecatorH16_ReactionExpression.xlsx');
exp_f26_RPKM = reac_expression_data(:,2);
exp_f16_RPKM = reac_expression_data(:,1);

m = readCbModel('../Model/iCN1361.xml');
m = changeRxnBounds(m, 'R_QUIN_OXIDASE_b03', 0, 'b');
m = changeRxnBounds(m, 'R_QUIN_OXIDASE_bd', 0, 'b');
m = changeRxnBounds(m, 'R_Biomass', 0.2, 'b');
m = changeRxnBounds(m, 'EX_BETA_D_FRUCTOSE_e', -2.105, 'b');
m = changeRxnBounds(m, 'R_FRUPTS', 0.0, 'u');



m_nh3 = changeRxnBounds(m, 'R_Biomass', 0.009, 'l');
m_nh3 = changeRxnBounds(m_nh3, 'EX_AMMONIUM_e', -0.1, 'b');
m_nh3 = changeRxnBounds(m_nh3, 'EX_BETA_D_FRUCTOSE_e', -2.105, 'b');
m_nh3 = changeRxnBounds(m_nh3, 'R_FRUPTS', 0.0, 'u');


iMAT_model_f16_212 = iMAT(m, exp_f16_RPKM, 212, 11439, 1e-8, [], '', 1800);
iMAT_model_f26_21 = iMAT(m_nh3, exp_f26_RPKM, 21, 17535, 1e-8, [], '', 1800);




