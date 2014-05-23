function [VAR] = INITSTRUCT(VAR,Type)
if strcmp(Type,'Sample')
    VAR.sample_i = [];
    VAR.imbal_s = [];
    VAR.Pi_st = [];
    VAR.actuals_s = [];
    VAR.price_actual_s = [];
    VAR.price_posActual_s = [];
    VAR.price_negActual_s = [];
    VAR.price_gas_s = [];
    VAR.price_elecC_s = [];
    VAR.price_elecS_s = [];
    VAR.elecD_s = [];
    VAR.elecWF_s = [];
    VAR.elecWP_s = [];
    VAR.heatD_s = [];
    return
end

if strcmp(Type,'CHPdata')
    VAR.R_b = [];
    VAR.R_ir = [];
    VAR.FC_bc = [];
    VAR.m_fCHP = [];
    VAR.m_fB = [];
    VAR.Q_CHP = [];
    VAR.Q_B = [];
    VAR.DeltaQ_S = [];
    VAR.Q_S = [];
    VAR.E_CHP = [];
    VAR.E_ir = [];
    VAR.E_b = [];
    VAR.ON = [];
    return
end
error('INITSTRUCT: Wrong type');
end