# This function is not used currently.
# To enable it please change the line 61 and 86-87 in simplemcgui\pages\secondary_page.py


def write_base_paramDefs(base_param_value, base_param_error, base_param_bound):
    """
    Writing the base parameters to the paramDefs.py file.

    Args:
        base_param_value [list]: Value of base parameters
        base_param_error [list]: Error in base parameters
        base_param_bound [list]: Bounds of base parameters
    """
    with open('simplemc\cosmo\paramDefs.py', 'r+') as param_file:
        data = param_file.read()
        param_file.seek(0)
        param_file.write('from simplemc.cosmo.Parameter import Parameter\n\n')
        param_file.write('Om_par = Parameter("Om", {0}, {1}, {2}, "\Omega_m")\n'.format(base_param_value[0], base_param_error[0], base_param_bound[0]))
        param_file.write('Obh2_par = Parameter("Obh2", {0}, {1}, {2}, "\Omega_{{b}}h^2")\n'.format(base_param_value[1], base_param_error[1], base_param_bound[1]))
        param_file.write('h_par = Parameter("h", {0}, {1}, {2}, "h")\n'.format(base_param_value[2], base_param_error[2], base_param_bound[2]))
        param_file.write('mnu_par = Parameter("mnu", {0}, {1}, {2}, "\Sigma m_{{\\nu}}")\n'.format(base_param_value[3], base_param_error[3], base_param_bound[3]))
        param_file.write('Nnu_par = Parameter("Nnu",  {0}, {1}, {2}, "N_{{\\rm eff}}")\n'.format(base_param_value[4], base_param_error[4], base_param_bound[4]))
        param_file.write('Ok_par = Parameter("Ok", {0}, {1}, {2}, "\Omega_k")\n'.format(base_param_value[5], base_param_error[5], base_param_bound[5]))
        param_file.write('w_par  = Parameter("w", {0}, {1}, {2}, "w_0")\n'.format(base_param_value[6], base_param_error[6], base_param_bound[6]))
        param_file.write('wa_par = Parameter("wa", {0}, {1}, {2}, "w_a")\n'.format(base_param_value[7], base_param_error[7], base_param_bound[7]))
        param_file.write('wb_par = Parameter("wb", {0}, {1}, {2}, "w_b")\n'.format(base_param_value[8], base_param_error[8], base_param_bound[8]))
        param_file.write('wc_par = Parameter("wc", {0}, {1}, {2}, "w_c")\n'.format(base_param_value[9], base_param_error[9], base_param_bound[9]))
        param_file.write('s8_par = Parameter("s8", {0}, {1}, {2}, "s8")\n'.format(base_param_value[10], base_param_error[10], base_param_bound[10]))
        param_file.write('Pr_par = Parameter("Pr", {0}, {1}, {2}, "c/(H_0r_d)")\n'.format(base_param_value[11], base_param_error[11], base_param_bound[11]))
        param_file.truncate()
    param_file.close()


def write_model_paramDefs(model, model_param_value, model_param_error, model_param_bound):
    """
    Writing the model parameters to the paramDefs.py file.

    Args:
        model_param_value [list]: Value of model parameters
        model_param_error [list]: Error in model parameters
        model_param_bound [list]: Bounds of model parameters
    """
    with open('simplemc\cosmo\paramDefs.py', 'a') as param_file:
        if model == 'Grad_Ok':
            param_file.write('ggama_par = Parameter("ggama", {0}, {1}, {2}, "\gamma")\n'.format(model_param_value[0], model_param_error[0], model_param_bound[0]))
            param_file.write('glambda_par = Parameter("glambda", {0}, {1}, {2}, "\lambda")\n'.format(model_param_value[1], model_param_error[1], model_param_bound[1]))
    param_file.close()

