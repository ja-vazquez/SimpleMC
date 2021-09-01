def write_paramDefs(param_value_settings, param_error_settings, param_bound_settings):
    with open('simplemc\cosmo\paramDefs.py', 'r+') as param_file:
        data = param_file.read()
        param_file.seek(0)
        param_file.write('from simplemc.cosmo.Parameter import Parameter\n\n')
        param_file.write('Om_par = Parameter("Om", {0}, {1}, {2}, "\Omega_m")\n'.format(param_value_settings[0], param_error_settings[0], param_bound_settings[0]))
        param_file.write('Obh2_par = Parameter("Obh2", {0}, {1}, {2}, "\Omega_{{b}}h^2")\n'.format(param_value_settings[1], param_error_settings[1], param_bound_settings[1]))
        param_file.write('h_par = Parameter("h", {0}, {1}, {2}, "h")\n'.format(param_value_settings[2], param_error_settings[2], param_bound_settings[2]))
        param_file.write('mnu_par = Parameter("mnu", {0}, {1}, {2}, "\Sigma m_{{\\nu}}")\n'.format(param_value_settings[3], param_error_settings[3], param_bound_settings[3]))
        param_file.write('Nnu_par = Parameter("Nnu",  {0}, {1}, {2}, "N_{{\\rm eff}}")\n'.format(param_value_settings[4], param_error_settings[4], param_bound_settings[4]))
        param_file.write('Ok_par = Parameter("Ok", {0}, {1}, {2}, "\Omega_k")\n'.format(param_value_settings[5], param_error_settings[5], param_bound_settings[5]))
        param_file.write('w_par  = Parameter("w", {0}, {1}, {2}, "w_0")\n'.format(param_value_settings[6], param_error_settings[6], param_bound_settings[6]))
        param_file.write('wa_par = Parameter("wa", {0}, {1}, {2}, "w_a")\n'.format(param_value_settings[7], param_error_settings[7], param_bound_settings[7]))
        param_file.write('wb_par = Parameter("wb", {0}, {1}, {2}, "w_b")\n'.format(param_value_settings[8], param_error_settings[8], param_bound_settings[8]))
        param_file.write('wc_par = Parameter("wc", {0}, {1}, {2}, "w_c")\n'.format(param_value_settings[9], param_error_settings[9], param_bound_settings[9]))
        param_file.write('s8_par = Parameter("s8", {0}, {1}, {2}, "s8")\n'.format(param_value_settings[10], param_error_settings[10], param_bound_settings[10]))
        param_file.write('Pr_par = Parameter("Pr", {0}, {1}, {2}, "c/(H_0r_d)")\n'.format(param_value_settings[11], param_error_settings[11], param_bound_settings[11]))
        param_file.truncate()
    param_file.close()


