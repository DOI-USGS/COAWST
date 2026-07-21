import matplotlib.pyplot as plt

def get_data_from_file(fname):

    fn = open(fname, "r")
    first_line = fn.readline()
    var_names = list(first_line.split())

    list_of_lists = []
    for line in fn:
        split_line = list(line.strip().split())
        split_line_to_float = [float(e) for e in split_line]
        list_of_lists.append(split_line_to_float)

    transpose_list = list(zip(*list_of_lists))
    output_data = {}
    for i, var_name in enumerate(var_names):
        output_data[var_name] = transpose_list[i]
    return output_data


init_data = get_data_from_file("graupel_sedi_init.txt")
baseline_data01 = get_data_from_file("graupel_sedi_dt_20_runtime1200.txt")
baseline_data01s = get_data_from_file("graupel_sedi_dt_20_semi_sedi_runtime1200.txt")

# plot
fig, ax = plt.subplots(1,1, figsize=(4,6))
ax.plot([i*1000. for i in init_data['mass']], [i for i in init_data['k']], c='k', label='init profile', linestyle='dashed', linewidth=1)
ax.plot([i*1000. for i in baseline_data01['mass']], [i for i in baseline_data01['k']], c='C0', label='dt=20s, integration=1200s', linestyle='solid', linewidth=0.5)
ax.plot([i*1000. for i in baseline_data01s['mass']], [i for i in baseline_data01s['k']], c='C0', label='dt=20s, integration=1200s, semi_sedi=true', linestyle='dashed', linewidth=0.5)
ax.legend()
ax.set_ylabel(r'k level')
ax.set_xlabel(r'q$_{g}$ [g kg$^{-1}$]')
plt.savefig('plot_test_graupel_sedimentation_mpas_59lev.pdf', bbox_inches='tight', dpi=25)
