import pandas as pd
import sys
import os

def add_column_postifix(df, postifix):
    dc = {}
    for x in df.columns:
        dc[x] = x+postifix
    return df.rename(columns=dc)

def mapper(x):
    if x == -1:
        return 0
    else:
        return 1

def get_stem(x):
    pos = x.rfind(".")
    name = x.rfind("/")
    if name == -1:
        return x[:pos]
    else:
        return x[name+1:pos]


if len(sys.argv) < 3:
    print("Usage: analyse_data.py [input_file] [output_folder]")
    sys.exit(-1)


data = pd.read_csv(sys.argv[1])
data["speed_ilp_s"] = data["speed_ilp_ms"]/1000
output_folder = sys.argv[2]
base_file = os.path.join(sys.argv[2],get_stem(sys.argv[1]))

data = data.drop(["build_hash","generator"], axis=1)

data_ilp_solved = pd.DataFrame()
data_ilp_solved["solved"] = data["size_ilp"].apply(mapper)
data_ilp_solved["n"] = data["n"]

data_ilp = data.drop(data[data["size_ilp"] == -1].index).drop("size_ilp_feasable",axis=1)
data_ilp_feasable = data.drop(data[data["size_ilp_feasable"] == -1].index).drop("size_ilp",axis=1)

data_ilp["m_blowup_ilp"] = (data_ilp["m"]+data_ilp["size_ilp"])/data_ilp["m"]
data_ilp["m_blowup_blelloch"] = (data_ilp["m"]+data_ilp["size_bl_dp"])/data_ilp["m"]

ilp_grouped_max = data_ilp.copy()

ilp_grouped = data_ilp.groupby("n")
ilp_grouped_mean = ilp_grouped.mean()

ilp_grouped_max['approximation_factor_blelloch'] = ilp_grouped_max['size_bl_dp']/ilp_grouped_max['size_ilp']
ilp_grouped_max['approximation_factor_bl_dp_pert_min_pc'] = ilp_grouped_max['size_bl_dp_pert_min_pc']/ilp_grouped_max['size_ilp']
ilp_grouped_max = ilp_grouped_max.groupby("n").max()

ilp_grouped = pd.concat([ilp_grouped_mean, data_ilp_solved.groupby("n").mean(), 2*add_column_postifix(ilp_grouped.sem(),"_error")],axis=1)

ilp_grouped['approximation_factor_blelloch'] = ilp_grouped['size_bl_dp']/ilp_grouped['size_ilp']
ilp_grouped['approximation_factor_bl_dp_pert_min_pc'] = ilp_grouped['size_bl_dp_pert_min_pc']/ilp_grouped['size_ilp']

ilp_grouped.to_csv(base_file+"-ilp-analysed.csv")
ilp_grouped_max.to_csv(base_file+"-ilp-max-analysed.csv")
data_ilp['approximation_factor_blelloch'] = data_ilp['size_bl_dp']/data_ilp['size_ilp']
data_ilp['approximation_factor_bl_dp_pert_min_pc'] = data_ilp['size_bl_dp_pert_min_pc']/data_ilp['size_ilp']
for x in data_ilp["n"].unique():
    data_ilp[data_ilp["n"] == x].to_csv(base_file+"-ilp-raw-n-"+str(x)+".csv",index=False)

data_ilp.to_csv(base_file+"-ilp-raw-all.csv",index=False)


grp = ilp_grouped["solved"]
grp.to_csv(base_file+"-ilp-analysed-solved.csv")




data_ilp_feasable_solved = pd.DataFrame()
data_ilp_feasable_solved["solved"] = data["size_ilp_feasable"].apply(mapper)
data_ilp_feasable_solved["n"] = data["n"]

data_ilp_feasable["m_blowup_ilp"] = (data_ilp_feasable["m"]+data_ilp_feasable["size_ilp_feasable"])/data_ilp_feasable["m"]
data_ilp_feasable["m_blowup_blelloch"] = (data_ilp_feasable["m"]+data_ilp_feasable["size_bl_dp"])/data_ilp_feasable["m"]
ilp_feasable_grouped = data_ilp_feasable.groupby("n")
ilp_feasable_grouped_mean = ilp_feasable_grouped.mean()


ilp_feasable_grouped = pd.concat([ilp_feasable_grouped_mean, data_ilp_feasable_solved.groupby("n").mean(), 2*add_column_postifix(ilp_feasable_grouped.sem(),"_error")],axis=1)

ilp_feasable_grouped.to_csv(base_file+"-ilp-feasable-analysed.csv")

