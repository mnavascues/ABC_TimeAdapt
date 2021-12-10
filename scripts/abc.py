#    TimeAdapt: joint inference of demography and selection
#    Copyright (C) 2021  Miguel de Navascués, Uppsala universitet, INRAE
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import dill
import pyabcranger

project = "test2"
batch = 1
ref_table_sumstats_file = "results/" + project + "/" + str(batch) + "/ref_table_sumstats.pkl"
ref_table_latent_variables_file = "results/" + project + "/" + str(batch) + "/ref_table_latent_variables.pkl"

ref_table_sumstats = dill.load(file = open(ref_table_sumstats_file, "rb"))
ref_table_latent_variables = dill.load(file = open(ref_table_latent_variables_file, "rb"))


ref_table_sumstats.shape
ref_table_sumstats.shape[0] # number of rows
ref_table_latent_variables.shape
ref_table_latent_variables.shape[0] # number of rows

ref_table_sumstats.columns.values # sumstats names
ref_table_latent_variables.columns.values # latent variables names

# for some reason not all simulation have latent variables, this need to be verified
# for the moment subsetting sumstats reftable
ref_table_sumstats = ref_table_sumstats.loc[ref_table_latent_variables.index.values]
ref_table_latent_variables.shape[0] # number of rows

reftable = pyabcranger.reftable(ref_table_latent_variables.shape[0], # nrec: number of rows (records?)
                                [ref_table_latent_variables.shape[0]], # nrecscen: number of recors per scenario
                                [4], # nparam: number of parameters
                                list(ref_table_latent_variables.columns.values), # params_names
                                list(ref_table_sumstats.columns.values), # stats_names
                                ref_table_sumstats.to_numpy(),
                                ref_table_latent_variables.to_numpy(),
                                [1]*ref_table_latent_variables.shape[0])

statobs = ref_table_sumstats.iloc[1]

ref_table_sumstats.columns.values

postres = pyabcranger.estimparam(reftable,list(statobs),"--ntree 500 --parameter Ne_0 --noob 2000 --chosenscen 1",False, False)




