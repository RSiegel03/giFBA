import cobra as cb
import numpy as np
import pandas as pd
from cobra.util.solver import linear_reaction_coefficients
from . import utils
from .config import GROWTH_MIN_OBJ, ROUND
from .summary import CommunitySummary

class gifbaObject:
    """_summary_

    Attributes:
        models: List[cobra.Model], (n_organisms_or_models, )
            A list of cobra.Model objects.
        media: Dict[str: float] 
            The media conditions for the models.
        rel_abund: np.ndarray[float], (n_organisms_or_models, 1)
            The relative abundance of the models, stored as a column vector. 
        id: str, (optional, default=None)
            An optional identifier for the giFBA analysis.
        size: int 
            The number of models in the community. Length of models list.
        objective_rxns: Dict[int: str]
            A dictionary mapping model indices to their objective reaction IDs.
        iters: int
            The number of iterations to run the giFBA analysis.
        method: str (optional, default="pfba")
            The method to use for flux balance analysis ("fba" or "pfba").
        early_stop: bool (optional, default=True)
            A boolean indicating whether to stop early if convergence is reached.
        v: bool (optional, default=False)
            A boolean indicating whether to print verbose output.
        m_vals: List[int], (2, ) (optional, default=[1,1])
            A list containing two integers that define the number of sample points and 
            start points for different runs per iteration. This variable is mainly 
            used for sampling via gifba_sampling, and should remain [1,1] for standard 
            giFBA.
        ex_to_met: Dict[str: str] 
            A dictionary mapping exchange reaction IDs to metabolite IDs.
        metid_to_name: Dict[str: str]
            A dictionary mapping metabolite IDs to their human readable names.
        exchange_metabolites: List[str]
            A list of all unique exchange metabolite IDs across the models. This is 
            a redundant variable, containing all values of ex_to_met.
        exchanges: List[str]
            A list of all unique exchange reaction IDs across the models.
        org_exs: List[str]
            A list of all unique exchange reaction IDs across the models.
        org_rxns: List[str]
            A list of all unique reaction IDs across the models.
        env_fluxes: pd.DataFrame, (n_iterations * m_vals[0] * m_vals[1] + 1, n_exchanges)
            A DataFrame storing the environmental fluxes for each iteration and run. Dataframe
            index is multi-indexed by iteration & run (giFBA drops run index, only necessary
            for sampling). Columns are unique exchange reaction IDs for the entire community.
        org_fluxes: pd.DataFrame, (n_iterations * m_vals[0] * m_vals[1] * n_orgs, n_reactions)
            A DataFrame storing the fluxes of all reactions for each model, iteration, and run. 
            Dataframe index is multi-indexed by model, iteration, and run (giFBA drops run index,
            only necessary for sampling). Columns are unique reaction IDs for the entire community.
        model_names: Dict[int: str]
            A dictionary mapping model indices to their names.
        summary: CommunitySummary
            A CommunitySummary object summarizing the results of the giFBA analysis, see 
            CommunitySummary for more details.

    Methods:
        __init__(self, models, media, rel_abund="equal", id=None):
            Initializes the gifbaObject with the given parameters.
        
        run_gifba(self, iters, method, early_stop=True, v=False):
            Runs the giFBA analysis for a specified number of iterations using the chosen method.
        
        create_vars(self, m_vals=[1,1]):
            Initializes variables for the community giFBA analysis and interpretation. This includes
            setting up DataFrames for environmental and organism fluxes, as well as storing model
            names and reaction mappings.

        update_media(self, iter):
            Updates the media conditions for each iteration based on the fluxes of the models. This 
            method wraps around the _flux_function and handles the media update logic.

        _flux_function(self, iter):
            Runs the flux function for each model in the community for the given iteration. This method
            wraps around the _set_env & _sim_fba methods and handles the overconsumption check.

        _set_env(self, iter, model_idx):
            Sets the exchange reactions of a model to match the environment fluxes for a given iteration and
            model index. This is mainly provided to ensure a cleaner wrapper function.

        _sim_fba(self, iter, model_idx):
            Runs Basic FBA or parsimonious FBA (pFBA) on a model and stores the resulting
            fluxes in the org_fluxes DataFrame. If the model's objective value is below a minimum
             threshold (entailing no growth), the fluxes are not updated (remain zero).

        _check_overconsumption(self, iter):
            Checks for over-consumption of environmental metabolites, scales down overconsumed reactions, and
            re-runs the flux function if necessary.
        
        summarize(self, iter_shown=None):
            Summarizes the results of the giFBA analysis in a CommunitySummary object, formatted to match
            a COBRApy model summary. Also formatting iteration information for a cytoscape-compatible
            node/edge table. 
    """
    def __init__(self, models, media, rel_abund="equal", id=None, step_size=1):
        self.models = utils.check_models(models)
        self.media = media
        self.media = utils.check_media(self)
        self.size = len(self.models)
        self.rel_abund = utils.check_rel_abund(rel_abund, self.size)
        self.id = id
        self.step_size = step_size
        self.iter_converged = None
        self.sim_type = "standard"
        self.flow = None

        # get obj rxn ids
        model_obj_rxns = []
        for model in self.models:
            obj_rxn = linear_reaction_coefficients(model).keys()
            model_obj_rxns.extend([rxn.id for rxn in obj_rxn])
        self.objective_rxns = dict(zip(range(self.size), 
                                      model_obj_rxns))

    def run_additive_gifba(self, iters, method, flow = 0, early_stop=True, v=False):
        """_summary_

        Args:
            iters (_type_): _description_
            method (_type_): _description_
            early_stop (bool, optional): _description_. Defaults to True.
            v (bool, optional): _description_. Defaults to False.

        Returns:
            env_fluxes: _description_
            org_fluxes: _description_
        """
        self.iters = utils.check_iters(iters)
        self.method = utils.check_method(method)
        self.early_stop = early_stop
        self.v = v
        self.flow = flow

        # create variables
        self.create_vars()

        # run iterations
        for iter in range(self.iters):
            self.iter= iter
            print("\nIteration:", iter)

            # update media for the iteration
            self._is_rerun = False # reset re-run flag for overconsumption
            self._update_media(iter)# maybe change name

            # check early stopping condition
            if self.early_stop:
                if self.v: print("Checking Convergence...")
                env_tmp = self.env_fluxes
                delta = env_tmp.loc[iter+1, 0] - env_tmp.loc[iter, 0]
                if np.all(np.abs(delta) < 1e-6):
                    # copy last iter to all future iters
                    self.env_fluxes.loc[(slice(iter+1, None),0), :] = self.env_fluxes.loc[(iter,0), :].values
                    
                    # copy last iter to all future iters
                    vals = self.org_fluxes.loc[(slice(None),iter,0), :].values
                    n_future = self.iters - (iter + 1)
                    if n_future > 0:
                        tiled = np.tile(vals, (n_future, 1))  # shape (n_future * n_models, n_rxns)
                        self.org_fluxes.iloc[-(n_future*self.size):] = tiled


                    if self.v: print("Converged at iteration", iter)
                    self.iter_converged = iter
                    break
        
        if self.iter_converged is None:
            self.iter_converged = self.iters - 1
        print("Total iterations run:", self.iter_converged)

        # drop run col
        self.org_fluxes = self.org_fluxes.droplevel("Run")
        self.env_fluxes = self.env_fluxes.droplevel("Run")
        
        # cumulative sum across iterations
        # self.org_fluxes = self.org_fluxes.groupby(level=["Model"]).cumsum()

        # normalize 
        # self.org_fluxes = self.org_fluxes.apply(
        #     lambda col: col / (1 + self.org_fluxes.index.get_level_values('Iteration') * self.flow),
        #     axis=0)
        
        # account for flow factor (remove flow for all iters after 0)
        # self.env_fluxes.loc[(slice(1, None)), :] -= self.flow * self.env_fluxes.loc[0, :].values
        # self.env_fluxes.loc[(slice(1, None)), :] *= 1/ (1-self.flow)

        # return results for total fluxes
        return self.env_fluxes.iloc[-1], self.org_fluxes.iloc[-self.size:]

    def create_vars(self, m_vals=[1,1]):
        """Initialize variables for giFBA.
        This function sets up the optimization variables for the giFBA analysis.

        """
        self.m_vals = m_vals # default to [1,1] for community giFBA, can be set to [n, m] for sampling via giFBA_sampling m_vals arg
        # get list of all unique rxns and exchanges
        self.ex_to_met = {}
        self.metid_to_name = {}
        self.exchange_metabolites = []
        self.exchanges = []
        self.org_exs = set()
        self.org_rxns = set()
        self.biomass_exs = set()

        # rxns/echanges/boundary mets per model
        for model in self.models:
            exs_set = set(model.exchanges.list_attr("id"))
            self.org_exs = self.org_exs | exs_set # exchanges

            rxns_set = set(model.reactions.list_attr("id"))
            self.org_rxns = self.org_rxns | rxns_set # reactions

            for rxn in model.exchanges:
                mets = list(rxn.metabolites.keys())
                if len(mets) == 1:
                    self.ex_to_met[rxn.id] = mets[0].id if pd.notnull(mets[0].id) else rxn.id
                    self.metid_to_name[mets[0].id] = mets[0].name if pd.notnull(mets[0].name) else mets[0].id
                    self.exchange_metabolites.extend(mets)
                    self.exchanges.append(rxn.id)

                    # add biomass exs to separate set
                    if "biomass" in list(rxn.metabolites.keys())[0].id.lower():
                        self.biomass_exs = self.biomass_exs | {rxn.id}
        
        # convert to attribute lists
        self.org_exs = list(self.org_exs)
        self.org_rxns = list(self.org_rxns)
        self.exchange_metabolites = list(set(self.exchange_metabolites))
        self.exchanges = list(set(self.exchanges))
        self.biomass_exs = list(self.biomass_exs)

        # initialize env
        self.media = utils.check_media(self)
        rows = (self.iters) * self.m_vals[0] * self.m_vals[1] + 1 # add one iteration for initial env
        cols = len(self.org_exs)
        self.env_fluxes = np.zeros((rows, cols))
        env0_masks = [np.array(self.org_exs) == rxn_id for rxn_id in list(self.media.keys())]
        for flux_idx, flux in enumerate(list(self.media.values())):
            self.env_fluxes[0][env0_masks[flux_idx]] = -flux

        #set columns for multi-indexing
        iters_col = np.repeat(np.arange(1, self.iters+1), self.m_vals[0] * self.m_vals[1]) 
        run_col = np.tile(np.arange(self.m_vals[0] * self.m_vals[1]), self.iters)
        iters_col = np.insert(iters_col, 0, 0) # add 0th iteration
        run_col = np.insert(run_col, 0, 0) # add 0th run 
        multi_idx = [iters_col , run_col]
        self.env_fluxes = pd.DataFrame(self.env_fluxes, columns=self.org_exs, index=multi_idx) # convert to interprettable df
        self.env_fluxes.index.names = ["Iteration", "Run"]

        # initialize org_fluxes
        rows = self.iters * self.m_vals[0] * self.m_vals[1] * len(self.models)
        cols = len(self.org_rxns)
        self.org_fluxes = np.zeros((rows, cols)) # pfba will drop run column
        
        # create unique multi-index for 
        models_col = np.tile(np.arange(self.size), self.iters * self.m_vals[0] * self.m_vals[1]) 
        iters_col = np.repeat(np.arange(self.iters), self.m_vals[0] * self.m_vals[1] * self.size) 
        run_col = np.tile(np.repeat(np.arange(self.m_vals[0] * self.m_vals[1]), self.size), self.iters) 
        multi_idx = [models_col, iters_col , run_col]
        self.org_fluxes = pd.DataFrame(self.org_fluxes, columns=self.org_rxns, index=multi_idx)	# convert to interprettable df
        self.org_fluxes.index.names = ["Model", "Iteration", "Run"]


        # store model names
        self.model_names = {model_idx: model.name for model_idx, model in enumerate(self.models)}

        return
        

    def _update_media(self, iter):
        """
        Update the media (f_n,j) for each iteration
        f_{n+1, j} =(1-flow)( f_{n,j} + sum(V_{n,i,j}) ) + flow*(f_{0,j})
        """


        # run organism flux function
        self._flux_function(iter)

        # update media: f_n+1 = f_n - sum(v_nij)
        env_tmp = self.env_fluxes.loc[iter, 0][:].to_numpy().reshape(-1, 1)   # (row, col) = (n_ex, 1)     # uptake = positive
        run_exs = self.org_fluxes.loc[:, iter, 0][self.env_fluxes.columns].to_numpy().T # (row, col) = (n_ex, n_org) # uptake = negative flux
        sum_org_flux = run_exs.sum(axis=1).reshape(-1, 1) # (n_ex, n_org) -> (n_ex, ) sum across orgs

        if self.sim_type == "standard":
            self.env_fluxes.loc[iter+1, 0] = ((1-self.flow) * (env_tmp +  sum_org_flux) + (self.flow * self.env_fluxes.loc[0,0].to_numpy().reshape(-1, 1))).flatten().round(ROUND) # (n_ex, 1) + (n_ex, 1) -> (n_ex, 1)

        elif self.sim_type == "consist_check":
            env_tmp = self.env_fluxes.loc[0, 0][:].to_numpy().reshape(-1, 1)
            run_exs = self.org_fluxes.loc[:, iter, 0][self.env_fluxes.columns].to_numpy().T # (row, col) = (n_ex, n_org) # uptake = negative flux
            run_exs[run_exs < 0] = 0 # only secretion counts
            sum_org_flux = run_exs.sum(axis=1).reshape(-1, 1)
            # sum_org_flux[sum_org_flux < 0] = 0 # no uptake into env
            self.env_fluxes.loc[iter+1, 0] = (env_tmp + sum_org_flux).flatten()#.round(ROUND) # (n_ex, 1) + (n_ex, 1) -> (n_ex, 1)
            biomass_mask = np.isin(self.env_fluxes.columns, self.biomass_exs) # biomass will not be cumulative
            self.env_fluxes.loc[(iter+1, 0), self.env_fluxes.columns[biomass_mask]] = sum_org_flux[biomass_mask].flatten()#.round(ROUND)
        return


    def _flux_function(self, iter):
        """
        run through flux function for organisms
        """
        # # define env bounds per organism for the current iteration
        if not(self._is_rerun): # if first run of iteration, just initialize scaled by rel abund only otherwise do nothing
            self._env_scaling_factors = np.ones((self.size, len(self.org_exs)))  # initialize update rate (used to scale ex flux bounds
            for model_idx in range(self.size):
                self._env_scaling_factors[model_idx, :] = self._env_scaling_factors[model_idx, :] / self.rel_abund[model_idx]

        # if self._is_rerun:
        #     print(self._env_scaling_factors[:, self.ex_over] * self.env_fluxes.loc[iter, 0][self.env_fluxes.columns[self.ex_over]])
        # else:
        #     for ii, ex in enumerate(self.env_fluxes.columns):
        #         for model_idx in range(self.size):
        #             ex_v = self.env_fluxes.loc[iter, 0][ex]
        #             print(self._env_scaling_factors[model_idx, ii] * ex_v, end = ", ")
        #         print("")
        # simulate each organism
        for model_idx in range(self.size):
            # if self.v: print(" Simulating model:", model_idx+1, " of ", self.size)
            # set media
            self._set_env(iter, model_idx)

            # simulate each org
            self._sim_fba(iter, model_idx)

        # check over consumption
        self._check_overconsumption(iter)

        # once all orgs have been simulated without overconsumption, update internal rxns
        if self.sim_type == "standard" and not(self._is_rerun):
            self._update_internal_reactions(iter)

        return

    def _set_env(self, iter, model_idx):
        """
        Function to set the exchange reactions of a model to match the environment fluxes
        for a given iteration and run. This is mainly provided to ensure a cleaner wrapper function.
        """
        for ex in self.models[model_idx].exchanges:
            mask = np.array(self.org_exs) == ex.id
            if mask.any():  # Check if the exchange reaction exists in org_exs
                ex.lower_bound = -self._env_scaling_factors[model_idx, mask] * self.env_fluxes.loc[iter, 0][ex.id]
    
                # print(ex.id, ex.lower_bound)
           

        return

    def _sim_fba(self, iter, model_idx):
        """General function to run parsimonious FBA (pFBA) on a model and store the results.
        This function runs pFBA on a given model, checks if the solution is above a minimum growth objective,
        and stores the resulting fluxes in the provided DataFrame.
        """
        # run pFBA
        sol1 = self.models[model_idx].slim_optimize()
        if model_idx == 0:
            print("#"*45)
            print("Run Info")
        print(f"Objective value (model {model_idx}): {sol1}")
        if sol1 > GROWTH_MIN_OBJ:
            if self.method == "pfba":
                sol = cb.flux_analysis.parsimonious.pfba(self.models[model_idx])
            elif self.method == "fba":
                sol = self.models[model_idx].optimize()
            
            self.org_fluxes.loc[(model_idx, iter, 0), list(sol.fluxes.index)] = self.rel_abund[model_idx] * sol.fluxes.values 
        # do nothing otherwise - already initiated as zeros!
        return
    
    def _check_overconsumption(self, iter):
        """
        Check over-consumption of env. mets. If over-consumption occurs, 
        re-run flux function (recursive subroutine)
        """
        #pull iter info and establish array shapes
        env_tmp = self.env_fluxes.loc[iter, 0][:].to_numpy().reshape(-1, 1)   # (row, col) = (n_ex, 1)     # uptake = positive
        run_exs = self.org_fluxes.loc[:, iter, 0][self.env_fluxes.columns].to_numpy().T # (row, col) = (n_ex, n_org) # uptake = negative flux
        #self.rel_abund  # (n_org, 1)

        # get org fluxes
        total_org_flux = run_exs.sum(axis=1).reshape(-1, 1) # (n_ex, n_org) -> (n_ex, 1) sum across orgs


        # check if environment fluxes are under-saturated
        is_overconsumed = np.zeros_like(total_org_flux)
        with np.errstate(divide='ignore', invalid='ignore'): # ignore division by zero warnings
            is_overconsumed[env_tmp != 0] = -total_org_flux[np.abs(env_tmp) >= 1e-6].astype(np.float64) / env_tmp[np.abs(env_tmp) >= 1e-6].astype(np.float64) # only check non-zero env fluxes

        print("\nenv fluxes (mmol/(gT/hr)):")
        print(self.env_fluxes.loc[iter, 0].T)
        print("\norg fluxes (mmol/(gT/hr)):")
        print(self.org_fluxes.loc[:, iter, 0][self.env_fluxes.columns])
        print("#"*45)
        print()

        if iter == 0 and not self._is_rerun:
            self.X_list = []
            self.OC_list = []
            self.rerun_list = []
            self.iter_list = []




        # check if iteration uses more flux than available in environment
        if not self._is_rerun:
            self.rerun_ct=0
        if self._is_rerun and is_overconsumed.max() != 1:
            ex_over = np.argmax(is_overconsumed) # index of flux causing over-consumed
            if ex_over != self.ex_over:
                self._is_rerun = False # reset re-run flag if different ex is overconsumed on re-run, will trigger new overconsumption loop if necessary
                self.X_list = []
                self.OC_list = []
                self.rerun_list = []
                self.iter_list = []
        if True:
            if is_overconsumed.max() > 1 or (self.rerun_ct !=0 and is_overconsumed.max() <1):
                # print("OC:", is_overconsumed)
                ex_over = np.argmax(is_overconsumed) # index of flux causing over-consumed
                if self._is_rerun and ex_over != self.ex_over:
                    self._is_rerun = False # reset re-run flag if different ex is overconsumed on re-run, will trigger new overconsumption loop if necessary
                    self.X_list = []
                    self.OC_list = []
                    self.rerun_list = []
                    self.iter_list = []


                print(self.env_fluxes.columns[ex_over], f"over-consumed by factor of {is_overconsumed.max():.8f} (rerun count: {self.rerun_ct})")
                print("v"*45)
                # adjust only over-consumed bound
                x_denom = 0
                for model_idx in range(self.size):
                    if self.env_fluxes.columns[ex_over] in self.models[model_idx].reactions:
                        lb_ij = self.models[model_idx].reactions.get_by_id(self.env_fluxes.columns[ex_over]).lower_bound
                        V_ij = run_exs[ex_over, model_idx]
                        a_i = self.rel_abund[model_idx]
                        x_denom += V_ij / lb_ij

                        print("Model idx", model_idx, "    (alpha =", a_i[0],")")
                        print(f"  big V: {V_ij: 3.6f}   mmol/(gT/hr)")
                        print(f"  lil v: {V_ij/a_i[0]: 3.6f}   mmol/(gi/hr)")
                        print(f"     lb: {lb_ij[0]: 3.6f}   mmol/(gi/hr)")
                
                if not self._is_rerun:
                    self.X_list.append(0)
                    self.OC_list.append(0)
                    self.rerun_list.append(-1)
                    self.iter_list.append(iter)

                    # assume n=1 uses this form
                    x_n = env_tmp[ex_over, 0] / x_denom[0]
                    self.X_list.append(x_n)
                if self._is_rerun:
                    # just use one LB for the ex over if we have already re-run
                    for model_idx in range(self.size):
                        if self.env_fluxes.columns[ex_over] in self.models[model_idx].reactions:
                            lb = self.models[0].reactions.get_by_id(self.env_fluxes.columns[ex_over]).lower_bound
                            break
                    self.X_list.append(-lb[0])

                self.OC_list.append(is_overconsumed[ex_over, 0])

                print(f"   X_n:    {self.X_list[-1]:>18.14f}")
                print(f"  OC_n:    {self.OC_list[-1]:>18.14f}")
                print(f" X_n-1:    {self.X_list[-2]:>18.14f}")
                print(f"OC_n-1:    {self.OC_list[-2]:>18.14f}")

                # infer next best X based on deg 1 Newton Method
                m = (self.X_list[-1] - self.X_list[-2]) / (self.OC_list[-1] - self.OC_list[-2])
                b = self.X_list[-1] - m * self.OC_list[-1]
                X_n_p_1 = m * 1 + b  # new env bound at OC = 1
                print(f" X_n+1:    {X_n_p_1:>18.14f}")

                # set new scaling factor for next run \
                # this is div by env bc gets re-multiplied in set_env
                self._env_scaling_factors[:, ex_over] = X_n_p_1 / env_tmp[ex_over, 0]

                
                print("^"*45)
                self.rerun_list.append(self.rerun_ct)
                self.iter_list.append(iter)
                self.ex_over = ex_over
                # self._env_scaling_factors[:, ex_over] = 1 / x_denom
                # self._env_scaling_factors[model_idx, ex_over] = (run_exs[ex_over, model_idx] / (run_exs[ex_over, :].T @ self.rel_abund))
                # re-run flux function with adjusted bounds
                self.rerun_ct = self.rerun_ct + 1
                self._is_rerun = True
                self.OC_n_min_1 = is_overconsumed[ex_over, 0]
                self._flux_function(iter)

        return
    
    def _update_internal_reactions(self, iter):
        """if is_overconsumed.max().round(ROUND) != 1:
            ex_over = np.argmax(is_overconsumed) # index of flux causing over-consumed

            print(self.env_fluxes.columns[ex_over], f"over-consumed by factor of {is_overconsumed.max().round(ROUND):.3f}")

            # adjust only over-consumed bound
            x_denom = 0
            for model_idx in range(self.size):
                lb_ij = self.models[model_idx].reactions.get_by_id(self.env_fluxes.columns[ex_over]).lower_bound
                V_ij = run_exs[ex_over, model_idx]
                a_i = self.rel_abund[model_idx]
                x_denom += V_ij / lb_ij

                print("Model idx", model_idx)
                print("a_i:", a_i)
                print("V_ij:", V_ij, "---- V_ij = v_ij * a_i")
                print("lb_ij:", lb_ij)
            
            X_now = env_tmp[ex_over, 0]/x_denom
            OC_now = is_overconsumed[ex_over, 0]
            print("X_now: ", X_now)
            print("OC_now: ", OC_now)

            if not(self._is_rerun):
                X_old = 0
                self.OC_old = 0
            else:
                X_old = self._env_scaling_factors[0, ex_over] * env_tmp[ex_over, 0]

            print("X_old: ", X_old)
            print("OC_old: ", self.OC_old)

            m = (X_now - X_old) / (OC_now - self.OC_old)
            b = X_now - m * OC_now
            X_new = m * 1 + b  # new env bound at OC = 1

            self._env_scaling_factors[:, ex_over] = X_new / env_tmp[ex_over, 0]

            

            
            print("X_new: ", X_new)
            self.ex_over = ex_over
            # self._env_scaling_factors[:, ex_over] = 1 / x_denom
            # self._env_scaling_factors[model_idx, ex_over] = (run_exs[ex_over, model_idx] / (run_exs[ex_over, :].T @ self.rel_abund))
            # re-run flux function 
        Update internal reactions based on total flux
        """
        for model_idx in range(self.size):
            for rxn in self.models[model_idx].reactions:
                if not(rxn in self.models[model_idx].exchanges):
                    # store previous bounds
                    lb_old = rxn.lower_bound
                    ub_old = rxn.upper_bound
                    
                    # change in flux for given iter
                    last_flux = self.org_fluxes.loc[(model_idx, iter, 0), rxn.id] / self.rel_abund[model_idx, 0]
                    
                    # update bounds
                    rxn.lower_bound = lb_old - last_flux
                    rxn.upper_bound = ub_old - last_flux
        return
    

    def __enter__(self):
        """Context manager entry point."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit point."""
        return False

    def summarize(self, iter_shown=None):
        return CommunitySummary(self, iter_shown)

    def run_gifba(self, iters, method, early_stop=True, v=False):
        """ After each iteration, add only the new fluxes, 
        and do not remove uptaken ones. If fluxes remains 
        the same, update the environment, otherwise- re-do this process """
        self.iters = utils.check_iters(iters)
        self.method = utils.check_method(method)
        self.early_stop = early_stop
        self.v = v

        # create variables
        self.create_vars()

        # run iterations
        for iter in range(self.iters):
            self.iter = iter
            print(f"\nIteration: {iter}")

            # update media for the iteration
            self._is_rerun = False # reset re-run flag for overconsumption
            self.sim_type = "consist_check"
            self._update_media(iter)# maybe change name

            # check early stopping condition
            if self.early_stop and self.iter > 0:
                if self.v: print("Checking Convergence...")
                for per in range(1, iter+1):
                    # check if last (-1) and per+1 iteration from end are the same (accounting for rounding) 
                    env_delta = self.env_fluxes.iloc[iter].values - self.env_fluxes.iloc[iter-per].values
                    org_delta = self.org_fluxes.iloc[self.size*iter:self.size*(iter+1)].values - self.org_fluxes.iloc[self.size*(iter-per):self.size*(iter-per+1)].values
                    if np.all(np.abs(org_delta) < 1e-6) and np.all(np.abs(env_delta) < 1e-6):
                        self.periodicity = per
                        # inside = True # just ensure findign the smallest periodicity if multiple exist
                        self.env_fluxes.loc[(slice(iter+1, None),0), :] = self.env_fluxes.loc[(iter,0), :].values
                        self.iter_converged = iter

                        break
                if self.iter_converged is not None:
                    if self.v: print("Converged at iteration", iter)
                    break
                        

                # deltas = self.env_fluxes.loc[(iter+1, 0), :] - self.env_fluxes.loc[(iter, 0), :]
                # org_flux_tmp = self.org_fluxes.groupby(level="Model").diff().loc[:, iter, 0]
                # if np.all(np.abs(deltas) < 1e-6) and np.all(np.abs(org_flux_tmp) < 1e-6):
                #     # copy last iter to all future iters     
                #     self.env_fluxes.loc[(slice(iter+1, None),0), :] = self.env_fluxes.loc[(iter,0), :].values



                #     if self.v: print("Consistent at iteration", iter)
                #     self.iter_converged = iter
                #     break

        # drop run col
        self.org_fluxes = self.org_fluxes.droplevel("Run")
        self.env_fluxes = self.env_fluxes.droplevel("Run")

        # check periodic/adjust
        env_final, org_final = self.average_periodicity()
        
        # return results for total fluxes
        return env_final, org_final
    
    def average_periodicity(self):
        """Calculate the average periodicity of the system based on the environmental fluxes."""
        # Calculate the difference in environmental fluxes between iterations
        # env_flux_diff = self.env_fluxes.groupby(level="Iteration").diff().iloc[1:]  # Skip the first iteration (initial conditions)

        # if convergence, just return last iteration
        if self.iter_converged is not None:
            periodicity = 1

        # if no convergence, check for periodicity
        else:
            inside = False
            for per in range(1, self.iters+1):
                print(f"Checking periodicity of {per} iterations...")
                # check if last (-1) and per+1 iteration from end are the same (accounting for rounding) 
                env_delta = self.env_fluxes.iloc[-1].values - self.env_fluxes.iloc[-(per+1)].values
                org_delta = self.org_fluxes.iloc[-self.size:].values - self.org_fluxes.iloc[-(per+1)*self.size:-(per)*self.size].values
                if np.all(np.abs(org_delta) < 1e-6) and np.all(np.abs(env_delta) < 1e-6) and not inside:
                    periodicity = per
                    inside = True # just ensure findign the smallest periodicity if multiple exist
                    break
        if periodicity == self.iters:
            print("Model did not converge or show periodicity within the iteration limit, results may be unreliable.")
            print("All iterations will be used for flux calculations, but consider increasing the number of iterations or checking model setup.")
        
        print("Model is periodic and average of the last", periodicity, "iterations will be used for flux calculations.")
        
        # calculate average for the period size
        env_flux_avg = self.env_fluxes.loc[(slice(self.iters - periodicity, self.iters -1)), :].mean()
        org_flux_avg = self.org_fluxes.iloc[-periodicity * self.size:].groupby(level="Model").mean()

        return env_flux_avg, org_flux_avg