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

        # get obj rxn ids
        model_obj_rxns = []
        for model in self.models:
            obj_rxn = linear_reaction_coefficients(model).keys()
            model_obj_rxns.extend([rxn.id for rxn in obj_rxn])
        self.objective_rxns = dict(zip(range(self.size), 
                                      model_obj_rxns))

    def run_gifba(self, iters, method, early_stop=True, v=False):
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

        # create variables
        self.create_vars()

        # run iterations
        for iter in range(self.iters):
            if self.v: print("Iteration:", iter)

            # update media for the iteration
            self._is_rerun = False # reset re-run flag for overconsumption
            self._update_media(iter)# maybe change name

            # check early stopping condition
            if self.early_stop:
                if self.v: print("Checking Convergence...")
                delta = self.env_fluxes.loc[iter+1, 0] - self.env_fluxes.loc[iter, 0]
                if np.all(np.abs(delta) < 1e-6):
                    # copy last iter to all future iters
                    self.env_fluxes.loc[(slice(iter+1, None),0), :] = self.env_fluxes.loc[(iter,0), :].values

                    if self.v: print("Converged at iteration", iter)
                    self.iter_converged = iter
                    break

        # drop run col
        self.org_fluxes = self.org_fluxes.droplevel("Run")
        self.env_fluxes = self.env_fluxes.droplevel("Run")

        # return results for total fluxes
        return self.env_fluxes.iloc[-1], self.org_fluxes.groupby(level=["Model"]).sum()

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
        
        # convert to attribute lists
        self.org_exs = list(self.org_exs)
        self.org_rxns = list(self.org_rxns)
        self.exchange_metabolites = list(set(self.exchange_metabolites))
        self.exchanges = list(set(self.exchanges))

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
        

    def _update_media(self, iter ):
        """
        Update the media (f_n,j) for each iteration
        f_{n+1, j} = f_{n,j} + sum(V_{n,i,j})
        """


        # run organism flux function
        self._flux_function(iter)

        # update media: f_n+1 = f_n - sum(v_nij)
        env_tmp = self.env_fluxes.loc[iter, 0][:].to_numpy().reshape(-1, 1)   # (row, col) = (n_ex, 1)     # uptake = positive
        run_exs = self.org_fluxes.loc[:, iter, 0][self.env_fluxes.columns].to_numpy().T # (row, col) = (n_ex, n_org) # uptake = negative flux
        sum_org_flux = run_exs.sum(axis=1).reshape(-1, 1) # (n_ex, n_org) -> (n_ex, 1) sum across orgs

        self.env_fluxes.loc[iter+1, 0] = (env_tmp + sum_org_flux).flatten().round(ROUND) # (n_ex, 1) + (n_ex, 1) -> (n_ex, 1)


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

        # simulate each organism
        for model_idx in range(self.size):
            if self.v: print(" Simulating model:", model_idx+1, " of ", self.size)
            # set media
            self._set_env(iter, model_idx)

            # simulate each org
            self._sim_fba(iter, model_idx)

        # check over consumption
        self._check_overconsumption(iter)

        # once all orgs have been simulated without overconsumption, update internal rxns
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

        return

    def _sim_fba(self, iter, model_idx):
        """General function to run parsimonious FBA (pFBA) on a model and store the results.
        This function runs pFBA on a given model, checks if the solution is above a minimum growth objective,
        and stores the resulting fluxes in the provided DataFrame.
        """
        # run pFBA
        sol1 = self.models[model_idx].slim_optimize()
        if sol1 > GROWTH_MIN_OBJ:
            if self.method == "pfba":
                sol = cb.flux_analysis.parsimonious.pfba(self.models[model_idx])
            elif self.method == "fba":
                sol = self.models[model_idx].optimize()
                
            self.org_fluxes.loc[(model_idx, iter, 0), list(sol.fluxes.index)] = self.rel_abund[model_idx] * sol.fluxes.values * self.step_size
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
            is_overconsumed[env_tmp != 0] = -total_org_flux[np.abs(env_tmp) >= 1e-6].astype(np.longlong) / env_tmp[np.abs(env_tmp) >= 1e-6].astype(np.longlong) # only check non-zero env fluxes

        # check if iteration uses more flux than available in environment
        if is_overconsumed.max().round(ROUND) > 1:
            ex_over = np.argmax(is_overconsumed) # index of flux causing over-consumed

            if self.v: print(self.env_fluxes.columns[ex_over], "over-consumed by factor of", is_overconsumed.max().round(ROUND))

            # adjust only over-consumed bound
            for model_idx in range(self.size):
                self._env_scaling_factors[model_idx, ex_over] = (run_exs[ex_over, model_idx] / (run_exs[ex_over, :].T @ self.rel_abund))
            # re-run flux function with adjusted bounds
            self._is_rerun = True
            self._flux_function(iter)
        

        return
    
    def _update_internal_reactions(self, iter):
        """
        Update internal reactions based on total flux
        """
        for model_idx in range(self.size):
            for rxn in self.models[model_idx].reactions:
                if not(rxn.boundary):
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
