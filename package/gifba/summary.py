import cobra as cb
import numpy as np
import pandas as pd
from cobra.util.solver import linear_reaction_coefficients
from . import utils
from .config import GROWTH_MIN_OBJ, ROUND

class CommunitySummary:
    """Class to summarize the results of giFBA analysis.
    
    Attributes:
        
    """
    def __init__(self, community, iter_shown=None, element="C"):
        
        # initialize attributes
        self.iter_shown = None
        self.method = None
        self.objective_rxns = None
        self.objective_vals = None
        self.objective_total = None
        self.uptake = None
        self.secretion = None
        self.element = element
        
        self._generate(community, iter_shown)

    def _generate(self, community, iter_shown):
        # check iter_shown is valid
        if self.iter_shown is not None:
            if not isinstance(self.iter_shown, (int, float)) or self.iter_shown < 0 or self.iter_shown >= self.iters:
                raise ValueError("iter_shown must be a non-negative integer less than the number of iterations.")
            else:
                self.iter_shown = int(iter_shown)
        else: 
            self.iter_shown = community.iters - 1

        # pull cumulative organism fluxes
        self.community = community
        self.flux = self.community.org_fluxes.copy()
        self.flux = self.flux.groupby(level='Model').cumsum()
        self.flux = self.flux.xs(self.iter_shown, level='Iteration')

        # extract objectives and create expressions to print
        self.method = self.community.method
        self.objective_rxns = self.community.objective_rxns
        self.objective_vals = [self.flux.loc[model, rxn] for model, rxn in self.objective_rxns.items()]
        self.objective_expressions = [f"1.0 * {rxn} = {self.objective_vals[model]}" for model, rxn in self.objective_rxns.items()]

        # calculate total objective value
        self.objective_total = np.array(self.objective_vals).sum()
        self.objective_total_expression = f"Sum(Model_i Biomass) = {self.objective_total}"

        # create summary dataframe for overall community
        self.env_flux = self.community.env_fluxes.loc[self.iter_shown].copy() - self.community.env_fluxes.loc[0]
        self.env_flux = self.env_flux.T.reset_index()
        self.env_flux.columns = ["Exchange", "Flux"]

        # add metabolite to env_flux
        self.env_flux["Metabolite"] = self.env_flux["Exchange"].map(self.community.ex_to_met)
        self.env_flux = self.env_flux.set_index("Metabolite")


        # add element information
        metabolites = {m.id: m for m in self.community.exchange_metabolites}
        self.env_flux[f"{self.element}-Number"] = [
            metabolites[met_id].elements.get(self.element, 0) if met_id in metabolites else 0
            for met_id in self.env_flux.index
        ]
        self.env_flux[f"{self.element}-Flux"] = self.env_flux[f"{self.element}-Number"] * self.env_flux["Flux"].abs()

        # remove unused fluxes
        self.env_flux = self.env_flux[self.env_flux['Flux'] != 0] # remove zero fluxes

        # create dfs for organisms
        self.flux = self.flux[self.community.org_exs].copy()
        self.flux = self.flux.reset_index()
        self.flux.columns = ["Model"] + list(self.flux.columns[1:])
        self.flux = pd.melt(
            self.flux, 
            id_vars=["Model"], 
            var_name="Exchange", 
            value_name="Flux"
        )
        self.flux["Metabolite"] = self.flux["Exchange"].map(self.community.ex_to_met)
        self.flux["Metabolite"] = self.flux["Metabolite"].fillna(self.flux["Exchange"])
        self.flux = self.flux.set_index(["Model", "Metabolite"])

        # add element information
        metabolites = {m.id if pd.notnull(m.id) else m: m for m in self.community.exchange_metabolites}
        self.flux[f"{self.element}-Number"] = [
            metabolites[met_id].elements.get(self.element, 0) if met_id in metabolites else 0
            for met_id in self.flux.index.get_level_values("Metabolite")
        ]
        self.flux[f"{self.element}-Flux"] = self.flux[f"{self.element}-Number"] * self.flux["Flux"].abs()        

        # remove unused fluxes
        self.flux = self.flux[self.flux['Flux'] != 0] # remove zero fluxes

        return
    
    def to_cytoscape(self):
        # pull pertinent info for cytoscape edge table
        self.cyto_edge = self.flux.reset_index()
        self.cyto_edge["Source"] = self.cyto_edge["Model"].map(self.community.model_names)
        self.cyto_edge["Target"] = self.cyto_edge["Metabolite"]
        self.cyto_edge["Type"] = ["Uptake" if flux < 0 else "Secretion" for flux in self.cyto_edge["Flux"]]
        self.cyto_edge["Value"] = self.cyto_edge["Flux"].abs()

        # drop all other info
        self.cyto_edge = self.cyto_edge[["Source", "Target", "Type", "Value"]]

        # create cytoscape node table
        self.cyto_node = pd.DataFrame()
        self.cyto_node["ID"] = pd.concat([self.cyto_edge["Source"], self.cyto_edge["Target"]]).unique()
        self.cyto_node["Name"] = [self.community.metid_to_name.get(id, id) for id in self.cyto_node["ID"]]
        self.cyto_node["Type"] = ["Organism" if id in self.community.model_names.values() else "Metabolite" for id in self.cyto_node["ID"]]

        return self.cyto_edge, self.cyto_node

    def to_string(self):
        """Display the summary of the community."""
        output = []
        output.append(f"Community Summary (Cumulative through Iteration {self.iter_shown}):\n")
        output.append(f"Optimization Type: {self.method}\n")
        output.append(f"{self.objective_total_expression}\n\n")
        output.append("Uptake:\n")
        uptake = self.env_flux[self.env_flux['Flux'] < 0].copy()
        uptake[f"{self.element}-Flux"] = uptake[f"{self.element}-Flux"] / uptake[f"{self.element}-Flux"].sum()
        uptake["Flux"] = uptake["Flux"].abs()
        output.append(f"{uptake.reset_index().to_string(index=False)}\n\n")
        output.append("Secretion:\n")
        self.secretion = self.env_flux[self.env_flux['Flux'] > 0].copy()
        self.secretion[f"{self.element}-Flux"] = secretion[f"{self.element}-Flux"] / secretion.reset_index()[f"{self.element}-Flux"].sum()
        output.append(f"{secretion.reset_index().to_string(index=False)}\n\n")

        for model in self.flux.index.get_level_values(0).unique():
            output.append("-----------------------------------------------------------------\n")
            output.append(f"{self.community.model_names[model]} (Model {model}) Summary:\n")
            output.append(f"{self.objective_expressions[model]}\n\n")
            output.append(f"{self.community.model_names[model]} Uptake:\n")
            uptake = self.flux.loc[model][self.flux.loc[model]['Flux'] < 0].copy()
            uptake["Flux"] = uptake["Flux"].abs()
            output.append(f"{uptake.reset_index().to_string(index=False)}\n\n")
            output.append(f"Model {model} Secretion:\n")
            secretion = self.flux.loc[model][self.flux.loc[model]['Flux'] > 0].copy()
            output.append(f"{secretion.reset_index().to_string(index=False)}\n\n")
        output.append("\nAccessible at summary.flux or summary.env_flux for cumulative fluxes.\n")

        return "".join(output)

    def __str__(self):
        """Return the string representation of the summary."""
        return self.to_string()
    
    def __repr__(self):
        """Return the string representation of the summary."""
        return self.to_string()

    def _repr_html_(self):
        html = f"<h3>Community Summary (Cumulative through Iteration {self.iter_shown})</h3>"
        html += f"<b>Optimization Type:</b> {self.method}<br>"
        html += f"{self.objective_total_expression}<br>"

        # Community Uptake Table
        html += "<h4>Community Uptake</h4>"
        uptake = self.env_flux[self.env_flux['Flux'] < 0].copy()
        uptake['Flux'] = uptake['Flux'].abs()
        if not uptake.empty:
            html += uptake.reset_index().to_html(index=False)
        else:
            html += "<i>No uptake fluxes</i>"

        # Community Secretion Table
        html += "<h4>Community Secretion</h4>"
        secretion = self.env_flux[self.env_flux['Flux'] > 0].copy()
        if not secretion.empty:
            html += secretion.reset_index().to_html(index=False)
        else:
            html += "<i>No secretion fluxes</i>"

        # Organism-level tables
        for model in self.flux.index.get_level_values(0).unique().sort_values():
            html += f"<hr><h4>{self.community.model_names[model]} (Model {model}) Summary</h4>"
            html += f"{self.objective_expressions[model]}<br>"

            # Organism Uptake Table
            org_uptake = self.flux.loc[model][self.flux.loc[model]['Flux'] < 0].copy()
            org_uptake['Flux'] = org_uptake['Flux'].abs()
            html += f"<b>Model {model} Uptake:</b>"
            if not org_uptake.empty:
                html += org_uptake.reset_index().to_html(index=False)
            else:
                html += "<i>No uptake fluxes</i>"

            # Organism Secretion Table
            org_secretion = self.flux.loc[model][self.flux.loc[model]['Flux'] > 0].copy()
            html += f"<b>Model {model} Secretion:</b>"
            if not org_secretion.empty:
                html += org_secretion.reset_index().to_html(index=False)
            else:
                html += "<i>No secretion fluxes</i>"

        html += "<p>Accessible at summary.flux or summary.env_flux for cumulative fluxes.</p>"

        return html

    

