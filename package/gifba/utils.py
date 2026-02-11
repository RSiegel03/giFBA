import matplotlib.pyplot as plt
import cobra as cb
from cobra.util.solver import linear_reaction_coefficients
from importlib.resources import files
import numpy as np
import pandas as pd


def load_simple_models(number):
	if number == "2b":
		number = "2a"  # same models as 2a
	elif number == "5b" or number == "5c":
		number = "5a"  # same models as 5a
		
	situation_models = {
		"1a": ["sim1a.json"],
		"1b": ["sim1b.json"],
		"1c": ["sim1c_1.json", "sim1c_2.json"],
		"2a": ["sim2a_1.json", "sim2a_2.json"],
		"2c": ["sim2c_1.json", "sim2c_2.json"],
		"2d": ["sim2d_1.json", "sim2d_2.json"],
		"2e": ["sim2e_1.json", "sim2e_2.json"],
		"3a": ["sim3a_1.json", "sim3a_2.json"],
		"3b": ["sim3b_1.json", "sim3b_2.json"],
		"3c": ["sim3c_1.json", "sim3c_2.json"],
		"4a": ["sim4a_1.json", "sim4a_2.json"],
		"5a": ["sim5a_1.json", "sim5a_2.json"],
		"test2": ["test2_sim_1.json", "test2_sim_2.json"],
		"test3": ["test3_sim_1.json", "test3_sim_2.json", "test3_sim_3.json"]

	}

	situation_media = None
	if number in ["1a", "1b", "2a", "2b", "2c", "2d", "2e", "3b", "3c", "5a"]: # A only in media 
		situation_media = {"Ex_A": -10}
	elif number in ["1c"]:
		situation_media = {"Ex_A": -10, "Ex_B": -10}
	elif number in ["3a"]:
		situation_media = {"Ex_A": -10, "Ex_C": -10}
	elif number in ["4a"]:
		situation_media = {"Ex_A": -10, "Ex_B": -10, "Ex_D": -10}
	elif number in ["test2"]:
		situation_media = {"Ex_A": -10, "Ex_D": -10}
	
	models = []
	for file_name in situation_models[number]:
		model_path = files("gifba").joinpath("Simple_Models", file_name)
		models.append(cb.io.load_json_model(str(model_path)))
	
	return models, situation_media

def find_min_medium(community=None, models=None, base_media=None, min_growth=None):
	"""result = {k: max(dict1.get(k, float('-inf')), dict2.get(k, float('-inf')))
          for k in set(dict1) | set(dict2)}"""
	
	if community is not None:
		if isinstance(community.media, (list)):
			base_media = {ex: np.abs(flux) for ex, flux in community.media[0].items()}
			min_growth = community.media[1]
		
		models = community.models
	else:
		models = models.deepcopy()
		base_media = {ex: np.abs(flux) for ex, flux in base_media.items()}
		min_growth = min_growth if min_growth is not None else 0.1

	min_medium = []
	for model in models:
		with model as model_t:
			for rxn_id, uptake in base_media.items():
				if rxn_id in model_t.exchanges:
					met = list(model_t.exchanges.get_by_id(rxn_id).metabolites.keys())[0]
					model_t.add_boundary(met, type="sink", reaction_id=rxn_id+'_tmp',lb=-1*uptake,ub=1000)
				
			for ex in model_t.exchanges:
				ex.lower_bound = -1000
				ex.upper_bound = 1000

			mm = cb.medium.minimal_medium(model_t, min_growth,minimize_components=True)

			model_min_med = mm.to_dict()
			min_medium.append(pd.Series(model_min_med))
	min_medium.append(pd.Series(base_media)) # add base media to ensure all components are included
	
	min_medium = pd.concat(min_medium, axis=1).fillna(0)
	min_medium = (- 1* min_medium.max(axis=1)).to_dict() # convert to uptake and dict

	return min_medium


def check_rel_abund(rel_abund, n_models):
	if rel_abund is None:
		rel_abund = np.ones(n_models) / n_models
	elif isinstance(rel_abund, str):
		rel_abund = np.ones(n_models) / n_models
	elif not isinstance(rel_abund, np.ndarray):
		rel_abund = np.array(rel_abund)
	if rel_abund.ndim != 1:
		rel_abund = rel_abund.flatten()
	if rel_abund.shape[0] != n_models:
		raise ValueError(f"Relative abundances must be a 1D array of length {n_models}.")
	if np.any(rel_abund < 0) or np.sum(rel_abund) == 0:
		raise ValueError("Relative abundances must be non-negative and sum to a positive value.")
	if rel_abund.sum() != 1:
		rel_abund = rel_abund / rel_abund.sum()
		print("Relative abundances set to:", rel_abund)

	rel_abund = rel_abund.astype(float).reshape(-1, 1)
	return rel_abund

def check_iters(iters):
	if iters is None:
		iters = 10
	elif not isinstance(iters, int):
		iters = int(iters)
	if iters < 1:
		iters = 1
		print("Iterations set to:", iters)
	
	return iters

def check_media(community):
	"""None, complete, [min, 0.10], dict"""

	# None or "complete" == Set all exchanges to -1000
	community.media = "complete" if community.media is None else community.media
	if isinstance(community.media, str):
		if community.media.lower() == "complete":
			community.media = dict(zip(community.org_exs, np.full(len(community.org_exs), -1000)))
		else:
			raise ValueError("Media must be None, 'complete', float, or a dict with reaction IDs as keys and flux values as values.")
	
	if isinstance(community.media, (list)):
		community.media = find_min_medium(community)
	elif not isinstance(community.media, (dict, str)):
		raise ValueError("Media must be None, 'complete', float, or a dict with reaction IDs as keys and flux values as values.")

	for rxn_id, flux in community.media.items():
		if not isinstance(rxn_id, str):
			raise ValueError(f"Reaction ID {rxn_id} must be a string.")
		if not isinstance(flux, (int, float)):
			raise ValueError(f"Flux value for reaction {rxn_id} must be a number.")

	return community.media.copy()

def check_models(models):
	if models is None:
		raise ValueError("Models must be provided as a list of cobra.Model objects or single cobra.Model.")

	elif not isinstance(models, (list, cb.Model)):
		raise ValueError("Models must be provided as a list of cobra.Model objects or single cobra.Model.")
	else:
		if isinstance(models, cb.Model):
			models = [models]

	models_list = []
	for model in models:
		if not isinstance(model, cb.Model):
			raise ValueError(f"Model {model} is not a valid cobra.Model object.")
		
		models_list.append(model.copy())
	
	return models_list

def check_method(method):
	if method is None:
		method = "pfba"
	elif not isinstance(method, str):
		raise ValueError("Method must be a string, either 'pfba' or 'fba'.")
	else:
		if isinstance(method, str):
			if method.lower() == "pfba":
				method = "pfba"
			elif method.lower() == "fba":
				method = "fba"
			else:
				raise ValueError("method must be either 'pfba' or 'fba'.")

	return method

