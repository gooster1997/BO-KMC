import torch
import Monte
import numpy as np
import pandas as pd
import warnings
from multiprocessing import Pool
from botorch.models import SingleTaskGP, ModelListGP
from gpytorch.mlls.sum_marginal_log_likelihood import SumMarginalLogLikelihood
from botorch.acquisition.objective import ConstrainedMCObjective
from botorch.optim import optimize_acqf
from botorch import fit_gpytorch_model
from botorch.acquisition.monte_carlo import qNoisyExpectedImprovement
from botorch.sampling.samplers import SobolQMCNormalSampler
from botorch.exceptions import BadInitialCandidatesWarning
from botorch.utils.sampling import draw_sobol_samples
from botorch.models.transforms.outcome import Standardize
from botorch.utils.transforms import normalize, unnormalize

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
dtype = torch.double

# set the value ranges of the parameters
bounds = torch.tensor([[-6,-2,-2,-2,-2,1,1,-5,-5,-5,-5,-3,-1], [-4,5,5,2,2,7,7,-2,-1,-1,-2,-1,1]], device=device, dtype=dtype)
standard_bounds = torch.zeros(2, 13, device=device, dtype=dtype)
standard_bounds[1] = 1

# define the batch_size and the hyperparameters
BATCH_SIZE = 5
NUM_RESTARTS = 10
RAW_SAMPLES = 512

def KMC(i):
    ki = 10 ** float(i[0])
    kia = 5.0e3
    kib = 5.0e3
    kpaa = 10 ** float(i[1])
    kpbb = 10 ** float(i[2])
    r1 = 10 ** float(i[3])
    r2 = 10 ** float(i[4])
    kpab = kpaa/r1
    kpba = kpbb/r2
    ktaa = 10 ** float(i[1] + i[5])
    ktbb = 10 ** float(i[2] + i[6])
    ktab = (ktaa * ktbb) ** 0.5
    ktraa = 2 * kpaa * (10 ** float(i[7]))
    ktrab = 2 * kpab * (10 ** float(i[8]))
    ktrba = 2 * kpba * (10 ** float(i[9]))
    ktrbb = 2 * kpbb * (10 ** float(i[10]))
    ratio = 2 * (10 ** float(i[11]))
    conc = 10 ** float(i[12])
    
    xlist = [ki,kia,kib,kpaa,kpab,kpba,kpbb,ktab,ktaa,ktbb,ktraa,ktrab,ktrba,ktrbb,ratio,conc]
    Mn, MWD, conv = Monte.therm_coup(xlist)
    return Mn, MWD, conv
    

def KMC_multi(X):
    X = X.to(torch.device('cpu'))
    Mns, MWD_res, convs = [], [], []
    
    n = len(X)
    p = Pool(n)
    results = []
    for i in range(n):
        r = p.apply_async(KMC, args=(X[i],))
        results.append(r)
    p.close()
    p.join()
    
    for i in results:
        r = i.get()
        Mns.append(r[0])
        MWD_res.append(1/r[1])
        convs.append(r[2])
        
    return torch.Tensor(Mns).to(device, dtype), torch.Tensor(MWD_res).to(device, dtype), torch.Tensor(convs).to(device, dtype)

def generate_initial_data(n=28):
    train_x = draw_sobol_samples(bounds=bounds,n=1, q=n, seed=torch.randint(1000000, (1,)).item()).squeeze(0)
    Mn, MWD_re, conv, r12 = KMC_multi(train_x)
    obj = (MWD_re).unsqueeze(-1)
    con1 = (1000-Mn).unsqueeze(-1)
    con2 = (0.6-conv).unsqueeze(-1)
    return train_x, obj, con1, con2

def initialize_model(train_x, obj, con1, con2, state_dict=None):
    train_x = normalize(train_x, bounds)
    # three gaussian models for MWD, Mn and conversion
    model_obj = SingleTaskGP(train_x, obj, outcome_transform=Standardize(m=1))
    model_con1 = SingleTaskGP(train_x, con1, outcome_transform=Standardize(m=1))
    model_con2 = SingleTaskGP(train_x, con2, outcome_transform=Standardize(m=1))
    
    model = ModelListGP(model_obj, model_con1, model_con2)
    mll = SumMarginalLogLikelihood(model.likelihood, model)
    
    if state_dict is not None:
        model.load_state_dict(state_dict)
    return mll, model

def obj_callable(Z):
    return Z[..., 0]

def constraint1_callable(Z):
    return Z[..., 1]

def constraint2_callable(Z):
    return Z[..., 2]

# define a feasibility-weighted objective for optimization
constrained_obj = ConstrainedMCObjective(
    objective=obj_callable,
    constraints=[constraint1_callable, constraint2_callable],
)

def optimize_acqf_and_get_observation(acq_func):
    """Optimizes the acquisition function, and returns a new candidate and a noisy observation."""
    # optimize
    candidates, _ = optimize_acqf(
        acq_function=acq_func,
        bounds=standard_bounds,
        q=BATCH_SIZE,
        num_restarts=NUM_RESTARTS,
        raw_samples=RAW_SAMPLES,  # used for intialization heuristic
        options={"batch_limit": 5, "maxiter": 200},
    )
    # observe new values 
    new_x = unnormalize(candidates.detach(), bounds=bounds)
    Mn, MWD_re, conv, r12 = KMC_multi(new_x)
    new_obj = (MWD_re).unsqueeze(-1)
    new_con1 = (1000-Mn).unsqueeze(-1)
    new_con2 = (0.6-conv).unsqueeze(-1)
    return new_x, new_obj, new_con1, new_con2

def output(x, obj, con1, con2, filename):
    lis = []
    for i in x:
        lis.append(i.cpu().numpy())
    df = pd.DataFrame(lis)
    df['MWD'] = [1/(i.item()) for i in obj]
    df['Mn'] = [1000 - i.item() for i in con1]
    df['conv'] = [0.6 - i.item() for i in con2]
    df.to_csv(filename, index=False)

warnings.filterwarnings('ignore', category=BadInitialCandidatesWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

N_BATCH = 200
MC_SAMPLES = 128
    
# generate initial training data and initialize model
train_x, train_obj, train_con1, train_con2, train_con3 = generate_initial_data(n=28)
mll, model = initialize_model(train_x, train_obj, train_con1, train_con2, train_con3)

# run N_BATCH rounds of BayesOpt after the initial random batch
for iteration in range(N_BATCH):    

    # fit the models
    fit_gpytorch_model(mll)

    # define the qNEI and qNEI acquisition modules using a QMC sampler
    qmc_sampler = SobolQMCNormalSampler(num_samples=MC_SAMPLES)

    # for best_f, we use the best observed noisy values as an approximation
    qNEI = qNoisyExpectedImprovement(
        model=model, 
        X_baseline=normalize(train_x, bounds),
        sampler=qmc_sampler, 
        objective=constrained_obj,
    )

    # optimize and get new observation
    new_x, new_obj, new_con1, new_con2, new_con3 = optimize_acqf_and_get_observation(qNEI)

    # update training points
    train_x = torch.cat([train_x, new_x])
    train_obj = torch.cat([train_obj, new_obj])
    train_con1 = torch.cat([train_con1, new_con1])
    train_con2 = torch.cat([train_con2, new_con2])

    # reinitialize the models so they are ready for fitting on next iteration
    # use the current state dict to speed up fitting
    mll, model = initialize_model(
        train_x, 
        train_obj, 
        train_con1,
        train_con2,
        model.state_dict(),
    )

    print('iteration %d' % iteration)

# output the iteration results
filename = 'therm_coup.csv'
output(train_x, train_obj, train_con1, train_con2, filename)



