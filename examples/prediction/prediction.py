import torch
import numpy as np
import pandas as pd
from multiprocessing import Pool
from botorch.models import SingleTaskGP, ModelListGP
from gpytorch.mlls.sum_marginal_log_likelihood import SumMarginalLogLikelihood
from botorch.models.transforms.outcome import Standardize
from botorch.utils.transforms import normalize

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

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
dtype = torch.double
bounds = torch.tensor([[-8,-2,-2,-2,-2,1,1,-5,-5,-5,-5,-3,-1], [-6,5,5,2,2,7,7,-2,-1,-1,-2,-1,1]], device=device, dtype=dtype)

df = pd.read_csv('dataset.csv')

train_x = df.iloc[:,:13]
train_x = torch.Tensor(train_x.values).to(device, dtype)

train_obj = torch.Tensor(df['MWD_re']).to(device, dtype).unsqueeze(-1)
train_con1 = torch.Tensor(df['Mn']).to(device, dtype).unsqueeze(-1)
train_con2 = torch.Tensor(df['conv']).to(device, dtype).unsqueeze(-1)

# using the BO iteration results to train the gaussian models
mll, model = initialize_model(
    train_x, 
    train_obj, 
    train_con1,
    train_con2,
    state_dict=None,
)

# define the combinatorial parameter space for prediction
dataset = [[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13]
          for x1 in np.linspace(-8,-6,2)
           for x2 in np.linspace(-2,5,4)
          for x3 in np.linspace(-2,5,4)
          for x4 in np.linspace(-2,2,9)
          for x5 in np.linspace(-2,2,9)
          for x6 in np.linspace(1,7,4)
          for x7 in np.linspace(1,7,4)
           for x8 in np.linspace(-5,-2,2)
           for x9 in np.linspace(-5,-1,2)
           for x10 in np.linspace(-5,-1,2)
           for x11 in np.linspace(-5,-2,2)
           for x12 in np.linspace(-3,-1,2)
          for x13 in np.linspace(-1,1,2)
          ]

model_obj, model_Mn, model_conv = model.models

class predictor:
    def __call__(self, arg_tuple):
        test_X = arg_tuple
        test_X = torch.Tensor([test_X,]).to(device, dtype)
        test_X = normalize(test_X, bounds)
        MWD = 1 / (model_obj.posterior(test_X).mean)
        Mn = 1000 - model_Mn.posterior(test_X).mean
        conv = 0.6 - model_conv.posterior(test_X).mean
        return MWD[0].item(), Mn[0].item(), conv[0].item()
    
pred = predictor()

# predict the polymerization results using the gaussian model
MWD_list, Mn_list, conv_list = [], [], []
processes = 100
with Pool(processes=processes) as pool:
    output_iter = pool.imap_unordered(pred, dataset)
    for output in output_iter:
        MWD, Mn, conv = output
        MWD_list.append(MWD)
        Mn_list.append(Mn)
        conv_list.append(conv)
        
df = pd.DataFrame(dataset)    
df['MWD'] = MWD_list
df['Mn'] = Mn_list
df['conv'] = conv_list
df.to_csv('predict_results.csv')
