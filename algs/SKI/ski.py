import time 
import torch
import gpytorch
import numpy as np


def gpr(train_x, train_y, test_x, grid_size, sigma2, kern_family, l, double=True):
    train_x = torch.tensor(train_x).double()
    train_y = torch.tensor(train_y).double()
    test_x = torch.tensor(test_x).double()
    grid_size = int(grid_size)
    double = int(double)

    # initialize model for ski
    likelihood = gpytorch.likelihoods.GaussianLikelihood()
    model = gpr_model(train_x, train_y, likelihood, grid_size, kern_family)

    if double:
        model.double()
        likelihood.double()

    # the following two lines appear necessary to avoid errors being raised
    # but we are not training any hyperparameters
    model.train()
    likelihood.train()

    set_hyperparams(model, sigma2, l)

    model.eval()
    likelihood.eval()

    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        tt1 = time.time()
        if double:
            mean_ski = likelihood(model(test_x.double())).mean.numpy()
        else:
            mean_ski = likelihood(model(test_x)).mean.numpy()
        tt2 = time.time()

    return mean_ski, (tt2-tt1)



class gpr_model(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, grid_size, kern_family):
        super(gpr_model, self).__init__(train_x, train_y, likelihood)

        # SKI requires a grid size hyperparameter. This util can help with 
        # that. Here we are using a grid that has the same number of points 
        # as the training data (a ratio of 1.0). Performance can be sensitive 
        # to this parameter, so you may want to adjust it for your own problem 
        # on a validation set.
        ###grid_size = gpytorch.utils.grid.choose_grid_size(train_x, 1.0)

        if len(train_x.size()) > 1:
            [_, dim] = train_x.size()
        else:
            dim = 1

        # choice of kernels: 'squared-exponential', 'matern12', 'matern32'
        if kern_family == 'squared-exponential':
            gpy_kern = gpytorch.kernels.RBFKernel()
        elif kern_family == 'matern12':
            gpy_kern = gpytorch.kernels.MaternKernel(nu=0.5)
        elif kern_family == 'matern32':
            gpy_kern = gpytorch.kernels.MaternKernel(nu=1.5)

        self.mean_module = gpytorch.means.ZeroMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.GridInterpolationKernel(
                gpy_kern, grid_size=grid_size, num_dims=dim
            )
        )

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


def set_hyperparams(model, sigma2, l):
    hypers = {
        'likelihood.noise_covar.noise': torch.tensor(sigma2),
        'covar_module.base_kernel.base_kernel.lengthscale': torch.tensor(l),
        'covar_module.outputscale': torch.tensor(1.0),
    }
    model.initialize(**hypers)
