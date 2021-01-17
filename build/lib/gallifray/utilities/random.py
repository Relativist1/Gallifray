import numpy as np
import matplotlib.pyplot as plt


def obs_datapoints_comp(obs_uvdist,obs_amp,obs_err, n=40):
    """Choose n random visibility amplitudes datapoints from the observation and model
        data over same baselines.
    """
    import random
    index = []
    final_new_uv = []
    obs_new_amp = []
    obs_new_sigma = []

    c1 = list(enumerate(zip(obs_uvdist,obs_amp, obs_err)))

    m = random.sample(c1,n)
    for i in np.arange(len(m)):
        index.append(m[i][0])
        uvv1 = m[i][1][0]
        ampp1 = m[i][1][1]
        err1 = m[i][1][2]
        final_new_uv.append(uvv1)
        obs_new_amp.append(ampp1)
        obs_new_sigma.append(err1)
        
    return np.asarray(index), np.asarray(final_new_uv), np.asarray(obs_new_amp), np.asarray(obs_new_sigma)
    
def model_datapoints(index, obs_uvdist, model_amp):
    """Choose n random visibility amplitudes datapoints from the model
        dataset over same baselines.
    """
    
    model_new_amp = []
    c1 = list(enumerate(zip(obs_uvdist,model_amp)))
    for i in index:
        ampp2 = c1[i][1][1]
        model_new_amp.append(ampp2)
    return np.asarray(model_new_amp)
