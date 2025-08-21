import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from scipy import stats

def synthetic_data_generation(num_zirc,eruption,saturation,dist,x_fit_norm,y_fit_norm,sigma_fit_norm,popt,data_dict,filter=False): 
    '''
    Generate synthetic zircon ages based on a chosen distribution and apply uncertainties.
    Parameters:
        num_zirc (int): Number of zircon ages to generate.
        eruption (float): Synthetic eruption age.
        saturation (float): Synthetic saturation age.
        dist (str): Distribution from which the synthetic zircon ages are generated, must be one of the keys in data_dict.
        x_fit_norm (np.ndarray): x-values (residual range) used to evaluate the fitted normal distribution of residuals after the exponential fit.
        y_fit_norm (np.ndarray): Probability density values of the fitted normal distribution on the residuals after the exponential fit at x_fit_norm.
        sigma_fit_norm (float): Standard deviation of the fitted normal distribution of residuals (spread around exponential fit).
        popt (np.ndarray): Array containing the best-fit parameters [a, b] of the exponential function uncer = a * exp(b * age).
        data_dict (dict): Dictionary containing the available distributions.
    
    Returns:
        data_zrc_synthetic (pd.DataFrame): DataFrame containing the synthetic zircon ages, uncertainties, and other parameters.
        '''

    # Extract the distribution from the data_dict and interpolate it
    variable_name = dist
    dist_choice = data_dict[variable_name][1]
    x_old = np.linspace(0, 1, len(dist_choice))
    x_new = np.linspace(0, 1, 10000)
    dist_choice = np.interp(x_new, x_old, dist_choice)

    # Use the distribution as probabilities for generating synthetic ages and translate them to the eruption-saturation range
    probabilities = np.array(dist_choice)
    probabilities = probabilities / np.sum(probabilities)
    lower_bound = 0  
    upper_bound = len(dist_choice) 

    synthetic_data = np.random.choice(np.arange(lower_bound, upper_bound), size=num_zirc, p=probabilities)
    synthetic_data = synthetic_data/len(dist_choice)
    synthetic_ages_initital = (synthetic_data*((saturation)-(eruption)))+(eruption)

    # Generate uncertainties for the ages based on the behaviour of the Torfajökull dataset
    lower_bound = np.abs(x_fit_norm + sigma_fit_norm).argmin()
    upper_bound = np.abs(x_fit_norm - sigma_fit_norm).argmin()
    probabilities = np.array(y_fit_norm[lower_bound:upper_bound])
    probabilities = probabilities / np.sum(probabilities)

    synthetic_norm = np.random.choice(np.arange(lower_bound, upper_bound),size=num_zirc, p=probabilities)
    synthetic_norm = synthetic_norm/len(y_fit_norm)*(max(x_fit_norm)-min(x_fit_norm))+min(x_fit_norm)

    synthetic_uncer = exponential_function(synthetic_ages_initital, *popt)
    synthetic_uncer = synthetic_uncer * synthetic_norm + synthetic_uncer

    # Add Gaussian noise within 2σ to the ages depending on their uncertainties
    noise = np.random.normal(loc=0, scale=synthetic_uncer, size=num_zirc)
    within_bounds = np.logical_and(noise >= -2 * synthetic_uncer, noise <= 2 * synthetic_uncer)
    while not np.all(within_bounds):
        noise[~within_bounds] = np.random.normal(loc=0, scale=synthetic_uncer[~within_bounds], size=np.sum(~within_bounds))
        within_bounds = np.logical_and(noise >= -2 * synthetic_uncer, noise <= 2 * synthetic_uncer)

    synthetic_ages = synthetic_ages_initital + noise

    # Calculate the deltaT, mean_sigma, and deltaT_sigma for the synthetic dataset
    min_age = min(synthetic_ages)
    max_age = max(synthetic_ages)
    deltaT = round(max_age-min_age,3)
    mean_sigma =  round(np.mean(synthetic_uncer),3)
    deltaT_sigma = round(deltaT/mean_sigma,3)

    # Create a DataFrame to store the synthetic zircon ages, uncertainties, and other parameters
    data_zrc_synthetic = pd.DataFrame(columns=['Unit','synthetic_ages','synthetic_uncer','num_zirc','eruption','saturation','deltaT_sigma'])
    data_zrc_synthetic['synthetic_ages'] = synthetic_ages
    data_zrc_synthetic['synthetic_uncer'] = synthetic_uncer
    data_zrc_synthetic['synthetic_ages_inital'] = synthetic_ages_initital
    data_zrc_synthetic['num_zirc'] = num_zirc
    data_zrc_synthetic['eruption'] = eruption
    data_zrc_synthetic['saturation'] = saturation
    data_zrc_synthetic['Unit'] = variable_name
    data_zrc_synthetic['deltaT'] = deltaT
    data_zrc_synthetic['mean_sigma'] = mean_sigma
    data_zrc_synthetic['deltaT_sigma'] = deltaT_sigma
    data_zrc_synthetic['MSWD'] = calculate_mswd(synthetic_ages, synthetic_uncer)
    data_zrc_synthetic = data_zrc_synthetic.sort_values(by='synthetic_ages')
    
    return data_zrc_synthetic

def plot_synthetic_ages(data_zrc_synthetic,data_dict,popt):
    '''
    Plot the synthetic zircon ages, uncertainties, and distributions.
    Parameters:
        data_zrc_synthetic (pd.DataFrame): DataFrame containing the synthetic zircon ages, uncertainties, and other parameters.
        data_dict (dict): Dictionary containing the available distributions and their parameters.
        popt (np.ndarray): Array containing the best-fit parameters [a, b] of the exponential function uncer = a * exp(b * age).
    Returns:
        fig (plt.Figure): Figure containing the plots of synthetic zircon ages, uncertainties, and distributions.
    '''
    # Extract parameters from the synthetic data
    variable_name = data_zrc_synthetic['Unit'].iloc[0]
    num_zirc = len(data_zrc_synthetic['num_zirc'])
    eruption = data_zrc_synthetic['eruption'].iloc[0]
    saturation = data_zrc_synthetic['saturation'].iloc[0]
    deltaT_sigma = data_zrc_synthetic['num_zirc'].iloc[0]
    mean_sigma = data_zrc_synthetic['mean_sigma'].iloc[0]
    deltaT = data_zrc_synthetic['deltaT'].iloc[0]
    mswd = data_zrc_synthetic['MSWD'].iloc[0]

    # Create a figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.delaxes(axes[1, 2])
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(data_dict))) 

    # Plot histogram of the synthetic ages drawn from the choosen distribution
    num_bin = int(0.2*num_zirc)
    axes[0,0].hist(data_zrc_synthetic['synthetic_ages_inital'],density=True, bins=num_bin)
    axes[0,0].set_xlabel('Age ka')
    axes[0,0].set_ylabel('Density')
    axes[0,0].set_title(f'Synthetic Dataset ({variable_name}): \nN = {num_zirc}, age between {eruption}-{saturation} ka')

    # Plot the synthetic uncertainties on top of the exponential fit
    x_fit = np.linspace(0, 170, 100) 
    y_fit = exponential_function(x_fit, *popt)
    axes[0,1].scatter(data_zrc_synthetic['synthetic_ages'],data_zrc_synthetic['synthetic_uncer'])
    axes[0,1].plot(x_fit, y_fit, color='red', label='Fitted Log-Normal Distribution')
    axes[0,1].set_xlabel('Age ka')
    axes[0,1].set_ylabel('sigma ka')
    axes[0,1].set_title(f'Synthetic Dataset ({variable_name}): \nN = {num_zirc}, age between {eruption}-{saturation} ka')

    # Plot a the ranked zircon ages together with their uncertainties, including the initial ages without Gaussian noise
    axes[0, 2].errorbar(
        np.linspace(1, len(data_zrc_synthetic['synthetic_ages']), len(data_zrc_synthetic['synthetic_ages'])), 
        data_zrc_synthetic['synthetic_ages'], 
        data_zrc_synthetic['synthetic_uncer'], 
        fmt='o',
        capsize=3,
        color='black',
        markerfacecolor='black',
        markersize=4,
        ecolor='grey',
        label = 'Ages after applying Gaussian noise',
        zorder = 1
    )
    axes[0,2].scatter(np.linspace(1, len(data_zrc_synthetic['synthetic_ages']), len(data_zrc_synthetic['synthetic_ages'])), data_zrc_synthetic['synthetic_ages_inital'],color = 'red',s = 3, label = 'Ages sampled from distribution',zorder = 2)
    axes[0,2].set_xlabel('N of zircon')
    axes[0,2].set_ylabel('Age with sigma ka')
    axes[0,2].set_title(f'Synthetic Dataset ({variable_name}): \nN = {num_zirc}, age between {eruption}-{saturation} ka')
    axes[0,2].text(0.05, 0.95, f'$\Delta T/\sigma = {deltaT_sigma}$ \n $\sigma$  = {mean_sigma} \n $\Delta T$ = {deltaT} \n MSWD = {mswd:.1f}', transform=axes[0, 2].transAxes,
                    verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    axes[0,2].legend(fontsize=8, loc = 'lower right')

    # Plot the KDE
    eval_points, y_sp = plot_kde(data_zrc_synthetic['synthetic_ages'], data_zrc_synthetic['synthetic_uncer'], False)
    sorted_ages = np.sort(data_zrc_synthetic['synthetic_ages'])
    sorted_ages = normalize_vector(data_zrc_synthetic['synthetic_ages']) #x2
    cumulative_ages = np.linspace(0, 1, len(sorted_ages)) #y2
    eval_points = normalize_vector(eval_points) #x1
    y_sp = y_sp / simpson(y_sp, eval_points) #y1

    for index, (name, (x, y, z)) in enumerate(data_dict.items()):
        axes[1,0].plot(x, y, color=colors[index], label=name)
    axes[1,0].plot(eval_points, y_sp, color='red', linewidth=2, label=f'N={num_zirc}')
    axes[1,0].set_xlabel('Eruption - zircon saturation')
    axes[1,0].set_ylabel('Density')
    axes[1,0].set_title(f'KDE: sampled from {variable_name}')
    axes[1,0].legend(fontsize=8)

    # Plot the CDE
    for index, (name, (x, y, z)) in enumerate(data_dict.items()):
        axes[1,1].plot(x, z, color=colors[index], label=name)
    axes[1,1].step(sorted_ages, cumulative_ages, color='red', linewidth=2, label=f'N={num_zirc}')
    axes[1,1].set_xlabel('Eruption - zircon saturation')
    axes[1,1].set_ylabel('Cumulative Probability')
    axes[1,1].set_title(f'CDE: sampled from {variable_name}')
    axes[1,1].legend(fontsize=8)
    plt.close()
    return fig

def exponential_function(x, a, b):
    return a * np.exp(b * x)

def normalize_vector(vector):
    """
    Normalize a vector to values between 0 and 1.
    """
    min_val = np.min(vector)
    max_val = np.max(vector)
    normalized_vector = (vector - min_val) / (max_val - min_val)
    return normalized_vector

def plot_kde(ages, uncertainty=None, normalize = False):
    """
    Calculate the Kernel Density Estimate (KDE) of ages and return evaluation points and their corresponding density values.
    Parameters:
        ages (array-like): The ages for which to calculate the KDE.
        uncertainty (array-like, optional): The uncertainties associated with the ages. If provided, weights will be applied to the KDE.
        normalize (bool, optional): If True, normalize the ages to the range [0, 1].
    Returns:
        eval_points (numpy.ndarray): The x-values at which the KDE is evaluated.
        y_sp (numpy.ndarray): The density values corresponding to eval_points.
    """
    ages = ages.dropna()
    if normalize is True:
        ages = (ages - np.min(ages)) / (np.max(ages) - np.min(ages))
    if uncertainty is None:
        kde = stats.gaussian_kde(ages)
    else:
        weights = (1/uncertainty**2)/np.sum(1/uncertainty**2)
        kde = stats.gaussian_kde(ages, weights = weights)
    bw = kde.covariance_factor()*np.std(ages) # extract bw
    eval_points = np.linspace(np.min(ages), np.max(ages))
    y_sp = kde.pdf(eval_points)
    return(eval_points, y_sp)

def calculate_mswd(ages, uncertainties):
    """
    Calculate the Mean Squared Weighted Deviation (MSWD) for a set of age determinations.
    Parameters:
        ages (list or numpy array): The age measurements.
        uncertainties (list or numpy array): The uncertainties (standard deviations) of the age measurements.
    Returns:
        float: The calculated MSWD.
    """
    ages = np.array(ages)
    uncertainties = np.array(uncertainties)
    
    N = len(ages)
    weights = 1 / uncertainties**2
    weighted_mean = np.sum(weights * ages) / np.sum(weights)
    mswd = np.sum(weights * (ages - weighted_mean)**2) / (N - 1)
    return mswd