import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.stats import norm
from tqdm import tqdm
from math import log, sqrt, pi
from scipy import stats
from scipy.integrate import simpson
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

def normalize_vector(vector):
    """
    Normalizes a vector to the range [0, 1].
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
    eval_points = np.linspace(np.min(ages)-0.05, np.max(ages))
    y_sp = kde.pdf(eval_points)
    return(eval_points, y_sp)

def plot_mcmc_simulation_results(tmin, tmax, unit_of_interest,variable_name):
    """
    Plots the results of the Markov Chain Monte Carlo (MCMC) simulation.
    Parameters:
        tmin (array-like): tested eruption ages from the MCMC simulation.
        tmax (array-like): tested zircon saturation ages from the MCMC simulation.
        unit_of_interest (str): unit from which the ages are derived.
        variable_name (str): The name of the prior distribution used for the Bayesian age calculation.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(range(len(tmin)), tmin, label='Eruption Age')
    plt.plot(range(len(tmax)), tmax, label='Zircon Saturation Age')
    plt.xlabel('Model Step')
    plt.ylabel('Age (ka)')
    plt.legend()
    plt.title(f'Markov Chain Monte Carlo (MCMC) Simulation: {unit_of_interest}, Prior: {variable_name}')
    plt.show()

def plot_histogram(tmin, tmax, burnin, unit_of_interest,variable_name):
    """
    Plots histograms of the Markov Chain Monte Carlo (MCMC) simulations, after an initial cutoff.
    Parameters:
        tmin (array-like): potential eruption ages from the MCMC simulation.
        tmax (array-like): potential zircon saturation ages from the MCMC simulation.
        burnin (int): The number of initial samples to discard as burn-in.
        unit_of_interest (str): unit from which the ages are derived.
        variable_name (str): The name of the prior distribution used for the Bayesian age calculation.
    """
    saturation_age = tmax[burnin:].mean()
    saturation_age_sigma = tmax[burnin:].std()
    eruption_age = tmin[burnin:].mean()
    eruption_age_sigma = tmin[burnin:].std()

    plt.figure(figsize=(10, 6))
    plt.hist(tmin[burnin:], bins=100, density=True, color='blue', alpha=0.7, label='Eruption Age')
    plt.hist(tmax[burnin:], bins=100, density=True, color='red', alpha=0.7, label='Zircon Saturation Age')
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p_max = norm.pdf(x, saturation_age, saturation_age_sigma)
    p_min = norm.pdf(x, eruption_age, eruption_age_sigma)
    plt.plot(x, p_max,  color='red', linewidth=2)
    plt.plot(x, p_min,  color='blue', linewidth=2)
    plt.xlabel('Age (ka)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.title(f'Histogram of Eruption and Zircon Saturation Age: {unit_of_interest}, Prior: {variable_name}')
    plt.show()

def plot_eruption_and_saturation_age(tmin, tmax, burnin, ages, uncer, unit_of_interest,variable_name):
    """
    Plots the eruption age and zircon saturation age, on top of the ranked zircon ages with their uncertainties.
    Parameters:
        tmin (array-like): potential eruption ages from the MCMC simulation.
        tmax (array-like): potential zircon saturation ages from the MCMC simulation.
        ages (array-like): The ages for which to plot the uncertainties.
        uncer (array-like): The uncertainties associated with the ages.
        unit_of_interest (str): unit from which the ages are derived.
        variable_name (str): The name of the prior distribution used for the Bayesian age calculation.
    """
    saturation_age = tmax[burnin:].mean()
    saturation_age_sigma = tmax[burnin:].std()
    eruption_age = tmin[burnin:].mean()
    eruption_age_sigma = tmin[burnin:].std()

    plt.figure(figsize=(10, 6))
    df = pd.concat([pd.Series(ages, name='Ages'), pd.Series(uncer, name='1s')], axis=1)
    df = df.sort_values('Ages', ascending=True)
    # plot zircon ages
    plt.vlines(x=np.arange(len(df)), ymax=df['Ages']+2*df['1s'], ymin=df['Ages']-2*df['1s'])
    # plot and name saturation and eruption age estimate
    plt.hlines(y=eruption_age, xmin=0, xmax=len(ages), linestyle='dashed', color='black')
    plt.annotate(xy=(len(df)/4, eruption_age + 0.05), text='Eruption Age = %1f ka' % eruption_age)
    plt.hlines(y=saturation_age, xmin=0, xmax=len(ages), linestyle='dashed', color='black')
    plt.annotate(xy=(len(df)/4, saturation_age - 0.05), text='Zircon Saturation Age = %1f ka' % saturation_age)
    # plot uncertainty rectangle
    plt.fill_between([0, len(ages)], saturation_age - saturation_age_sigma, saturation_age + saturation_age_sigma, color='0.7', alpha=0.2)
    plt.fill_between([0, len(ages)], eruption_age - eruption_age_sigma, eruption_age + eruption_age_sigma, color='0.2', alpha=0.2)
    plt.title(f'Estimated Eruption and Zircon Saturation Age: {unit_of_interest}, Prior: {variable_name}')
    plt.ylabel('Age (ka)')
    plt.xticks([])

    plt.show()


def calculate_bayesian_age(data_zrc, unit_of_interest,n_chain, dist, model_age,mean_sigma, do_plot = False):
    # Dynamically name the new column
    column_name = f"{model_age}_{dist}"
    
    # Extract data for the specified unit of interest
    ages = data_zrc[data_zrc['Unit'] == unit_of_interest][model_age]
    uncer = data_zrc[data_zrc['Unit'] == unit_of_interest][mean_sigma]

    eval_points, y_sp = plot_kde(ages, uncer, False)
    eval_points = normalize_vector(eval_points) #x1
    y_sp = y_sp / simpson(y_sp, eval_points) #y1

    ages = ages.tolist()
    uncer = uncer.tolist()

    melts_volc = [0.54933, 0.556409, 0.563488, 0.570567, 0.577653, 0.584759, 0.591912, 0.599251, 0.606793, 0.614519, 0.622425, 0.630421, 0.63852, 0.646681, 0.654972, 0.663533, 0.672274, 0.681233, 0.690399, 0.699787, 0.709334, 0.719174, 0.729157, 0.739423, 0.749935, 0.760644, 0.771726, 0.782974, 0.794507, 0.806296, 0.818297, 0.830517, 0.842957, 0.855411, 0.866744, 0.878127, 0.889792, 0.901792, 0.914121, 0.926689, 0.939557, 0.952834, 0.966425, 0.980333, 0.994521, 1.00914, 1.02403, 1.03928, 1.05487, 1.0705, 1.08587, 1.10097, 1.11608, 1.13153, 1.1474, 1.16353, 1.18025, 1.19743, 1.21504, 1.23312, 1.25034, 1.26711, 1.28441, 1.30212, 1.32024, 1.33892, 1.35769, 1.37491, 1.3923, 1.41046, 1.42924, 1.44775, 1.46432, 1.48171, 1.49969, 1.51516, 1.53001, 1.54571, 1.5566, 1.56814, 1.57522, 1.58168, 1.58206, 1.57869, 1.56907, 1.55064, 1.51982, 1.4737, 1.40944, 1.32047, 1.21218, 1.09157, 0.965488, 0.834108, 0.697552, 0.558304, 0.418827, 0.279262, 0.139695, 0.000127237]
    x = np.linspace(0, 1, 100)
    half_norm = norm.pdf(x, loc=0, scale=0.4)
    HalfNormalDistribution= [0.797884560802865, 0.7974280009515494, 0.796068839341563, 0.7938117187969819, 0.7906643181527357, 0.7866373087567604, 0.7817442942234079, 0.7760017339139075, 0.7694288507458882, 0.7620475240542243, 0.753882168338377, 0.7449595988359892, 0.7353088849577002, 0.7249611927031423, 0.7139496172520718, 0.7023090069869697, 0.6900757802537523, 0.6772877362050655, 0.6639838610958632, 0.6502041314134702, 0.6359893152242642, 0.6213807731066313, 0.6064202600153831, 0.5911497293868101, 0.575611140746636, 0.5598462720260221, 0.5438965377242568, 0.5278028139817732, 0.5116052715446525, 0.4953432175127526, 0.4790549466692079, 0.46277760309034227, 0.44654705263311023, 0.43039776679321695, 0.4143627183221232, 0.3984732888863076, 0.38275918894850053, 0.3672483899490719, 0.35196706876734657, 0.3369395643481251, 0.32218834628896875, 0.3077339950994956, 0.2935951937657152, 0.2797887301807786, 0.2663295099388795, 0.25323057893174533, 0.2405031551374514, 0.22815666894929357, 0.21619881135822686, 0.20463558927590253, 0.19347138726644494, 0.1827090349436518, 0.1723498792859793, 0.16239386112415657, 0.1528395950651687, 0.14368445213123268, 0.13492464441275206, 0.1265553110595796, 0.11857060496470585, 0.11096377952814968, 0.10372727492581779, 0.09685280334781428, 0.09033143271260066, 0.08415366840695701, 0.07830953264633381, 0.0727886410954193, 0.06758027643409335, 0.0626734585989266, 0.05805701147461269, 0.05371962585280789, 0.04964991851744419, 0.04583648735539818, 0.042267962429155315, 0.038933052983614144, 0.03582059039222886, 0.03291956707816806, 0.030219171473970206, 0.027708819108233737, 0.025378179930176256, 0.023217202002425913, 0.021216131709211297, 0.019365530641250912, 0.01765628933019398, 0.016079638014541365, 0.014627154625694378, 0.013290770187286376, 0.012062771823391886, 0.01093580357174178, 0.009902865196867894, 0.008957309195320616, 0.00809283618092745, 0.007303488832652165, 0.006583644581141927, 0.005928007202683858, 0.005331597481178707, 0.0047897430900379016, 0.004298067836758348, 0.003852480403458773, 0.003449162707002214, 0.0030845579925822]
    MeltsTZircDistribution= [0.279932598276178, 0.2867031776658, 0.293473757055415, 0.30024433644503, 0.30701715521514, 0.313845271864816, 0.320830884764955, 0.327796643262725, 0.334822491050461, 0.341882681922588, 0.349118200705518, 0.356549119619268, 0.364322271587003, 0.372391900054805, 0.380749356382449, 0.389449671179037, 0.398511598892507, 0.407979196700429, 0.417817958928498, 0.427997270253091, 0.438327785548115, 0.44905544732541, 0.460045602813051, 0.471302277398808, 0.482751168373747, 0.494356584164449, 0.506357798936355, 0.518528309854683, 0.531081757193001, 0.544044974206469, 0.557289930703817, 0.570974188360921, 0.58528593556987, 0.599677772774101, 0.615079934324004, 0.630552743575445, 0.646723719544709, 0.663358513619274, 0.680538233670485, 0.69828385999386, 0.716614802353539, 0.7356555221792, 0.755412283398949, 0.775811959280348, 0.796806321622935, 0.818607007784804, 0.841090957191169, 0.86439730980358, 0.888432981912822, 0.912547254896304, 0.937935447447592, 0.964597559566676, 0.991307020548356, 1.01891846438717, 1.04755288935075, 1.07702403102922, 1.1077018143346, 1.13930634959523, 1.17185978304267, 1.20536150622499, 1.23881668637306, 1.27238055574041, 1.3069433362983, 1.34262162411021, 1.37956742674343, 1.41776712500444, 1.45616189914587, 1.49515411445647, 1.53439000675662, 1.57521575240718, 1.61746239633646, 1.66040197025845, 1.70247028643145, 1.74610988807671, 1.79113931931073, 1.8350940564619, 1.87915814526374, 1.92474422194058, 1.96673182187778, 2.01041078899474, 2.0496007359692, 2.08857186945567, 2.11932089388414, 2.14446347280944, 2.15970102681196, 2.16049597963371, 2.14214951000388, 2.09922776245012, 2.02685847965723, 1.9166107908296, 1.77688486888092, 1.61129844550338, 1.43183209158837, 1.24322922919823, 1.04374357155004, 0.836373869090645, 0.627246359520935, 0.418247843028168, 0.209249326535399, 0.0002508100426302]
    TriangularDistribution= [1.5, 1.4848484848484849, 1.4696969696969697, 1.4545454545454546, 1.4393939393939394, 1.424242424242424, 1.409090909090909, 1.3939393939393938, 1.3787878787878787, 1.3636363636363635, 1.3484848484848484, 1.3333333333333333, 1.3181818181818181, 1.303030303030303, 1.2878787878787878, 1.2727272727272727, 1.2575757575757576, 1.2424242424242424, 1.2272727272727273, 1.2121212121212122, 1.196969696969697, 1.1818181818181819, 1.1666666666666667, 1.1515151515151516, 1.1363636363636365, 1.1212121212121213, 1.1060606060606062, 1.0909090909090908, 1.0757575757575757, 1.0606060606060606, 1.0454545454545454, 1.0303030303030303, 1.0151515151515151, 0.9999999999999999, 0.9848484848484848, 0.9696969696969696, 0.9545454545454546, 0.9393939393939394, 0.9242424242424242, 0.9090909090909091, 0.8939393939393939, 0.8787878787878788, 0.8636363636363636, 0.8484848484848485, 0.8333333333333334, 0.8181818181818181, 0.803030303030303, 0.7878787878787878, 0.7727272727272727, 0.7575757575757576, 0.7424242424242423, 0.7272727272727271, 0.712121212121212, 0.6969696969696969, 0.6818181818181817, 0.6666666666666665, 0.6515151515151514, 0.6363636363636362, 0.6212121212121211, 0.606060606060606, 0.5909090909090908, 0.5757575757575756, 0.5606060606060606, 0.5454545454545452, 0.5303030303030302, 0.5151515151515151, 0.49999999999999983, 0.48484848484848475, 0.4696969696969695, 0.4545454545454544, 0.43939393939393917, 0.4242424242424241, 0.409090909090909, 0.39393939393939376, 0.3787878787878787, 0.3636363636363635, 0.3484848484848484, 0.3333333333333333, 0.3181818181818181, 0.30303030303030304, 0.2878787878787878, 0.2727272727272727, 0.25757575757575746, 0.24242424242424238, 0.2272727272727273, 0.21212121212121204, 0.19696969696969696, 0.1818181818181817, 0.16666666666666663, 0.15151515151515138, 0.1363636363636363, 0.12121212121212119, 0.10606060606060594, 0.09090909090909086, 0.07575757575757561, 0.06060606060606051, 0.04545454545454543, 0.03030303030303018, 0.015151515151515083, 0.0]
    TruncatedNormalDistribution= [0.241971, 0.251742, 0.261481, 0.271153, 0.280724, 0.29016, 0.299423, 0.308478, 0.317288, 0.325817, 0.334031, 0.341892, 0.349368, 0.356425, 0.363032, 0.369157, 0.374774, 0.379856, 0.384378, 0.38832, 0.391662, 0.394389, 0.396487, 0.397946, 0.398759, 0.398922, 0.398434, 0.397297, 0.395518, 0.393104, 0.390067, 0.386423, 0.382188, 0.377383, 0.372031, 0.366156, 0.359787, 0.352951, 0.345681, 0.338008, 0.329966, 0.32159, 0.312916, 0.303978, 0.294815, 0.285461, 0.275953, 0.266327, 0.256617, 0.246858, 0.237083, 0.227324, 0.21761, 0.207972, 0.198437, 0.18903, 0.179775, 0.170694, 0.161808, 0.153134, 0.144689, 0.136486, 0.128538, 0.120856, 0.113448, 0.10632, 0.0994771, 0.092923, 0.0866592, 0.0806857, 0.0750015, 0.069604, 0.0644895, 0.0596534, 0.05509, 0.0507926, 0.0467541, 0.0429665, 0.0394214, 0.0361097, 0.0330223, 0.0301496, 0.0274819, 0.0250094, 0.0227222, 0.0206105, 0.0186646, 0.0168748, 0.0152318, 0.0137263, 0.0123494, 0.0110926, 0.00994734, 0.00890582, 0.00796034, 0.00710363, 0.00632878, 0.00562925, 0.00499887, 0.00443185]
    VolcanicZirconDistribution= [0.545335133570427, 0.552362549285121, 0.55938996499982, 0.566417381408414, 0.573451694478375, 0.58050661396188, 0.587607589519472, 0.594892672217814, 0.602380243936629, 0.610049502743681, 0.617898448430796, 0.625836373705725, 0.633875778399037, 0.641977872085481, 0.650208560617337, 0.658707640374667, 0.667385068616828, 0.676278206907816, 0.685378323765953, 0.694697443892042, 0.704174797556201, 0.71394407676406, 0.723854066251978, 0.734045338222697, 0.744481086182191, 0.755111910752756, 0.766113420680074, 0.777279410925614, 0.788729129991484, 0.800432420550727, 0.812345940043929, 0.824476857182168, 0.836826698333265, 0.849189937637448, 0.860440380738532, 0.871740160686158, 0.883320760619869, 0.895233778569715, 0.907472933652423, 0.919949056056459, 0.932723654268405, 0.945904001861541, 0.95939652020957, 0.973203816081034, 0.987288063694791, 1.00179914992514, 1.01658666090717, 1.03171712350443, 1.04720311991123, 1.06271754003706, 1.07797063027672, 1.09296239013487, 1.10796490309432, 1.12329852032616, 1.13905814837756, 1.15507107892548, 1.17166291592868, 1.18871940548722, 1.20619847793269, 1.22415364593791, 1.24124868413795, 1.25789091287708, 1.27506639070296, 1.29264559362098, 1.3106406234172, 1.32918673600548, 1.34781339259003, 1.36491496000665, 1.38217173159262, 1.4002028490809, 1.4188495344159, 1.43722401954538, 1.45367531631901, 1.47093626533172, 1.48878532213817, 1.50413627708598, 1.51888561033322, 1.53446510415179, 1.54528204105819, 1.55673917794689, 1.56375900729345, 1.57017608482977, 1.57055368507083, 1.56721340549175, 1.55766152364649, 1.53936076666456, 1.50876500985997, 1.46297882113932, 1.39918689258733, 1.31086441368344, 1.20336837717604, 1.08363596175037, 0.958465872031531, 0.828041603755384, 0.692478823982408, 0.554244033653465, 0.415780749218191, 0.277231126612473, 0.138678719074772, 0.0001263115370711]
    VolcanicZirconLowXDistribution= [1.57017608482977, 1.5702447394190535, 1.5703133940083371, 1.5703820485976208, 1.5704507031869044, 1.570519357776188, 1.5702500232909136, 1.5696426997310808, 1.569035376171248, 1.5684280526114154, 1.5678207290515827, 1.56721340549175, 1.5654766997017027, 1.5637399939116554, 1.5620032881216082, 1.560266582331561, 1.5585298765415136, 1.5559978184663144, 1.5526704081059635, 1.5493429977456126, 1.5460155873852617, 1.5426881770249108, 1.53936076666456, 1.5337979017909982, 1.5282350369174362, 1.5226721720438745, 1.5171093071703128, 1.511546442296751, 1.5046026290671837, 1.496277867481611, 1.4879531058960382, 1.4796283443104654, 1.4713035827248926, 1.46297882113932, 1.4513802886753218, 1.4397817562113235, 1.4281832237473255, 1.4165846912833273, 1.404986158819329, 1.39115757632334, 1.37509894379536, 1.35904031126738, 1.3429816787394, 1.3269230462114199, 1.31086441368344, 1.2913196797730035, 1.2717749458625671, 1.2522302119521307, 1.2326854780416945, 1.213140744131258, 1.1924836121373426, 1.170714082059948, 1.1489445519825534, 1.127175021905159, 1.1054054918277643, 1.08363596175037, 1.0608777636196718, 1.0381195654889739, 1.015361367358276, 0.9926031692275777, 0.9698449710968798, 0.9466091203700627, 0.9228956170471271, 0.8991821137241911, 0.8754686104012555, 0.8517551070783197, 0.8280416037553838, 0.8033938256148428, 0.7787460474743015, 0.7540982693337606, 0.7294504911932194, 0.7048027130526784, 0.6799120248615949, 0.6547784266199687, 0.6296448283783428, 0.6045112301367166, 0.5793776318950906, 0.5542440336534648, 0.5290688910288693, 0.5038937484042741, 0.47871860577967873, 0.4535434631550836, 0.4283683205304882, 0.4031853289813072, 0.3779944885075404, 0.3528036480337734, 0.3276128075600066, 0.30242196708623953, 0.27723112661247273, 0.252039779787436, 0.22684843296239962, 0.2016570861373632, 0.1764657393123265, 0.15127439248729008, 0.1260830456622534, 0.10089169883721699, 0.07570035201218059, 0.05050900518714391, 0.025317658362107504, 0.0001263115370711]
    TruncatedLowNormalDistribution = [0.12865812,0.13620658,0.14385367,0.15158655,0.15939159,0.16725439,0.17515984,0.18309212,0.19103478,0.19897077,0.20688249,0.21475189,0.22256048,0.23028941,0.2379196,0.24543171,0.25280634,0.260024,0.26706528,0.27391087,0.28054171,0.28693903,0.29308443,0.29896003,0.3045485,0.30983314,0.31479802,0.319428,0.32370885,0.32762728,0.33117105,0.33432902,0.33709118,0.33944875,0.34139419,0.34292128,0.34402509,0.34470207,0.34495004,0.34476818,0.3441571,0.34311875,0.3416565,0.33977504,0.33748042,0.33477996,0.33168227,0.32819715,0.32433559,0.32010966,0.31553249,0.31061817,0.3053817,0.29983891,0.29400637,0.28790132,0.28154158,0.27494546,0.26813166,0.26111924,0.25392745,0.24657571,0.23908348,0.2314702,0.2237552,0.21595764,0.20809638,0.20018998,0.19225657,0.18431381,0.17637885,0.16846823,0.16059786,0.15278298,0.14503809,0.13737694,0.12981251,0.12235693,0.11502152,0.10781677,0.10075229,0.09383684,0.08707833,0.08048381,0.07405949,0.06781075,0.06174214,0.05585745,0.05015967,0.04465108,0.0393332,0.03420692,0.02927244,0.02452936,0.01997672,0.01561299,0.01143615,0.00744371,0.00363276,0.0]
    NormalDistribution= [0.0, 0.017985079999999987, 0.03705455000000002, 0.057232039999999984, 0.07853747, 0.10098657, 0.12459051999999998, 0.14935543, 0.175282, 0.20236510000000002, 0.23059337000000002, 0.25994889, 0.29040685000000005, 0.32193526, 0.35449471000000005, 0.38803817999999995, 0.42251088, 0.45785018000000005, 0.49398557000000004, 0.53083868, 0.5683234500000001, 0.6063462, 0.64480597, 0.68359481, 0.72259815, 0.76169529, 0.8007599300000001, 0.8396607900000002, 0.87826227, 0.9164251800000001, 0.9540075700000001, 0.99086552, 1.02685407, 1.06182815, 1.0956435100000002, 1.12815769, 1.1592310700000001, 1.18872779, 1.21651678, 1.24247269, 1.26647687, 1.2884182400000002, 1.30819416, 1.32571123, 1.34088601, 1.35364569, 1.3639287, 1.3716851300000001, 1.37687724, 1.37947965, 1.37947965, 1.37687724, 1.3716851300000001, 1.3639287, 1.35364569, 1.34088601, 1.32571123, 1.30819416, 1.2884182400000002, 1.26647687, 1.24247269, 1.21651678, 1.18872779, 1.1592310700000001, 1.12815769, 1.0956435100000002, 1.06182815, 1.02685407, 0.99086552, 0.9540075700000001, 0.9164251800000001, 0.87826227, 0.8396607900000002, 0.8007599300000001, 0.76169529, 0.72259815, 0.68359481, 0.64480597, 0.6063462, 0.5683234500000001, 0.53083868, 0.49398557000000004, 0.45785018000000005, 0.42251088, 0.38803817999999995, 0.35449471000000005, 0.32193526, 0.29040685000000005, 0.25994889, 0.23059337000000002, 0.20236510000000002, 0.175282, 0.14935543, 0.12459051999999998, 0.10098657, 0.07853747, 0.057232039999999984, 0.03705455000000002, 0.017985079999999987, 0.0]
    UniformDistribution = np.linspace(1,1,100)
    Bootstrapped = y_sp.tolist()

    if dist == "melts_volc":
        dist = melts_volc
        variable_name = "metls_volc"
    elif dist == "half_norm":
        dist = half_norm
        variable_name = "half_norm"
    elif dist == "HalfNormalDistribution":
        dist = HalfNormalDistribution
        variable_name = "HalfNormalDistribution"
    elif dist == "MeltsTZircDistribution":
        dist = MeltsTZircDistribution
        variable_name = "MeltsTZircDistribution"
    elif dist == "TriangularDistribution":
        dist = TriangularDistribution
        variable_name = "TriangularDistribution"
    elif dist == "TruncatedNormalDistribution":
        dist = TruncatedNormalDistribution
        variable_name = "TruncatedNormalDistribution"
    elif dist == "VolcanicZirconDistribution":
        dist = VolcanicZirconDistribution
        variable_name = "VolcanicZirconDistribution"
    elif dist == "VolcanicZirconLowXDistribution":
        dist = VolcanicZirconLowXDistribution
        variable_name = "VolcanicZirconLowXDistribution"
    elif dist == "NormalDistribution":
        dist = NormalDistribution
        variable_name = "NormalDistribution"
    elif dist == "TruncatedLowNormalDistribution":
        dist = TruncatedLowNormalDistribution
        variable_name = "TruncatedLowNormalDistribution"
    elif dist == "UniformDistribution":
        dist = UniformDistribution
        variable_name = "UniformDistribution"
    elif dist == "Bootstrapped":
        dist = Bootstrapped
        variable_name = "Bootstrapped"
    else:
        raise ValueError("Invalid distribution specified.")


    # avoid outliers
    # Step 1: Calculate z-scores for ages
    ages_array = np.array(ages)
    z_scores = np.abs((ages_array - np.mean(ages_array)) / np.std(ages_array))

    # Step 2: Define threshold for outliers (e.g., z-score > 3)
    threshold = 2.6

    # Step 3: Remove outliers from ages and corresponding uncertainties
    ages = [ages[i] for i in range(len(ages)) if z_scores[i] < threshold]
    uncer = [uncer[i] for i in range(len(uncer)) if z_scores[i] < threshold]

    
    # Metropolis algorithm for Bayesian estimation
    def dist_ll(dist, mu, sigma, tmin, tmax):
        if tmax < tmin or any(np.isnan(mu)) or any(x <= 0 for x in sigma):
            return float('nan')
        mu = np.array(mu)
        sigma = np.array(sigma)
        sort_indices = np.argsort(mu)
        mu = mu[sort_indices]
        sigma = sigma[sort_indices]
        mu_n, sigma_n = mu[0], sigma[0]
        mu_n1, sigma_n1 = mu[-1], sigma[-1]
        nbins = len(dist)-1
        dt = abs(tmax - tmin)

        loglikelihood = 0.0

        for j in range(len(mu)):
            mu_j, sigma_j = mu[j], sigma[j]

            ix = (mu_j - tmin) / dt * nbins + 1

            if (sigma_j < dt / nbins) and 0 < ix < len(dist)-1:
                f = int(ix)
                d = ix - f
                likelihood = (dist[f + 1] * d + dist[f] * (1 - d)) / dt
            else:
                i_range = range(len(dist))
                likelihood = 0.0
                normconst = 1 / (len(dist) * sigma_j * sqrt(2 * pi))
                for i in range(len(dist)):
                    distx = tmin + dt * (i_range[i] - 1) / nbins
                    likelihood += dist[i] * normconst * np.exp(-(distx - mu_j) ** 2 / (2 * sigma_j * sigma_j))
            
            loglikelihood += np.log(likelihood)

        return loglikelihood

    def metropolis_min(nsteps, dist, mu, sigma):
        tmindist = np.zeros(nsteps)
        tmaxdist = np.zeros(nsteps)
        
        stepfactor = 2.9

        # Sort the dataset from youngest to oldest
        mu = np.array(mu)
        sigma = np.array(sigma)
        sort_indices = np.argsort(mu)
        mu_sorted = mu[sort_indices]
        sigma_sorted = sigma[sort_indices]
        youngest, oldest = mu_sorted[0], mu_sorted[-1]

        dt = oldest - youngest + sigma_sorted[0] + sigma_sorted[-1]
        tmin_step = tmax_step = dt / len(mu)

        tminp = tmin = youngest - sigma_sorted[0]
        tmaxp = tmax = oldest + sigma_sorted[-1]

        llp = ll = dist_ll(dist, mu_sorted, sigma_sorted, tmin, tmax)

        #for i in tqdm(range(nsteps), desc="Running metropolis_min", leave=True):
        for i in range(nsteps):        
            tminp, tmaxp = tmin, tmax
            r = np.random.rand()
            if r < 0.5:
                tmaxp += tmin_step * np.random.randn()
            else:
                tminp += tmax_step * np.random.randn()

            if tminp > tmaxp:
                tminp, tmaxp = tmaxp, tminp

            llp = dist_ll(dist, mu_sorted, sigma_sorted, tminp, tmaxp)
            if np.log(np.random.rand()) < (llp - ll):
                if tminp != tmin:
                    tmin_step = abs(tminp - tmin) * stepfactor
                if tmaxp != tmax:
                    tmax_step = abs(tmaxp - tmax) * stepfactor
                
                ll = llp
                tmin = tminp
                tmax = tmaxp

            tmindist[i] = tmin
            tmaxdist[i] = tmax

        return tmindist, tmaxdist

    nsteps = n_chain  # Length of Markov chain
    burnin = int(n_chain/5)   # Number of steps to discard at beginning of Markov chain
    tmin, tmax = metropolis_min(nsteps=nsteps, dist=dist, mu=ages, sigma=uncer)
    
    saturation_age = tmax[burnin:].mean()
    saturation_age_sigma = tmax[burnin:].std()
    eruption_age = tmin[burnin:].mean()
    eruption_age_sigma = tmin[burnin:].std()

    saturation_age_MAD = abs(tmax[burnin:] - saturation_age).mean()
    eruption_age_MAD = abs(tmin[burnin:] - eruption_age).mean()

    data_zrc.loc[data_zrc['Unit'] == unit_of_interest, f"Bayesian_eruption_{column_name}"] = eruption_age
    data_zrc.loc[data_zrc['Unit'] == unit_of_interest, f"Bayesian_saturation_{column_name}"] = saturation_age
    data_zrc.loc[data_zrc['Unit'] == unit_of_interest, f"Bayesian_eruption_sigma_{column_name}"] = eruption_age_sigma 
    data_zrc.loc[data_zrc['Unit'] == unit_of_interest, f"Bayesian_saturation_sigma_{column_name}"] = saturation_age_sigma 
    data_zrc.loc[data_zrc['Unit'] == unit_of_interest, f"Bayesian_eruption_MAD_{column_name}"] = eruption_age_MAD
    data_zrc.loc[data_zrc['Unit'] == unit_of_interest, f"Bayesian_saturation_MAD_{column_name}"] = saturation_age_MAD 

    # Generate plots
    if do_plot == True:
        plot_mcmc_simulation_results(tmin, tmax, unit_of_interest,variable_name)
        plot_histogram(tmin, tmax, burnin, unit_of_interest,variable_name)
        plot_eruption_and_saturation_age(tmin, tmax, burnin,ages, uncer, unit_of_interest,variable_name)

    return data_zrc