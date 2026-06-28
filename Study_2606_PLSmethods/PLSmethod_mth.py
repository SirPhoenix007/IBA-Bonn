#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import time
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import os
import sys
import json
import tqdm
import warnings
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import numpy as np
import pandas as pd
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import multiprocessing as mup
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
from pybaselines import Baseline
from sklearn.metrics import root_mean_squared_error

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
from colors import load_colors
from functions_pixe import *
from polygauss import *
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times"],
    "text.usetex": True,
    "font.size": 8,
    "pgf.rcfonts": False
})


plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": "\n".join([
          r'\usepackage{amsmath}',
     ]),
})

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
color_schemes = load_colors()

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
def dual_print(*args, name, **kwargs):
    print(*args)
    with open(name, 'a') as out:
        print(*args, **kwargs, file=out)

def job_report(start_setup, job_time_start, durations, job_dict, report_file_name, job_id, job_ids, report_id_long):
    
    n_job = len(job_dict)
    
    job_time_stop = time.time()
    delta_t_job = (job_time_stop - job_time_start)/60
    delta_t_sinceStart = (job_time_stop - start_setup)/60
    durations.append(delta_t_job)
    job_ids.append(job_id)
    print(f'Job {len(durations)-1} took {delta_t_job:.3f} min.')
    print(f'Since the start it took {delta_t_sinceStart:.3f} min.')
    
    if (n_job == len(durations)):
        dual_print('     <<< ================== JOB REPORT ================== >>>', name=report_file_name)
        dual_print('| ' + f'Job no.'.ljust(8) + ' | ' + f'Job duration'.ljust(12) + ' | ' + f'Job parameter'.ljust(20) + ' | ' + f'Job file'.ljust(14) + ' |', name=report_file_name)
        dual_print('-------------------------------------------------------------------', name=report_file_name)
        for i in range(0,n_job):
            dual_print('| ' + f'Job {i}'.ljust(8) + ' | ' + f'{durations[i]:.2f} min'.ljust(12) + ' | ' + f'{job_dict[str(i)]}'.ljust(20) + ' | ' + f'{job_ids[i]}'.ljust(14) + ' |', name=report_file_name)
            dual_print('|-----------------------------------------------------------------|', name=report_file_name)
        
        dual_print(f'\n >>> Start of Job queue: {report_id_long}', name=report_file_name)
        dual_print(f' >>>>> End of Job queue: {time.strftime("%d/%m/%y %H:%M:%S", time.localtime())}\n', name=report_file_name)
        dual_print(f' >>> Total Runtime: {delta_t_sinceStart:.3f} min / {delta_t_sinceStart/60:.3f} h.\n', name=report_file_name)
        dual_print('     <<< =============== END of JOB REPORT =============== >>>', name=report_file_name)
        
    return durations


def rmse(data:list, fit:list, length:int):
    '''
    Root-Mean-Square-Error
    '''
    if data is None or fit is None:
        # A failed baseline fit (e.g. caught LinAlgError) is represented
        # as None upstream; propagate that as a NaN result instead of
        # crashing inside sklearn's input validation.
        return np.nan
    
    data = np.array(data)
    fit = np.array(fit)
    return root_mean_squared_error(data, fit)

def arpls_baseline(bin_data:list, bins:list, lam:int = 1e2):
    bsl_fitter = Baseline(x_data=bins)

    try:
        baseline, params = bsl_fitter.arpls(bin_data, lam=lam, max_iter=500)
    except np.linalg.LinAlgError:
        # Extremely large lam can make the banded Whittaker system
        # numerically non-positive-definite (floating point round-off).
        # Treat this lam as a failed fit rather than crashing the worker.
        return None, None
    subtracted = bin_data - baseline
    return baseline, subtracted

def aspls_baseline(bin_data:list, bins:list, lam:int = 1e2):
    bsl_fitter = Baseline(x_data=bins)

    try:
        baseline, params = bsl_fitter.aspls(bin_data, lam=lam, max_iter=500)
    except np.linalg.LinAlgError:
        return None, None
    subtracted = bin_data - baseline
    return baseline, subtracted

def toy_model(N:int, snr:int, baseline_type:str, plot_flag:bool, D:int=0):
    rng = np.random.default_rng()
    
    x = np.linspace(0,8191,8192)
    # GAUSS
    heights = rng.integers(low=100, high=5e4, size=N)
    widths = rng.integers(low=10, high=100, size=N)
    centers = rng.integers(low=100, high=8100, size=N)
    # print(heights, widths, centers)
    gauss_info = []
    for i in range(0,N):
        gauss_info.append(heights[i])
        gauss_info.append(centers[i])
        gauss_info.append(widths[i])
        
    peak_data = multi_gauss(gauss_info, x)
    P_signal = np.mean(peak_data**2)
    
    # BASELINE
    if (baseline_type == 'sine'):
        if (N > 5):
            sine_amp = rng.uniform(low=10000, high=2000*N)
        else:
            sine_amp = rng.uniform(low=10000, high=10000 + 2000*N)
        sine_freq = rng.uniform(low=0.00025, high=0.0006)
        sine_shift = rng.uniform(low=0, high=10000)
        bsl_param = [sine_amp, sine_freq, sine_shift]
        baseline = sine_amp*np.sin(sine_freq*x + sine_shift)
    elif (baseline_type == 'lin'):
        offset = rng.uniform(low=-1000, high=10000)
        steep = rng.uniform(low=-4, high=4)
        bsl_param = [offset,steep]
        baseline = offset + x*steep
    elif (baseline_type == 'exp'):
        # a*e^{bx+c}
        exp_a = rng.uniform(low=10, high=20*N)
        exp_b = rng.uniform(low=-1e-3, high=-1e-4)
        exp_c = rng.uniform(low=4, high=8)
        bsl_param = [exp_a, exp_b, exp_c]
        baseline = exp_a*np.exp(exp_b*x + exp_c)
    elif (baseline_type == 'poly'):
        degree = D
        poly_param = rng.uniform(low=-0.2, high=0.2, size=D+1)
        poly = np.polynomial.polynomial.Polynomial(poly_param)
        baseline = poly(x/100)
        bsl_param = poly_param
    
       
    # NOISE
    P_noise = P_signal / (10**(snr/10))
    sigma_noise = np.sqrt(P_noise)
    noise = rng.normal(loc=0.0, scale=sigma_noise, size=8192)
    # print(noise)
    # print(f'SNR: {snr:.3f}')
    
    if (plot_flag == True):
        plt.figure(figsize=(20,4), dpi=250)
        plt.plot(x, peak_data, zorder=3, color='black', label='pure data')
        plt.plot(x, peak_data + noise + baseline, zorder=2, color='firebrick', label='full synthetic data')
        # plt.plot(x, peak_data + baseline + 10000)
        plt.plot(x, baseline, color='teal', zorder=1, label='pure baseline')
        plt.grid(axis='y')
        plt.legend()
        plt.show()
    
    toy_model_data = {
        "Bins": x,
        "NumberOfPeaks": N,
        "SNR": snr,
        "GaussianPeaks": gauss_info,
        "Noise": noise,
        "BaselineType": baseline_type,
        "BaselineParameter": bsl_param,
        "Baseline": baseline,
        "SyntheticData": peak_data+noise+baseline       
    }
    return toy_model_data

def evaluate_baseline(toy_model_data:dict):
    
    # plt.figure(figsize=(20,4), dpi=250)
    # plt.grid(axis='y', which='both')
    
    lambda_range = np.linspace(5,13,1601)
    
    bins = toy_model_data['Bins']
    data = toy_model_data['SyntheticData']
    real_bsl = toy_model_data['Baseline']
    lam_min = [0,0]
    minimum = [1e5,1e5]
    res_list = [[],[]]
    for lam in tqdm.tqdm(lambda_range):
        arpls_bsl,_ = arpls_baseline(data, bins, 10**lam)
        if arpls_bsl is None:
            res_list[0].append(np.nan)
        else:
            arpls_result = rmse(real_bsl,arpls_bsl, 8192)
            res_list[0].append(arpls_result)
            if arpls_result < minimum[0]:
                minimum[0] = arpls_result
                lam_min[0] = lam
        
        aspls_bsl,_ = aspls_baseline(data, bins, 10**lam)
        if aspls_bsl is None:
            res_list[1].append(np.nan)
        else:
            aspls_result = rmse(real_bsl,aspls_bsl, 8192)
            res_list[1].append(aspls_result)
            if aspls_result < minimum[1]:
                minimum[1] = aspls_result
                lam_min[1] = lam
            
        # print(f'{lam:.4f}: {result:.4f}')
        # plt.plot(bins,arpls_bsl)
    # plt.plot(bins,real_bsl, color='black')
    # plt.plot(lambda_range,res_list[0],label='arPLS')
    # plt.plot(lambda_range,res_list[1],label='asPLS')
    # plt.legend()
    # plt.ylim(1,20000)
    # plt.yscale('log')
    # print(f'arPLS: Minimum: {lam_min[0]:.2f}: {minimum[0]:.2f} / {(minimum[0]/(data.mean())):.3f}')
    # print(f'asPLS: Minimum: {lam_min[1]:.2f}: {minimum[1]:.2f} / {(minimum[1]/(data.mean())):.3f}')
    
    
    # as_min,_ = aspls_baseline(data, bins, 10**lam_min[1])
    # ar_min,_ = arpls_baseline(data, bins, 10**lam_min[0])
    
    # plt.figure(figsize=(20,4), dpi=250)
    # plt.grid(axis='y', which='both')
    # plt.plot(bins, real_bsl, label='real Baseline')
    # plt.plot(bins, as_min, label='asPLS Baseline')
    # plt.plot(bins, ar_min, label='arPLS Baseline')
    # plt.legend()
    return {
        "arpls_lambda_min": lam_min[0],
        "arpls_rmse_min": minimum[0],
        "arpls_rmse_all": res_list[0],
        "aspls_lambda_min": lam_min[1],
        "aspls_rmse_min": minimum[1],
        "aspls_rmse_all": res_list[1],}   
    
    
def single_experiment(args):
    N, snr, baseline_type, D = args
    
    model = toy_model(N=N, snr=snr, baseline_type=baseline_type, plot_flag=False, D=D)
    
    result = evaluate_baseline(model)
    
    # model["Bins"] = model["Bins"].tolist()
    # model["Noise"] = model["Noise"].tolist()
    # model["Baseline"] = model["Baseline"].tolist()
    # model["SyntheticData"] = model["SyntheticData"].tolist()
    
    # with open(f"./simulation/toy_model_0{i}.json", "w") as tm:
    #     json.dump(model, tm, indent=2)
    # with open(f"./simulation/pls_results_toy_model_0{i}.json", "w") as rtm:
    #     json.dump(result, rtm, indent=2)
        
    return model, result

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        return super().default(obj)


if __name__ == "__main__":
    
    report_id = time.strftime("%d%m%y_%H%M%S", time.localtime())
    report_id_long = time.strftime("%d/%m/%y %H:%M:%S", time.localtime())
    report_file_name = f'./simulation/reports/job_report_{report_id}.txt'
    
    start_setup = time.time()
    durations, jobs_ids = [],[]
    
    jobs_dict = {
        '0':[8,30,'sine',0],
        '1':[24,30,'sine',0],
        '2':[40,30,'sine',0],
        '3':[56,30,'sine',0],
        '4':[8,20,'sine',0],
        '5':[24,20,'sine',0],
        '6':[40,20,'sine',0],
        '7':[56,20,'sine',0],
        '8':[8,15,'sine',0],
        '9':[24,15,'sine',0],
        '10':[40,15,'sine',0],
        '11':[56,15,'sine',0],
        '12':[8,30,'exp',0],
        '13':[24,30,'exp',0],
        '14':[40,30,'exp',0],
        '15':[56,30,'exp',0],
        '16':[8,20,'exp',0],
        '17':[24,20,'exp',0],
        '18':[40,20,'exp',0],
        '19':[56,20,'exp',0],
        '20':[8,15,'exp',0],
        '21':[24,15,'exp',0],
        '22':[40,15,'exp',0],
        '23':[56,15,'exp',0],
        '24':[8,30,'lin',0],
        '25':[24,30,'lin',0],
        '26':[40,30,'lin',0],
        '27':[56,30,'lin',0],
        '28':[8,20,'lin',0],
        '29':[24,20,'lin',0],
        '30':[40,20,'lin',0],
        '31':[56,20,'lin',0],
        '32':[8,15,'lin',0],
        '33':[24,15,'lin',0],
        '34':[40,15,'lin',0],
        '35':[56,15,'lin',0]
    }
    
    print(f'Queue started with {len(jobs_dict)} jobs to run.')
    print(f'Jobs in queue:')
    for i in range(0,len(jobs_dict)):
        print(i, jobs_dict[str(i)])
        
        
    for i in range(0,len(jobs_dict)):
        job_time_start = time.time()
        print(f'Starting on job {i} with the parameters {jobs_dict[str(i)]}')
        idnr = time.strftime("%d%m%y_%H%M%S", time.localtime())
        
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            
        jobs = [jobs_dict[str(i)] for _ in range(100)]
        results = []
        with ProcessPoolExecutor(max_workers=11, mp_context=mup.get_context("spawn")) as executor:
            futures = [executor.submit(single_experiment, job) for job in jobs]

            for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
                results.append(future.result())
        
        warning_count = len(caught)
        
        print(f'Completed job {i} with {len(results)} experiments!')
        
        print(f'Supression of {warning_count} warnings!')
        
        with open(f"./simulation/toy_model_run_{idnr}.json", "w") as dumptruck:
            json.dump(results, dumptruck, cls=NumpyEncoder, indent=2)
        
        durations = job_report(start_setup, job_time_start, durations, jobs_dict, report_file_name, idnr, jobs_ids , report_id_long)

        