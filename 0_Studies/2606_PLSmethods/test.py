import time

def dual_print(*args, name, **kwargs):
    print(*args)
    with open(name, 'a') as out:
        print(*args, **kwargs, file=out)

def job_report(start_setup, job_time_start, durations, job_dict, report_file_name, job_id, report_id_long):
    
    n_job = len(job_dict)
    
    job_time_stop = time.time()
    delta_t_job = (job_time_stop - job_time_start)/60
    delta_t_sinceStart = (job_time_stop - start_setup)/60
    durations.append(delta_t_job)
    print(f'Job {len(durations)-1} took {delta_t_job:.3f} min.')
    print(f'Since the start it took {delta_t_sinceStart:.3f} min.')
    
    if (n_job == len(durations)):
        dual_print(f'\n >>> Start of Job queue: {report_id}. \n', name=report_file_name)
        dual_print(f'\n >>> End of Job queue: {time.strftime("%d%m%y_%H%M%S", time.localtime())}\n', name=report_file_name)
        dual_print('     <<< ================== JOB REPORT ================== >>>', name=report_file_name)
        dual_print('| ' + f'Job no.'.ljust(8) + ' | ' + f'Job duration'.ljust(12) + ' | ' + f'Job parameter'.ljust(20) + ' | ' + f'Job file'.ljust(14) + ' |', name=report_file_name)
        dual_print('-------------------------------------------------------------------', name=report_file_name)
        for i in range(0,n_job):
            dual_print('| ' + f'Job {i}'.ljust(8) + ' | ' + f'{durations[i]:.2f} min'.ljust(12) + ' | ' + f'{job_dict[str(i)]}'.ljust(20) + ' | ' + f'{job_id}'.ljust(14) + ' |', name=report_file_name)
            dual_print('|-----------------------------------------------------------------|', name=report_file_name)
        
        dual_print(f'\n >>> Start of Job queue: {report_id_long}', name=report_file_name)
        dual_print(f' >>>>> End of Job queue: {time.strftime("%d/%m/%y %H:%M:%S", time.localtime())}\n', name=report_file_name)
        dual_print(f' >>> Total Runtime: {delta_t_sinceStart:.3f} min / {delta_t_sinceStart/60:.3f} h.\n', name=report_file_name)
        dual_print('     <<< =============== END of JOB REPORT =============== >>>', name=report_file_name)
        
    return durations


if __name__ == '__main__':
    report_id = time.strftime("%d%m%y_%H%M%S", time.localtime())
    report_id_long = time.strftime("%d/%m/%y %H:%M:%S", time.localtime())
    report_file_name = f'./simulation/reports/job_report_{report_id}.txt'
    
    start_setup = time.time()
    durations = []
    job_dict = {'0':0,
                '1':1,
                '2':2,
                '3':3,
                '4':4,
                '5':5,
                '6':6,
                '7':7,}

    for i in range(len(job_dict)):
        idnr = time.strftime("%d%m%y_%H%M%S", time.localtime())
        job_time_start = time.time()
        print(job_dict[str(i)])
        time.sleep(0.2)
        durations = job_report(start_setup, job_time_start, durations, job_dict, report_file_name, idnr, report_id_long)
