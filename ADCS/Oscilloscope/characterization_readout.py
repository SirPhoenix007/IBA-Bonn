import pyvisa
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
import struct
import math
import warnings
import gc
import time
from datetime import datetime
import numpy as np
import seaborn as sns
import scipy
from scipy.optimize import brentq
import threading
import queue
import os
import h5py


ch_colors = ['green', 'purple', 'blue']
ch_fit_colors = ['orange', 'red', 'cyan']

# Define the error functions --- taken from Niclas code
def errorfunc(x, a, b, c, d):
    x = np.array(x)
    return a*np.array(scipy.special.erf((x-b)*c))+d

def gaussian(x,A,mu,sigma):
    return A * np.exp(-(x-mu)**2 / (2* sigma**2))

# Define the difference of the error functions --- taken from Niclas code
def diff_error(x, a, b, c, d, e, f, g, h):
    x = np.array(x)
    return errorfunc(x, a, b, c, d) + errorfunc(x, e, f, g, h)

def onoff_errorfunc(x, A, t_on, a_on, t_off, a_off, C):
    '''
    1/2 because erf from -1 to 1, A: amplitude above base
    C: baseline (const),    t_on: center time of rising edge, t_off: center time of falling edge
                            a_on: steepness of rising edge, a_off: steepness of falling edge
    '''
    # x = np.array(x)
    return C + 0.5 * A * (
        scipy.special.erf((x - t_on) * a_on)
        - scipy.special.erf((x - t_off) * a_off)
    )

def amplifier_func(t, A, tau_rise, tau_decay, t0, C_amp):
    '''
    A: Amplitude (negative here), t: x-values
    tau_rise: rise const, tau_decay: decay const.
    t0: time from which the function rises
    '''
    t = np.array(t)
    t_pulse = t-t0
    volt = np.full_like(t,C_amp,dtype=float)
    mask = t_pulse >= 0
    volt[mask] = A* (np.exp(-t_pulse[mask]/tau_decay) - np.exp(-t_pulse[mask]/tau_rise)) + C_amp
    return volt

def csa_large_signal(t, A, tau2, B, t1, C):
    '''
    from the papers on CSA in smaller CMOS technology:
        They state, that the two exponentials only work for small signals, but
        in saturation, the decay is linear (-> Niclas' approach). 
        Although the plots in the paper do not look like it (but look like our signal) 
        it is worth the try. If not, we go on with the comparator function.
    '''
    t = np.array(t)
    volt = np.full_like(t, C, dtype=float)

    t2 = t1 + (2.0 * A / B)
    
    # pulse window
    mask = (t >= t1) & (t <= t2)
    t_pulse = t[mask] - t1
    
    volt[mask] = (A * (1 - np.exp(-t_pulse / tau2)) - 0.5 * B * t_pulse) + C
    
    return volt

fit_param_names = {
    'comp': ["A_comp", "t_on", "a_on", "t_off", "a_off", "C_comp"],
    'amp': ["A_amp", "t_on", "a_on", "t_off", "a_off", "C_amp"],
    # 'amp': ["A_amp", "tau", "B"],  # from the function from the papers
    'inj': ['A_inj', 't0_inj', 'rise_inj', 'C_inj']
}


"""Modify the following global variables according to the model""" # --- taken from Niclas code
SDS_RSC = 'TCPIP::131.220.221.73::inst0::INSTR' # this includes the IP adress of the scope
CHANNEL = ["C1","C2","C3","C4"] # available channels
HORI_NUM = 10 #default 10
vertical_divisions = 8
tdiv_enum = [100e-12, 200e-12, 500e-12, 1e-9,
2e-9, 5e-9, 10e-9, 20e-9, 50e-9, 100e-9, 200e-9, 500e-9, \
1e-6, 2e-6, 5e-6, 10e-6, 20e-6, 50e-6, 100e-6, 200e-6, 500e-6, \
1e-3, 2e-3, 5e-3, 10e-3, 20e-3, 50e-3, 100e-3, 200e-3, 500e-3, \
1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]#remove 100e-12, for slow osci and add for fast --- taken from Niclas code


"""The following code realizes the process of waveform reconstruction with slice""" # --- taken from Niclas code


def main_desc(recv):

#receiving raw data from osci
    WAVE_ARRAY_1 = recv[0x3c:0x3f + 1]
    wave_array_count = recv[0x74:0x77 + 1]
    first_point = recv[0x84:0x87 + 1]
    sp = recv[0x88:0x8b + 1]
    v_scale = recv[0x9c:0x9f + 1]
    v_offset = recv[0xa0:0xa3 + 1]
    interval = recv[0xb0:0xb3 + 1]
    code_per_div = recv[0xa4:0Xa7 + 1]#
    adc_bit = recv[0xac:0Xad + 1]
    delay = recv[0xb4:0xbb + 1]
    tdiv = recv[0x144:0x145 + 1]
    probe = recv[0x148:0x14b + 1]

#converting raw data into values
    data_bytes = struct.unpack('i', WAVE_ARRAY_1)[0]
    point_num = struct.unpack('i', wave_array_count)[0]
    fp = struct.unpack('i', first_point)[0]
    sp = struct.unpack('i', sp)[0]
    interval = struct.unpack('f', interval)[0]
    delay = struct.unpack('d', delay)[0]
    tdiv_index = struct.unpack('h', tdiv)[0]
    probe = struct.unpack('f', probe)[0]
    vdiv = struct.unpack('f', v_scale)[0] * probe
    offset = struct.unpack('f', v_offset)[0] * probe
    code = struct.unpack('f', code_per_div)[0]
    adc_bit = struct.unpack('h', adc_bit)[0]
    tdiv = tdiv_enum[tdiv_index]

    return vdiv, offset, interval, delay, tdiv, code, adc_bit

def main_wf_data(number_ch):

#setting up to read waveform data
    sds.write(":WAVeform:STARt 0")
    sds.write("WAV:SOUR {}".format(CHANNEL[number_ch]))
    sds.write("WAV:PREamble?")
    recv_all = sds.read_raw()
    recv = recv_all[recv_all.find(b'#') + 11:]
    vdiv, ofst, interval, trdl, tdiv, vcode_per, adc_bit = main_desc(recv)
    points = sds.query(":ACQuire:POINts?").strip()
    points = float(sds.query(":ACQuire:POINts?").strip())
    one_piece_num = float(sds.query(":WAVeform:MAXPoint?").strip())
    if points > one_piece_num:
        sds.write(":WAVeform:POINt {}".format(one_piece_num))
    if adc_bit > 8:
        sds.write(":WAVeform:WIDTh WORD")
    read_times = math.ceil(points / one_piece_num)
    recv_all = []

#reading raw waveform data for a single channel
    for i in range(0, read_times):
        start = i * one_piece_num
        sds.write(":WAVeform:STARt {}".format(start))
        sds.write("WAV:DATA?")
        recv_rtn = sds.read_raw().rstrip()
        block_start = recv_rtn.find(b'#')
        data_digit = int(recv_rtn[block_start + 1:block_start + 2])
        data_start = block_start + 2 + data_digit
        recv = list(recv_rtn[data_start:])
        recv_all += recv

    convert_data = []
    if adc_bit > 8:
        for i in range(0, int(len(recv_all) / 2)):
            data = recv_all[2 * i + 1] * 256 + recv_all[2 * i]
            convert_data.append(data)
    else:
        convert_data = recv_all
    volt_value = []
    for data in convert_data:
        if data > pow(2, adc_bit - 1) - 1:
            data = data - pow(2, adc_bit)
        else:
            pass
        volt_value.append(data)
    del recv, recv_all, convert_data
    gc.collect()

#converting raw data into values with units
    time_value = []
    for idx in range(0, len(volt_value)):
        volt_value[idx] = volt_value[idx] / vcode_per * float(vdiv) - float(ofst) #add *vertical_divisions/2 for slow osci !!!!!!!!! problem with original code: full measurement range is from -vcode_per to + vcode_per -> need to scale voltage per division with total number of divisions
        time_data = -(float(tdiv) * HORI_NUM / 2) + idx * interval + 1*float(trdl) #interval should be divided by 2? not anymore????????????? was maybe broken because of incorrect query for trigger status? trdl had wrong sign
        time_value.append(time_data)

    return time_value, volt_value, vcode_per, vdiv

####################################
def measurement_loop():
    global run_loop
    global valid_fit_count

    measurement_number = 0
    total_measurement_number = 0
    broken_measurement_number = 0
    all_measurement_data = []

    try:
        trigger_type = sds.query('TRIG:TYPE?')
        trigger_level = sds.query('TRIG:EDGE:LEV?')
        trigger_source = sds.query('TRIG:EDGE:SOUR?')
        print(f'Trigger type: {trigger_type} Trigger level: {trigger_level} Trigger Source: {trigger_source}')
        
    except:
        print(f'could not get trigger status')
    connected_channels.append(len(connected_channels))
    channel_use.append("empty")
    print(f"Channels: {connected_channels}{channel_use}")

    raw_data_file = None
    if store_raw_data:
        raw_data_file = h5py.File(data_folder+filename+file_suffix+".h5", "w")
        raw_data_file.attrs["description"] = "Raw waveform data"
        raw_data_file.attrs["time_unit"] = "ns"
        raw_data_file.attrs["voltage_unit"] = "V"
        raw_data_file.attrs["trigger_info"] = [trigger_type,trigger_level,trigger_source]
    
    while not stop_event.is_set():

        if data_taking_number is not None and measurement_number >= data_taking_number:
            print(f'Maximum number of measurements ({data_taking_number}) reached. Ending measurement!')
            stop_event.set()
            break
        if data_taking_time is not None and (time.time()-start_time) > data_taking_time:
            print(f'Maximum time of measuring ({data_taking_time} s) reached. Ending measurement!')
            stop_event.set()
            break

        if stop_event.is_set():
            run_loop=False
            print("Ending run loop!")
            break


        sds.write("TRIG:MODE SING")
        wait_time = 0
        measurement_number += 1
        total_measurement_number += 1
        print(f"\nThis is Measurement {measurement_number}.")
        measurement_data = {
            'measurement_number' : measurement_number,
            'channels' : {}  
        }

        while not stop_event.is_set():
            time.sleep(0.001)
            wait_time += 1
            trigger_stat = sds.query('TRIG:STAT?')
            if trigger_stat == "Stop\n":
                # print(f"Trigger detected! Waited {wait_time / 1000:.3f} second(s) for trigger.")
                
                break

        if take_screenshot:
                time.sleep(0.01)
                sds.write("SCDP") # https://telonic.co.uk/programming-example-sds-oscilloscope-save-a-copy-of-a-screen-image-via-python-pyvisa/
                screenshot = sds.read_raw()
                with open(screenshot_folder+f'{filename}{file_suffix}_{measurement_number:03d}.png','wb') as f:
                    f.write(bytearray(screenshot))


        # result = {}
        # result['measurement_number'] = measurement_number
        write_row = [str(measurement_number)]
        measurement_valid = True
        comparator_anomaly = False

        for channel in connected_channels:
            ch_use = channel_use[channel]
            
            # print(f"This is channel {channel}, {ch_use}")
            time.sleep(0.1)
            time_data, volt_data, vcode_per, vdiv = main_wf_data(channel)
            time_data = np.array(time_data) * 1e9 # convert to ns for better readability of results
            # print(f"Length of data: {len(time_data)}, {len(volt_data)}, mean volt={np.mean(np.array(volt_data))}")
            if len(time_data) == 0 or len(volt_data) == 0:
                print(f"\033[31mNo data for channel {channel}. Excluding measurement.\033[0m")
                measurement_valid = False
                measurement_number -=1
                broken_measurement_number += 1
                break

            if ch_use != "empty":
                
                if ch_use == 'amp':
                    time_data = np.asarray(time_data)
                    volt_data = np.asarray(volt_data)
                    time_data = time_data.reshape(-1,5,*time_data.shape[1:]).mean(axis=1)
                    volt_data = volt_data.reshape(-1,5,*volt_data.shape[1:]).mean(axis=1)
                
                measurement_data['channels'][f'ch{channel}'] = {
                    'type':ch_use,
                    'time_data':time_data,
                    'volt_data':np.array(volt_data),
                    'vcode_per':vcode_per,
                    'vdiv':vdiv
                }

                amplitude = np.max(volt_data) - np.min(volt_data)
                diff = np.diff(volt_data)
                i_on, i_off = np.argmax(diff), np.argmin(diff)
                t_on0, t_off0 = time_data[i_on], time_data[i_off]
                if t_off0 < t_on0:
                    t_on0, t_off0 = t_off0, t_on0

                if volt_data[-1] < volt_data[0]:
                    amplitude = -amplitude

                if ch_use == 'comp':
                    p0 = [amplitude, 30.0, 0.65, 2000.0, 0.01, 0.0]
                    fit_func = onoff_errorfunc
                    fit_bounds = (-np.inf, np.inf)
                elif ch_use == 'amp':
                    p0 = [2*amplitude, 10.0, 0.5, 100.0, 0.01, 0.5]
                    fit_func = onoff_errorfunc
                    fit_bounds = (-np.inf, np.inf)
                elif ch_use == "inj":
                    p0 = [amplitude,0,-0.1,1.4]
                    fit_func = errorfunc
                    fit_bounds = (-np.inf, np.inf)
                else:
                    print("Something went horribly wrong. Abort!!")
                    quit()

                try:
                    popt, pcov = scipy.optimize.curve_fit(
                        fit_func,
                        time_data,
                        volt_data,
                        p0=p0,
                        bounds=fit_bounds,
                        maxfev=1000
                    )
                except Exception as e:
                    print(f"\033[31mFit failed for channel {channel}: {e}\n Excluding measurement.\033[0m")
                    measurement_valid = False
                    measurement_number -=1
                    broken_measurement_number += 1
                    break

                popt_err = np.sqrt(np.diag(pcov))
                fit_line = fit_func(time_data, *popt)

                yerr = np.ones_like(volt_data) / vcode_per * float(vdiv) * vertical_divisions / 2
                reduced_chi_sq = (1.0 / (len(volt_data) - len(popt))) * np.sum(
                    ((volt_data - fit_line) / yerr) ** 2
                )

                measured_amplitude = np.max(volt_data) - np.min(volt_data)
                if ch_use == 'comp':
                    fitted_amp = abs(popt[0])
                    # if popt[3] > 3000:
                        # comparator_anomaly = True
                elif ch_use == 'amp':
                    fitted_amp = abs(popt[0])
                elif ch_use == 'inj':
                    fitted_amp = 2*abs(popt[0])
                else:
                    fitted_amp = np.nan
                

                if (fitted_amp < 0.3 * measured_amplitude or fitted_amp > 3.0 * measured_amplitude or popt_err[0] > 3*fitted_amp or reduced_chi_sq> chi_barrier):
                    print(f"\033[31mFit amplitude unrealistic for channel {channel} {ch_use}: fitted={fitted_amp:.3f}+/-{popt_err[0]:.3f}, measured={measured_amplitude:.3f}, chi2={reduced_chi_sq:.3f}.\nExcluding measurement.\033[0m")
                    measurement_valid = False
                    measurement_number -=1
                    broken_measurement_number += 1
                    break

                # time_start, time_end = int(0.19*len(time_data)), int(0.21*len(time_data))
                # result[f"type_ch{channel}"] = ch_use
                write_row.append(ch_use)  # these rows need to be in the correct order for the header
                # result[f'time_data_ch{channel}'] = time_data#[time_start:time_end]
                # result[f'volt_data_ch{channel}'] = volt_data#[time_start:time_end]
                # result[f'volt_error_ch{channel}'] = yerr#[time_start:time_end]
                # result[f'popt_ch{channel}'] = popt
                write_row.extend(str(val) for val in popt)
                # result[f'fit_err_ch{channel}'] = popt_err
                write_row.extend(str(val) for val in popt_err)
                # result[f'chi_red_ch{channel}'] = reduced_chi_sq
                write_row.append(str(reduced_chi_sq))
                measurement_data['channels'][f'ch{channel}']['volt_error'] = yerr
                measurement_data['channels'][f'ch{channel}']['popt'] = popt
                measurement_data['channels'][f'ch{channel}']['popt_err'] = popt_err
                measurement_data['channels'][f'ch{channel}']['chi2_red'] = reduced_chi_sq
            # else:
            #     # print("empty channel")
            #     print("")

        if measurement_valid:   
            with open(fit_file_path, 'a') as f:
                f.write('\t'.join(write_row)+'\n')

            valid_fit_count += 1

            if store_raw_data and raw_data_file is not None:
                meas_group_name = f"measurement_{measurement_number:05d}"
                meas_group = raw_data_file.create_group(meas_group_name)
                meas_group.attrs["measurement_number"] = measurement_number
                meas_group.attrs["total_measurement_number"] = total_measurement_number
                for ch_name, ch_data in measurement_data["channels"].items():

                    ch_group = meas_group.create_group(ch_name)

                    ch_group.attrs["type"] = ch_data["type"]
                    ch_group.attrs["vcode_per"] = ch_data["vcode_per"]
                    ch_group.attrs["vdiv"] = ch_data["vdiv"]

                    ch_group.create_dataset(
                        "time_data",
                        data=ch_data["time_data"],
                        compression="gzip",
                        compression_opts=4
                    )

                    ch_group.create_dataset(
                        "volt_data",
                        data=ch_data["volt_data"],
                        compression="gzip",
                        compression_opts=2 # compression level from 0 (none) to 9(max)
                    )

                    ch_group.create_dataset(
                        "volt_err",
                        data=ch_data['volt_error'],
                        compression="gzip",
                        compression_opts=4
                    )

                raw_data_file.flush()

            if comparator_anomaly: # Thanks to Google gemini
                comp_ch = None
                amp_ch = None
                for ch in connected_channels[:-1]:
                    # if result.get(f"type_ch{ch}") == "comp": comp_ch = ch
                    # if result.get(f"type_ch{ch}") == "amp": amp_ch = ch
                    if measurement_data["channels"][f'ch{ch}']['type'] == "comp": comp_ch = ch
                    elif measurement_data["channels"][f'ch{ch}']['type'] == "amp": amp_ch = ch
                
                # 2. If both exist, check the turn-off time and plot
                if comp_ch is not None and amp_ch is not None:
                    # popt_comp = result[f'popt_ch{comp_ch}']
                    comp_channel_data = measurement_data["channels"][f"ch{comp_ch}"] 
                    amp_channel_data = measurement_data["channels"][f"ch{amp_ch}"]
                    popt_comp = comp_channel_data['popt']
                    popt_amp = amp_channel_data['popt']
                    
                    if popt_comp[3] > 3000:
                        print(f"\033[33mFun Fact: The comparator turn-off time is > 3000 ns: {popt_comp[3]:.2f} ns.\nI will make you a nice plot :)\033[0m")
                        
                        fig, ax1 = plt.subplots(figsize=(15, 6))
                        ax2 = ax1.twinx()  # Create a secondary y-axis that shares the same x-axis
                        
                        # --- Comparator Data (Left Y-Axis) ---
                        # t_comp = result[f'time_data_ch{comp_ch}']
                        # v_comp = result[f'volt_data_ch{comp_ch}']
                        # err_comp = result[f'volt_error_ch{comp_ch}']
                        # chi2_comp = result[f'chi_red_ch{comp_ch}']
                        # popt_err_comp = result[f'fit_err_ch{comp_ch}']
                        t_comp = comp_channel_data['time_data']
                        v_comp = comp_channel_data['volt_data']
                        err_comp = comp_channel_data['volt_error']
                        chi2_comp = comp_channel_data['chi2_red']
                        popt_err_comp = comp_channel_data['popt_err']
                        
                        ax1.errorbar(t_comp, v_comp, yerr=err_comp, fmt='.', ms=0.5, color='black', label='Comparator Data')
                        ax1.plot(t_comp, onoff_errorfunc(t_comp, *popt_comp), color='red', ls='-', label=f'Comp Fit\n' + r'$\chi^2_{{red}}$ = ' + f'{chi2_comp:.3f}')
                        ax1.axvline(popt_comp[1], color='blue', ls='--', linewidth=2, label=r'$t_{on}$ = ' + f'{popt_comp[1]:.3e} +/- {popt_err_comp[1]:.3e} ns')
                        ax1.axvline(popt_comp[3], color='blue', ls='--', linewidth=2, label=r'$t_{off}$ = ' + f'{popt_comp[3]:.3e} +/- {popt_err_comp[3]:.3e} ns')
                        
                        ax1.set_xlabel('Time (ns)')
                        ax1.set_ylabel('Comparator Voltage (V)', color='red')
                        ax1.tick_params(axis='y', labelcolor='red')
                        
                        # --- Amplifier Data (Right Y-Axis) ---
                        # t_amp = result[f'time_data_ch{amp_ch}']
                        # v_amp = result[f'volt_data_ch{amp_ch}']
                        # err_amp = result[f'volt_error_ch{amp_ch}']
                        # popt_amp = result[f'popt_ch{amp_ch}']
                        # chi2_amp = result[f'chi_red_ch{amp_ch}']
                        t_amp = amp_channel_data['time_data']
                        v_amp = amp_channel_data['volt_data']
                        err_amp = amp_channel_data['volt_error']
                        chi2_amp = amp_channel_data['chi2_red']
                        popt_err_amp = amp_channel_data['popt_err']
                        
                        ax2.errorbar(t_amp, v_amp, yerr=err_amp, fmt='.', ms=0.5, color='gray', alpha=0.5, label='Amplifier Data')
                        ax2.plot(t_amp, onoff_errorfunc(t_amp, *popt_amp), color='green', ls='-', label=f'Amp Fit\n' + r'$\chi^2_{{red}}$ = ' + f'{chi2_amp:.3f}')
                        # ax2.plot(t_amp, csa_large_signal(t_amp, *popt_amp), color='green', ls='-', label=f'Amp Fit\n' + r'$\chi^2_{{red}}$ = ' + f'{chi2_amp:.3f}')
                        
                        ax2.set_ylabel('Amplifier Voltage (V)', color='green')
                        ax2.tick_params(axis='y', labelcolor='green')
                        
                        # --- Combine Legends and Save ---
                        lines_1, labels_1 = ax1.get_legend_handles_labels()
                        lines_2, labels_2 = ax2.get_legend_handles_labels()
                        # Merging handles/labels so they all appear in one unified legend box
                        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')
                        
                        ax1.grid()
                        fig.tight_layout()
                        plt.savefig(plots_folder + filename + file_suffix + f'_{measurement_number:03d}' + '.png', dpi=200)
                        plt.close(fig)

            print("\033[32mMeasurement done.\033[0m")
        else:
            print("\033[33mMeasurement excluded\033[0m")
        
    if raw_data_file is not None:
        raw_data_file.close()
        print(f"Saved raw data as {data_folder+filename+file_suffix+".h5"}")
        # if store_raw_data: all_measurement_data.append(measurement_data)

    # if all_measurement_data and store_raw_data and len(stored_results)!= 0: 
    #     np.savez_compressed(data_folder+filename+file_suffix+'.npz', measurements=np.array(all_measurement_data, dtype=object))
    #     print(f"Saved data as '{filename}.npz'")

    if valid_fit_count == 0 and os.path.exists(fit_file_path):
        print("No valid fits, delete data files")
        os.remove(fit_file_path)
        try:
            os.remove(data_folder+filename+file_suffix+".h5")
        except Exception as e:
            print(f"Could not delete raw data file {data_folder+filename+file_suffix+".h5"}\n{e}")
    ###    
    print(f"Summary of Measurement: Total Tries: {total_measurement_number}\n\tWorking measurements: {measurement_number}, Broken measurements: {broken_measurement_number} = {(broken_measurement_number/total_measurement_number)*100:.1f} % of measurements")

def make_wave_plots(stored_results, chip=""):
    print("Plotting wave forms ...")
    if chip != "":
        chip = " | " + chip
    max_plots_per_figure = 6
    n_tot = len(stored_results)

    for fig_i, j in enumerate(range(0, n_tot, max_plots_per_figure), start=1):
        batch = stored_results[j:j+max_plots_per_figure]
        n = len(batch)

        ncols = 2
        nrows = int(np.ceil(n / ncols))

        fig, axes = plt.subplots(nrows, ncols, figsize=(20, 6 * nrows), sharex=True, sharey=False)
        axes = np.atleast_1d(axes).flatten()
        

        for ax, measurement in zip(axes, batch):
            ax_comp = ax
            ax_inj = ax
            ax_amp = ax.twinx()

            ax_amp.set_zorder(1)
            ax_comp.set_zorder(2)

            ax_comp.patch.set_alpha(0)

            comp_handles, comp_labels = [], []
            amp_handles, amp_labels = [], []
            inj_handles, inj_labels = [],  []

            meas_channels = measurement["channels"]
            measurement_number = measurement["measurement_number"]

            for i, channel in enumerate(connected_channels):
                ch_key = f"ch{channel}"

                if ch_key not in meas_channels:
                    continue

                ch_data = meas_channels[ch_key]
                ch_use = ch_data["type"]

                time_data = np.asarray(ch_data["time_data"])
                volt_data = np.asarray(ch_data["volt_data"])
                volt_error = np.asarray(ch_data["volt_error"])
                popt = np.asarray(ch_data["popt"])
                popt_err = np.asarray(ch_data["popt_err"])
                chi_red = ch_data["chi2_red"]

                plot_time = np.linspace(np.min(time_data), np.max(time_data), 1000)

                if ch_use == "comp":
                    use_ax = ax_comp

                    data_plot_comp = use_ax.errorbar(
                        time_data, volt_data, yerr=volt_error,
                        fmt='.', ms=0.5, color=ch_colors[channel], zorder=1
                    )

                    fit_plot_comp, = use_ax.plot(
                        plot_time, onoff_errorfunc(plot_time, *popt),
                        color=ch_fit_colors[channel], ls='-', zorder=2
                    )

                    if chi_red < chi_barrier:
                        use_ax.fill_between(
                            plot_time,
                            onoff_errorfunc(plot_time, *(popt - popt_err)),
                            onoff_errorfunc(plot_time, *(popt + popt_err)),
                            color=ch_fit_colors[channel], alpha=0.5, zorder=2
                        )

                    v1_comp = use_ax.axvline(popt[1], color=ch_colors[channel], ls='--', linewidth=0.5, zorder=0)
                    v2_comp = use_ax.axvline(popt[3], color=ch_colors[channel], ls='--', linewidth=0.5, zorder=0)

                    comp_handles += [data_plot_comp, fit_plot_comp, v1_comp, v2_comp]
                    comp_labels += [
                        'Comparator data',
                        f'Fit (comp)\n' + r'$\chi^2_{{red}}$ = ' + f'{chi_red:.3f}',
                        r'$t_{on}$ = ' + f'{popt[1]:.3e} +/- {popt_err[1]:.3e} ns',
                        r'$t_{off}$ = ' + f'{popt[3]:.3e} +/- {popt_err[3]:.3e} ns'
                    ]

                elif ch_use == "amp":
                    use_ax = ax_amp

                    data_plot_amp = use_ax.errorbar(
                        time_data, volt_data, yerr=volt_error,
                        fmt='.', ms=0.5, color=ch_colors[channel],
                        zorder=1, label='Amplifier data'
                    )

                    fit_plot_amp, = use_ax.plot(
                        plot_time, onoff_errorfunc(plot_time, *popt),
                        color=ch_fit_colors[channel], ls='-', zorder=2
                    )

                    fit_data_min = onoff_errorfunc(plot_time, *(popt - popt_err))
                    fit_data_max = onoff_errorfunc(plot_time, *(popt + popt_err))

                    if chi_red < chi_barrier and max(fit_data_min) < 10 and max(fit_data_max) < 10:
                        use_ax.fill_between(
                            plot_time, fit_data_min, fit_data_max,
                            color=ch_fit_colors[channel], alpha=0.3, zorder=2
                        )

                    amp_handles += [data_plot_amp, fit_plot_amp]
                    amp_labels += [
                        'Amplifier data',
                        f'Fit (amp)\n' + r'$\chi^2_{{red}}$ = ' + f'{chi_red:.3f}'
                    ]

                elif ch_use == "inj":
                    use_ax = ax_inj

                    data_plot_inj = use_ax.errorbar(
                        time_data, volt_data, yerr=volt_error,
                        fmt='.', ms=0.5, color=ch_colors[channel],
                        zorder=1, label='Injection Data'
                    )

                    fit_plot_inj, = use_ax.plot(
                        plot_time, errorfunc(plot_time, *popt),
                        color=ch_fit_colors[channel], ls='-', zorder=1
                    )

                    fit_data_min = errorfunc(plot_time, *(popt - popt_err))
                    fit_data_max = errorfunc(plot_time, *(popt + popt_err))

                    if chi_red < chi_barrier and max(fit_data_min) < 10 and max(fit_data_max) < 10:
                        use_ax.fill_between(
                            plot_time, fit_data_min, fit_data_max,
                            color=ch_fit_colors[channel], alpha=0.3, zorder=2
                        )

                    v1_inj = use_ax.axvline(popt[1], color=ch_colors[channel], ls='--', linewidth=0.5, zorder=0)

                    inj_handles += [data_plot_inj, fit_plot_inj, v1_inj]
                    inj_labels += [
                        'Injection Data',
                        f'Fit (inj)\n' + r'$\chi^2_{{red}}$ = ' + f'{chi_red:.3f}',
                        r'$t_{on}$ = ' + f'{popt[1]:.3e} +/- {popt_err[1]:.3e} ns'
                    ]

            ax_comp.set_xlabel("Time (ns)")
            ax_comp.set_ylabel("Voltage (V)", color="black")
            ax_amp.set_ylabel("Voltage (V)", color="red")

            ax_comp.tick_params(axis='y', labelcolor="black")
            ax_amp.tick_params(axis='y', labelcolor="red")

            ax_comp.xaxis.set_tick_params(labelbottom=True)
            ax_comp.set_title(f"Measurement {measurement_number}")
            ax_comp.grid(True, color='grey', zorder=0)
            ax_amp.grid(True, color='red', zorder=0)

            handles = comp_handles + amp_handles + inj_handles
            labels = comp_labels + amp_labels + inj_labels
            ax_comp.legend(handles, labels, loc='upper right', labelspacing=1.0, handlelength=1, fontsize=8)

        for ax in axes[len(batch):]:
            ax.set_visible(False)

        fig.suptitle(
            f"Measurements {batch[0]['measurement_number']} to {batch[-1]['measurement_number']} from {len(stored_results)}\n"
            f"({filename[12:]}{chip})",
            fontsize=16
        )

        fig.tight_layout(h_pad=3.0, w_pad=5.0, rect=[0, 0, 1, 0.98])
        plt.savefig(plots_folder + f"{filename}{file_suffix}_{fig_i:03d}.pdf")
        print(f"Saved Plot as '{filename}{file_suffix}_{fig_i:03d}.pdf'")
        plt.close(fig)

def make_hist_plots(stored_results, plot_correlation=False, plot_pulse_width = False, chip=""):
    print("Plotting histograms ...")
    if chip!= "":
            chip = " | "+ chip
    if elapsed_time != 0.0:
        elapsed_time_string = f' | {elapsed_time:.2f} s'
    else:
        elapsed_time_string = ' | offline'

    measurement_numbers = []

    comp_turn_on_times = {}
    comp_turn_off_times = {}
    comp_pulse_widths = {}
    comp_rise_times = {}
    comp_fall_times = {}
    comp_amplitude = {}

    amp_turn_on_times = {}
    amp_turn_off_times = {}
    amp_pulse_widths = {}
    amp_rise_times = {}
    amp_fall_times = {}
    amp_amplitude = {}

    inj_t0 = {}
    comp_t0_to_inj = {}
    amp_t0_to_inj = {}

    for channel in connected_channels:
        ch_use = channel_use[channel]

        if ch_use == "comp":
            comp_turn_on_times[channel] = []
            comp_turn_off_times[channel] = []
            comp_pulse_widths[channel] = []
            comp_rise_times[channel] = []
            comp_fall_times[channel] = []
            comp_amplitude[channel] = []

        elif ch_use == "amp":
            amp_turn_on_times[channel] = []
            amp_turn_off_times[channel] = []
            amp_pulse_widths[channel] = []
            amp_rise_times[channel] = []
            amp_fall_times[channel] = []
            amp_amplitude[channel] = []

        elif ch_use == "inj":
            inj_t0[channel] = []

    inj_channels = [ch for ch in connected_channels if channel_use[ch] == "inj"]
    comp_channels = [ch for ch in connected_channels if channel_use[ch] == "comp"]
    amp_channels = [ch for ch in connected_channels if channel_use[ch] == "amp"]

    for ch in comp_channels:
        comp_t0_to_inj[ch] = []

    for ch in amp_channels:
        amp_t0_to_inj[ch] = []

    for measurement in stored_results:
        all_inj_t0 = {}
        measurement_numbers.append(measurement["measurement_number"])

        meas_channels = measurement["channels"]

        for channel in connected_channels:
            ch_key = f"ch{channel}"

            if ch_key not in meas_channels:
                continue

            ch_data = meas_channels[ch_key]
            ch_use = ch_data["type"]
            popt = ch_data["popt"]

            if ch_use == "comp":
                A, t_on, a_on, t_off, a_off, C = popt

                if t_off < t_on:
                    t_on, t_off = t_off, t_on
                    a_on, a_off = a_off, a_on

                comp_turn_on_times[channel].append(t_on)
                comp_turn_off_times[channel].append(t_off)
                comp_pulse_widths[channel].append(t_off - t_on)
                comp_rise_times[channel].append(1.812 / abs(a_on))
                comp_fall_times[channel].append(1.812 / abs(a_off))
                comp_amplitude[channel].append(A)

            elif ch_use == "amp":
                A, t_on, a_on, t_off, a_off, C = popt

                if t_off < t_on:
                    t_on, t_off = t_off, t_on
                    a_on, a_off = a_off, a_on

                amp_turn_on_times[channel].append(t_on)
                amp_turn_off_times[channel].append(t_off)
                amp_pulse_widths[channel].append(t_off - t_on)
                amp_rise_times[channel].append(1.812 / abs(a_on))
                amp_fall_times[channel].append(1.812 / abs(a_off))
                amp_amplitude[channel].append(A)

            elif ch_use == "inj":
                t0_inj = popt[3]
                inj_t0[channel].append(t0_inj)
                all_inj_t0[channel] = t0_inj

        if len(inj_channels) > 0:
            ref_inj_ch = inj_channels[0]

            if ref_inj_ch in all_inj_t0:
                t0_ref = all_inj_t0[ref_inj_ch]

                for ch in comp_channels:
                    ch_key = f"ch{ch}"
                    if ch_key in meas_channels:
                        popt = meas_channels[ch_key]["popt"]
                        comp_t0_to_inj[ch].append(popt[1] - t0_ref)

                for ch in amp_channels:
                    ch_key = f"ch{ch}"
                    if ch_key in meas_channels:
                        popt = meas_channels[ch_key]["popt"]
                        amp_t0_to_inj[ch].append(popt[1] - t0_ref)

 #  -------- Comparator Figure -------------------- ####  
    if len(comp_channels) > 0:
        n_comp = len(comp_channels)
        n_cols = 3
        n_rows = 2 * n_comp

        fig_comp, axes_comp = plt.subplots(
            n_rows,
            n_cols,
            figsize=(4.6 * n_cols, 3.8 * n_rows),
            constrained_layout=True,
            squeeze=False
        )


        for row, ch in enumerate(comp_channels):
            base_row = 2 * row

            plots = [
                (
                    comp_rise_times[ch],
                    f"Comparator Ch{ch}: Rise time",
                    "Time", "ns",
                    "tab:green",
                ),
                (
                    comp_fall_times[ch],
                    f"Comparator Ch{ch}: Fall time",
                    "Time", "ns",
                    "tab:red",
                ),
                (
                    comp_turn_on_times[ch],
                    f"Comparator Ch{ch}: Turn-on moment",
                    "Time", "ns",
                    "tab:blue",
                ),
                # (
                #     comp_turn_off_times[ch],
                #     f"Comparator Ch{ch}: Turn-off moment",
                #     "Time", "ns",
                #     "tab:orange",
                # ),
                (
                    comp_pulse_widths[ch],
                    f"Comparator Ch{ch}: Pulse width",
                    "Time", "ns",
                    "tab:orange",
                ),
                (
                    comp_amplitude[ch],
                    f"Comparator Ch{ch}: Amplitude",
                    "Voltage", "V",
                    "tab:purple",
                ),
                
            ]

            if len(inj_channels) > 0:
                plots.append(
                    (
                        comp_t0_to_inj[ch],
                        f"Comp Ch {ch}: $t_{{on}} - t_{{0,inj}}$",
                        "Time difference", "ns",
                        "tab:cyan",
                    )
                )

            # Fill the 2 x 3 block for this channel
            for i, plot_args in enumerate(plots):
                ax = axes_comp[base_row + i // 3, i % 3]
                make_hist(ax, *plot_args)

            # Remove unused axes in this channel block
            for i in range(len(plots), 6):
                ax = axes_comp[base_row + i // 3, i % 3]
                fig_comp.delaxes(ax)
   
        fig_comp.suptitle(
            f"Comparator Timing Distributions\n"
            f"No. of meas.: {len(stored_results)}" + elapsed_time_string + f" | {filename[12:]}"+chip,
            fontsize=16
        )

        plt.savefig(plots_folder + filename + file_suffix + "_hist_comp.pdf")
        print(f"Saved Plot as '{filename}{file_suffix}_hist_comp.pdf'\n")
        plt.close(fig_comp)

 #  -------- Amplifier Figure -------------------- ####
    if len(amp_channels) > 0:
        n_amp = len(amp_channels)
        n_cols = 3
        n_rows = 2 * n_amp

        with open(data_folder + filename + file_suffix +"_amp_hist_fit.txt", "w") as f:
            f.write("Channel\tPlot\tUnit\tA_fit\tA_err\tmu_fit\tmu_err\tsigma_fit\tsigma_err\n")

        fig_amp, axes_amp = plt.subplots(
            n_rows,
            n_cols,
            figsize=(4.6 * n_cols, 3.8 * n_rows),
            constrained_layout=True,
            squeeze=False
        )


        for row, ch in enumerate(amp_channels):
            base_row = 2 * row

            plots = [
                (
                    amp_rise_times[ch],
                    f"Amplifier Ch{ch}: Rise time",
                    "Time", "ns",
                    "tab:green",
                ),
                (
                    amp_fall_times[ch],
                    f"Amplifier Ch{ch}: Fall time",
                    "Time", "ns",
                    "tab:red",
                ),
                (
                    amp_amplitude[ch],
                    f"Amplifier Ch{ch}: Amplitude",
                    "Voltage", "V",
                    "tab:blue",
                ),
                (
                    amp_pulse_widths[ch],
                    f"Amplifier Ch{ch}: Pulse width",
                    "Time", "ns",
                    "tab:purple",
                ),
                (
                    amp_turn_off_times[ch],
                    f"Amplifier Ch{ch}: Turn-off moment",
                    "Time", "ns",
                    "tab:orange",
                ),
                
            ]

            if len(inj_channels) > 0:
                plots.append(
                    (
                        amp_t0_to_inj[ch],
                        f"Amplifier Ch {ch}: $t_{{on}} - t_{{0,inj}}$",
                        "Time difference", "ns",
                        "tab:cyan",
                    )
                )
            else:
                plots.append(
                (
                    amp_turn_on_times[ch],
                    f"Amplifier Ch{ch}: Turn-on moment",
                    "Time", "ns",
                    "tab:cyan",
                )
                )

            # Fill the 2 x 3 block for this channel
            for i, plot_args in enumerate(plots):
                ax = axes_amp[base_row + i // 3, i % 3]
                A_fit, A_err, mu_fit, mu_err, sigma_fit, sigma_err = make_hist(ax, *plot_args)
                with open(data_folder+filename+file_suffix+"_amp_hist_fit.txt", "a") as f:
                    f.write(f"{ch}\t{i}\t{plot_args[3]}\t{A_fit}\t{A_err}\t{mu_fit}\t{mu_err}\t{sigma_fit}\t{sigma_err}\n")

            # Remove unused axes in this channel block
            for i in range(len(plots), 6):
                ax = axes_amp[base_row + i // 3, i % 3]
                fig_comp.delaxes(ax)
        
        fig_amp.suptitle(
            f"Amplifier Timing Distributions\n"
            f"No. of meas.: {len(stored_results)}" + elapsed_time_string + f" | {filename[12:]}"+chip,
            fontsize=16
        )

        plt.savefig(plots_folder + filename + file_suffix + "_hist_amp.pdf")
        print(f"Saved Plot as '{filename}{file_suffix}_hist_amp.pdf'\n")
        plt.close(fig_amp)

    # ------ Heatmap -------- ####
    if plot_correlation:
        ch_c = comp_channels[0]
        ch_a = amp_channels[0]
        # 2. Define the labels and build the 2D arrays (Rows = Variables, Columns = Observations)
        # Note: We only include t0_to_inj if injection channels exist
        comp_labels = ['Turn On', 'Turn Off', 'Pulse Width', 'Rise Time', 'Fall Time', 'Amplitude']
        comp_data_list = [
            comp_turn_on_times[ch_c], comp_turn_off_times[ch_c], comp_pulse_widths[ch_c],
            comp_rise_times[ch_c], comp_fall_times[ch_c], comp_amplitude[ch_c]
        ]
        
        amp_labels = ['Turn On', 'Turn Off', 'Pulse Width', 'Rise Time', 'Fall Time', 'Amplitude']
        amp_data_list = [
            amp_turn_on_times[ch_a], amp_turn_off_times[ch_a], amp_pulse_widths[ch_a],
            amp_rise_times[ch_a], amp_fall_times[ch_a], amp_amplitude[ch_a]
        ]

        # Append injection time differences if they were calculated
        if len(inj_channels) > 0:
            comp_labels.append('t0 to Inj')
            comp_data_list.append(comp_t0_to_inj[ch_c])
            
            amp_labels.append('t0 to Inj')
            amp_data_list.append(amp_t0_to_inj[ch_a])

        # Convert to NumPy arrays
        comp_data = np.array(comp_data_list)
        amp_data = np.array(amp_data_list)
        make_correlation_map(comp_data,comp_labels,amp_data,amp_labels, chip=chip)

    if plot_pulse_width and len(comp_channels)>0:
        print("Making Plot Comparator Width vs. Measurement number ...")
        for row, ch in enumerate(comp_channels):
            try:
                fig, ax = plt.subplots(1, len(comp_channels), figsize=(12,6*len(comp_channels)))

                ax.scatter(measurement_numbers, comp_pulse_widths[ch], s=1, color='black', label='Data from fits')
                ax.set_xlabel('Measuement Number')
                ax.set_ylabel('Comparator Pulse Width (ns)')
                ax.legend(title=f"Ch{ch}")
                ax.grid()
                fig.suptitle(
                    f"Comparator Pulse width per Measurement\n"
                    f"No. of meas.: {len(stored_results)}" + f" | {filename[12:]}"+chip,
                    fontsize=12
                )

                plt.savefig(plots_folder + filename + file_suffix + "_pulse_vs_n.pdf")
                print(f"Saved plot as {plots_folder + filename + file_suffix + "_pulse_vs_n.pdf"}\n")
            except Exception as e:
                print(f"Could not make plot Width vs. Measurement number\n{e}")

def make_correlation_map(comp_data, comp_labels, amp_data, amp_labels, chip=""):
    print("Making Correlation Heatmap...")
    corr_comp = np.corrcoef(comp_data)
    corr_amp = np.corrcoef(amp_data)

    combined_data = np.vstack((comp_data, amp_data))
    corr_combined = np.corrcoef(combined_data)

    n_comp_vars = len(comp_labels)
    corr_cross = corr_combined[0:n_comp_vars, n_comp_vars:] 

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    heatmap_kwargs = {
        'annot': True, 
        'fmt': ".2f", 
        'cmap': "coolwarm", 
        'vmin': -1, 
        'vmax': 1, 
        'square': True,
        'cbar_kws': {"shrink": .8}
    }

    sns.heatmap(corr_comp, ax=axes[0], xticklabels=comp_labels, yticklabels=comp_labels, **heatmap_kwargs)
    axes[0].set_title(f"Comparator Correlation")

    sns.heatmap(corr_amp, ax=axes[1], xticklabels=amp_labels, yticklabels=amp_labels, **heatmap_kwargs)
    axes[1].set_title(f"Amplifier Correlation")

    sns.heatmap(corr_cross, ax=axes[2], xticklabels=amp_labels, yticklabels=comp_labels, **heatmap_kwargs)
    axes[2].set_title(f"Cross-Correlation")
    axes[2].set_ylabel("Comparator Variables")
    axes[2].set_xlabel("Amplifier Variables")

    fig.suptitle(
            f"Fit parameter Correlations\n"
            f"No. of meas.: {len(stored_results)}" + f" | {filename[12:]}"+chip,
            fontsize=12
        )
    
    fig.subplots_adjust(left=0.02, right=0.98, top=0.85, bottom=0.2, wspace=0.1)
    plt.savefig(plots_folder + filename + file_suffix + "_corr.pdf")
    print(f"Saved Plot as '{filename}{file_suffix}_corr.pdf'\n")
    plt.close(fig)

# ---------- Helper for histogram plotting ----------  AI generated with KIWIS ChatGPT-5.5 Full
def make_hist(
    ax,
    data,
    title,
    xlabel="",
    unit="",
    color="blue",
    max_bins=50,
    min_bins=12,
    width_factor=0.75,
):
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]

    if len(data) == 0:
        ax.set_title(title, fontsize=13)
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        ax.grid(alpha=0.3)
        return

    if len(data) == 1 or np.all(data == data[0]):
        counts, bin_edges, _ = ax.hist(
            data,
            bins=1,
            color=color,
            alpha=0.85,
            edgecolor="black",
            linewidth=0.8,
        )

        ax.set_title(title, fontsize=13)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Count")
        ax.set_ylim(0, np.max(counts) + 1)
        ax.grid(alpha=0.3)
        return

    # Use full stored data range
    hist_min = np.min(data)
    hist_max = np.max(data)

    # Freedman-Diaconis bin width
    q25, q75 = np.percentile(data, [25, 75])
    iqr = q75 - q25

    if iqr > 0:
        bin_width = width_factor * 2 * iqr / (len(data) ** (1 / 3))
    else:
        data_range = hist_max - hist_min
        bin_width = data_range / min_bins if data_range > 0 else 1.0

    data_range = hist_max - hist_min

    if bin_width <= 0 or not np.isfinite(bin_width):
        n_bins = min_bins
    else:
        n_bins = int(np.ceil(data_range / bin_width))

    n_bins = max(min_bins, n_bins)
    n_bins = min(max_bins, n_bins)

    counts, bin_edges, _ = ax.hist(
        data,
        bins=n_bins,
        range=(hist_min, hist_max),
        color=color,
        alpha=0.85,
        edgecolor="black",
        linewidth=0.8,
        label="Data",
    )

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    count_err = np.sqrt(counts)
    ax.errorbar(bin_centers,counts,yerr=count_err,fmt='none',ecolor='black',capsize=2,linewidth=1,alpha=0.85, label=r"$\sqrt{N}$")

    fit_err = count_err.copy()
    fit_err[fit_err==0] = 1.0 # to avoid zero for fit

    A_fit, A_err, mu_fit, mu_err, sigma_fit, sigma_err = 0,0,0,0,0,0

    try:
        p0 = [np.max(counts), np.mean(data), np.std(data)]
        popt, pcov = scipy.optimize.curve_fit(
            gaussian,
            bin_centers,
            counts,
            p0=p0,
            sigma=fit_err,
            absolute_sigma=True,
            maxfev=1000
        )

        A_fit, mu_fit, sigma_fit = popt
        perr = np.sqrt(np.diag(pcov))
        A_err, mu_err, sigma_err = perr

        x_fit = np.linspace(hist_min, hist_max, 500)
        y_fit = gaussian(x_fit, *popt)

        model_counts = gaussian(bin_centers, *popt)
        chi2_ndf = np.sum(((counts - model_counts) / fit_err) ** 2)/(len(counts) - len(popt))

        ax.plot(
            x_fit,
            y_fit,
            color="black",
            linewidth=2,
            alpha=0.5,
            label=(
                "Gauss fit:\n"
                + rf"$\mu=({mu_fit:.3f}\pm{mu_err:.3f})$ {unit})"
                + "\n"
                + rf"$\sigma=({sigma_fit:.3f}\pm{sigma_err:.3f})$ {unit}"
                + "\n"
                + rf"$\chi^2/\mathrm{{ndf}}={chi2_ndf:.2f}$"
            ),
        )

    except Exception as e:
        ax.text(
            0.98, 0.65,
            "Gaussian fit failed!",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8,
            color="red",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
        )
        print(f"Gaussian fit to Histogram failed: {e}")

    max_labels = 10
    if len(bin_centers) <= max_labels:
        tick_positions = bin_centers
        tick_labels = [f"{x:.3f}" for x in bin_centers]
    else:
        step = int(np.ceil(len(bin_centers) / max_labels))
        tick_positions = bin_centers[::step]
        tick_labels = [f"{x:.3f}" for x in tick_positions]

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=45, fontsize=7)

    ax.set_ylim(0, max(np.max(counts) + 1, ax.get_ylim()[1]))

    ax.set_title(title, fontsize=13)
    ax.set_xlabel(xlabel+f" ({unit})")
    ax.set_ylabel("Count")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, loc="upper right")

    return A_fit, A_err, mu_fit, mu_err, sigma_fit, sigma_err

def load_measurements_from_h5(h5_path, txt_path, connected_channels, fit_param_names):

    stored_results = []
    
    # Load fit results from txt file
    fit_results_dict = {}
    try:
        with open(txt_path, "r") as f:
            lines = f.readlines()
        
        for line in lines[1:]:  # skip header
            line = line.strip()
            if not line:
                continue
            
            parts = line.split("\t")
            idx = 0
            
            meas_num = int(float(parts[idx]))
            fit_results_dict[meas_num] = {}
            idx += 1
            
            for channel in connected_channels:
                ch_use = parts[idx]
                idx += 1
                
                n_par = len(fit_param_names[ch_use])
                
                popt = np.array([float(x) for x in parts[idx:idx + n_par]])
                idx += n_par
                
                popt_err = np.array([float(x) for x in parts[idx:idx + n_par]])
                idx += n_par
                
                chi_red = float(parts[idx])
                idx += 1
                
                fit_results_dict[meas_num][f"ch{channel}"] = {
                    "popt": popt,
                    "popt_err": popt_err,
                    "chi2_red": chi_red
                }
    except Exception as e:
        print(f"Warning: Could not load fit results from txt file: {e}")
    
    # Load raw data from h5 file
    try:
        raw_data_file = h5py.File(h5_path, "r")
        measurement_ids = sorted(raw_data_file.keys())
        
        for meas_id in measurement_ids:
            meas_group = raw_data_file[meas_id]
            meas_num = int(meas_group.attrs["measurement_number"])
            
            measurement_data = {
                "measurement_number": meas_num,
                "channels": {}
            }
            
            channel_ids = sorted(meas_group.keys())
            
            for ch_id in channel_ids:
                ch_group = meas_group[ch_id]
                ch_type = ch_group.attrs["type"]
                
                if isinstance(ch_type, bytes):
                    ch_type = ch_type.decode()
                
                channel_data = {
                    "type": ch_type,
                    "time_data": ch_group["time_data"][:],
                    "volt_data": ch_group["volt_data"][:],
                    "volt_error": ch_group["volt_err"][:],
                    "vcode_per": ch_group.attrs["vcode_per"],
                    "vdiv": ch_group.attrs["vdiv"]
                }
                
                # Merge with fit results if available
                if meas_num in fit_results_dict and ch_id in fit_results_dict[meas_num]:
                    fit_data = fit_results_dict[meas_num][ch_id]
                    channel_data["popt"] = fit_data["popt"]
                    channel_data["popt_err"] = fit_data["popt_err"]
                    channel_data["chi2_red"] = fit_data["chi2_red"]
                
                measurement_data["channels"][ch_id] = channel_data
            
            stored_results.append(measurement_data)
        
        raw_data_file.close()
    except Exception as e:
        print(f"Error loading h5 file: {e}")
        return []
    
    return stored_results

def load_fit_results_from_txt(fit_file_path, connected_channels, fit_param_names):
    stored_results = []

    with open(fit_file_path, "r") as f:
        lines = f.readlines()

    # skip header
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue

        parts = line.split("\t")
        idx = 0

        measurement_number = int(float(parts[idx]))
        measurement_data = {
            "measurement_number": measurement_number,
            "channels": {}
        }
        idx += 1

        for channel in connected_channels:
            ch_use = parts[idx]
            idx += 1

            n_par = len(fit_param_names[ch_use])

            popt = np.array([float(x) for x in parts[idx:idx + n_par]])
            idx += n_par

            popt_err = np.array([float(x) for x in parts[idx:idx + n_par]])
            idx += n_par

            chi_red = float(parts[idx])
            idx += 1

            measurement_data["channels"][f"ch{channel}"] = {
                "type": ch_use,
                "popt": popt,
                "popt_err": popt_err,
                "chi2_red": chi_red
            }

        stored_results.append(measurement_data)

    return stored_results

def offline_fit(raw_data_path=""):
    raw_data_file = None
    loaded_data = []
    stored_results_offline = []
    valid_fit_count_offline = 0
    broken_measurement_number = 0

    if raw_data_path.endswith(".npz"):
        try: # load data
            with np.load(raw_data_path, allow_pickle=True) as data:
                measurements_array = data["measurements"]
                loaded_data = measurements_array.tolist()
        except Exception as e:
            print(f"Error loading file {raw_data_path}! Abort!\n{e}")
            return []
    elif raw_data_path.endswith(".h5"):
        try:
            raw_data_file = h5py.File(raw_data_path, "r")
            measurement_ids = sorted(raw_data_file.keys())

            for meas_id in measurement_ids:
                meas_group = raw_data_file[meas_id]
                measurement_data = {"measurement_number":int(meas_group.attrs["measurement_number"]),
                                    "channels": {}}
                channel_ids = sorted(meas_group.keys())

                for ch_id in channel_ids:
                    ch_group = meas_group[ch_id]
                    ch_type = ch_group.attrs["type"]

                    if isinstance(ch_type, bytes):
                        ch_type = ch_type.decode()

                    channel_data = {
                        "type": ch_type,
                        "time_data": ch_group["time_data"][:],
                        "volt_data": ch_group["volt_data"][:],
                        "vcode_per": ch_group.attrs["vcode_per"],
                        "vdiv": ch_group.attrs["vdiv"]
                    }
                    measurement_data["channels"][ch_id] = channel_data
                loaded_data.append(measurement_data)
        except Exception as e:
            print(f"Could not load file {raw_data_path}\n{e}")
            return []
        finally:
            if raw_data_file is not None:
                raw_data_file.close()
    else:
        print("ERROR: File format not recognized! Can not load file. Abort!\n")
        quit()

    if len(loaded_data) == 0:
        print("No measurements found in file.")
        return []

    channel_type_map = {}

    for measurement in loaded_data:
        meas_channels = measurement["channels"]

        for ch_key, ch_data in meas_channels.items(): # ch_key is ch0, ch1, ch2
            ch_num = int(ch_key.replace("ch", ""))
            ch_type = ch_data["type"]

            if ch_num not in channel_type_map:
                channel_type_map[ch_num] = ch_type

    offline_channels = sorted(channel_type_map.keys())

    print("Detected channels:")
    for ch in offline_channels:
        print(f"\tChannel {ch}: {channel_type_map[ch]}")

    base_path = raw_data_path
    if base_path.endswith(".npz"):
        base_path = base_path[:-4]
    elif base_path.endswith(".h5"):
        base_path = base_path[:-3]

    offline_fit_file_path = base_path + "_offline_fit.txt"

    header = ["meas_no"]

    for channel in offline_channels:
        ch_use = channel_type_map[channel]

        header.append(f"type_ch_{channel}")

        if ch_use in fit_param_names:
            par_names = fit_param_names[ch_use]
        else: print("Problem with header parameters")

        header.extend(f"{par}_ch{channel}" for par in par_names)
        header.extend(f"{par}_err_ch{channel}" for par in par_names)
        header.append(f"chi2_red_ch{channel}")

    with open(offline_fit_file_path, "w") as f:
        f.write("\t".join(header) + "\n")

    print(f"\nWriting offline fit results to:\n{offline_fit_file_path}\n")
    print(f"Starting fits for {len(loaded_data)} measurements ...")

    for measurement in loaded_data:
        meas_num = measurement["measurement_number"]
        meas_channels = measurement["channels"]

        # print(f"Processing measurement {meas_num}")

        measurement_data = {
            "measurement_number": meas_num,
            "channels": {}
        }

        write_row = [str(meas_num)]
        measurement_valid = True

        for channel in offline_channels:
            ch_key = f"ch{channel}"

            if ch_key not in meas_channels:
                print(f"\033[31mMeasurement {meas_num}: missing channel {channel}. Excluding measurement.\033[0m")
                measurement_valid = False
                broken_measurement_number += 1
                break

            ch_data = meas_channels[ch_key]
            ch_use = ch_data["type"]

            time_data = np.asarray(ch_data["time_data"], dtype=float)
            volt_data = np.asarray(ch_data["volt_data"], dtype=float)
            vcode_per = ch_data["vcode_per"]
            vdiv = ch_data["vdiv"]

            # print(f"  Channel {channel}: {ch_use}, len(data)={len(time_data)}")

            if len(time_data) == 0 or len(volt_data) == 0:
                print(f"\033[31mNo data for measurement {meas_num}, channel {channel}. Excluding measurement.\033[0m")
                measurement_valid = False
                broken_measurement_number += 1
                break

            amplitude = np.max(volt_data) - np.min(volt_data)

            diff = np.diff(volt_data)
            i_on, i_off = np.argmax(diff), np.argmin(diff)
            t_on0, t_off0 = time_data[i_on], time_data[i_off]

            if t_off0 < t_on0:
                t_on0, t_off0 = t_off0, t_on0

            if volt_data[-1] < volt_data[0]:
                amplitude = -amplitude

            if ch_use == "comp":
                p0 = [amplitude, 30.0, 0.65, 2000.0, 0.01, 0.0]
                fit_func = onoff_errorfunc
                fit_bounds = (-np.inf, np.inf)

            elif ch_use == "amp":
                p0 = [2 * amplitude, 10.0, 0.5, 100.0, 0.01, 0.5]
                fit_func = onoff_errorfunc
                fit_bounds = (-np.inf, np.inf)

            elif ch_use == "inj":
                p0 = [amplitude, 0.0, -0.1, 1.4]
                fit_func = errorfunc
                fit_bounds = (-np.inf, np.inf)

            else:
                print(f"\033[31mUnknown channel type '{ch_use}' for channel {channel} in measurement {meas_num}. Excluding measurement.\033[0m")
                measurement_valid = False
                broken_measurement_number += 1
                break

            try:
                popt, pcov = scipy.optimize.curve_fit(
                    fit_func,
                    time_data,
                    volt_data,
                    p0=p0,
                    bounds=fit_bounds,
                    maxfev=1000
                )
            except Exception as e:
                print(f"\033[31mFit failed for measurement {meas_num}, channel {channel}: {e}\033[0m")
                measurement_valid = False
                broken_measurement_number += 1
                break

            popt_err = np.sqrt(np.diag(pcov))
            fit_line = fit_func(time_data, *popt)

            yerr = np.ones_like(volt_data) / vcode_per * float(vdiv) * vertical_divisions / 2

            reduced_chi_sq = (1.0 / (len(volt_data) - len(popt))) * np.sum(
                ((volt_data - fit_line) / yerr) ** 2
            )

            measured_amplitude = np.max(volt_data) - np.min(volt_data)

            if ch_use == "comp":
                fitted_amp = abs(popt[0])
            elif ch_use == "amp":
                fitted_amp = abs(popt[0])
            elif ch_use == "inj":
                fitted_amp = 2 * abs(popt[0])
            else:
                fitted_amp = np.nan

            if (
                fitted_amp < 0.3 * measured_amplitude
                or fitted_amp > 3.0 * measured_amplitude
                or popt_err[0] > 3 * fitted_amp
                or reduced_chi_sq > chi_barrier
            ):
                print(
                    f"\033[31mFit amplitude unrealistic for measurement {meas_num}, "
                    f"channel {channel} {ch_use}: "
                    f"fitted={fitted_amp:.3f}+/-{popt_err[0]:.3f}, "
                    f"measured={measured_amplitude:.3f}, "
                    f"chi2={reduced_chi_sq:.3f}. Excluding measurement.\033[0m"
                )
                measurement_valid = False
                broken_measurement_number += 1
                break

            measurement_data["channels"][f"ch{channel}"] = {
                "type": ch_use,
                "time_data": time_data,
                "volt_data": volt_data,
                "volt_error": yerr,
                "popt": popt,
                "popt_err": popt_err,
                "chi2_red": reduced_chi_sq
            }

            write_row.append(ch_use)
            write_row.extend(str(val) for val in popt)
            write_row.extend(str(val) for val in popt_err)
            write_row.append(str(reduced_chi_sq))

        if measurement_valid:
            with open(offline_fit_file_path, "a") as f:
                f.write("\t".join(write_row) + "\n")

            stored_results_offline.append(measurement_data)
            valid_fit_count_offline += 1
            # print("\033[32mOffline fit successful.\033[0m")

        else:
            print("\033[33mMeasurement excluded in offline analysis.\033[0m")

    print("\n\033[32mOffline analysis finished.\033[0m")
    print(f"Valid fits: {valid_fit_count_offline}")
    print(f"Excluded measurements: {broken_measurement_number}")
    print(f"Fit results written to: {offline_fit_file_path}\n")

    return stored_results_offline

stop_event = threading.Event()   
def look_for_stop():
    global stop_event
    while not stop_event.is_set():
        command = input("\nEnter 'x' and press enter to stop the measurement.\n").lower()
        if command == 'x':
            print("Stopping after active measurement.")
            stop_event.set()
            break
##############################################################
'''
Start the actual measurement
'''
if __name__ == '__main__':
    chip = "chip029" # or ""  -> to be included in the plots
    common_folder = "Chip029/reference/effdacs/60V/Injection/"
    plots_folder = "Plots/" + common_folder
    data_folder = "Data/" + common_folder
    screenshot_folder = "Screenshots/" + common_folder
    file_suffix = '_chip029_inj_500mV' 
    os.makedirs(plots_folder, exist_ok=True)
    os.makedirs(data_folder, exist_ok=True)
    os.makedirs(screenshot_folder, exist_ok=True)
    start_time = time.time()
    chi_barrier = 15000

    run_loop = True # always true, set the measurement ON, If False: Analysis only (in the future)
    plot_wave = False
    plot_hist = True
    plot_corr_map = True # needs plot_hist to be activated
    plot_width_comparison = True # needs plot_hist to be activated
    take_screenshot = False
    store_raw_data = True # needs to be activated all the time as the wave forms are read from raw data
    data_taking_time = 2000 # time on which the measurement runs, if None: runs until stopped by keyboard or numbers of measurement
    data_taking_number = 1000 # Number of measurements after which the data taking ends, if None: Runs until stopped otherwise (keyboard or time)
    
    recreate_fit_offline = False
    analysis_fit_file = data_folder+"measurement_26-06-2026_12-12-38_chip026_inj_500mV.h5"


    connected_channels = [0,1,2] # channels in use, starting with 0 for the first channel(1)
    channel_use = ["inj", "amp", "comp"] # what the channels are used for in the correct order; "comp" for comparator, "amp" for amplifier, "empty" for the last channel to avoid empty measurements in the important channel

    stored_results = []
    if run_loop:
        #establishing connection to osci and setting length of eceived data
        _rm = pyvisa.ResourceManager()
        sds = _rm.open_resource(SDS_RSC)
        sds.timeout = 30000 # default value is 2000(2s)
        sds.chunk_size = 20 * 1024 * 1024 # default value is 20*1024(20k bytes)
        if data_taking_time is not None:
            print(f"Estimated End of Measurement({data_taking_time} s): {datetime.fromtimestamp(start_time+data_taking_time)}")
        if data_taking_number is not None:
            print(f"Estimated End of Measurement(No. {data_taking_number}): {datetime.fromtimestamp(start_time+(data_taking_number*1.2))}")
        '''
            Create file to store fit data and identifier
        '''
        dt_object = datetime.fromtimestamp(start_time)
        filename = dt_object.strftime("measurement_%d-%m-%Y_%H-%M-%S")
        header = ["meas_no"]
        for channel in connected_channels:
            header.append(f"type_ch_{channel}")
            header.extend(f"{par}_ch{channel}" for par in fit_param_names[channel_use[channel]])
            header.extend(f"{par}_err_ch{channel}" for par in fit_param_names[channel_use[channel]])
            header.append(f"chi2_red_ch{channel}")

        fit_file_path = data_folder + filename + file_suffix + "_fit.txt"
        valid_fit_count = 0
        with open(fit_file_path, 'w') as f:
            f.write('\t'.join(header)+'\n')
    else:
        print("Analysis Mode, no data taking")
        try:
            if recreate_fit_offline:
                print("Recreating fit to raw data ...")
                stored_results = offline_fit(analysis_fit_file)
                filename = os.path.basename(analysis_fit_file)
                if filename.endswith(".npz"):
                    filename = filename.replace(".npz","_re")
                elif filename.endswith(".h5"):
                    filename = filename.replace(".h5", "_re")
                file_suffix = ""

            else:
                fit_file_path = analysis_fit_file
                filename = os.path.basename(fit_file_path)
                filename = filename.replace("_fit.txt","_offline")
                file_suffix = ""
                plot_wave=False

                stored_results = load_fit_results_from_txt(fit_file_path, connected_channels, fit_param_names)
                valid_fit_count=len(stored_results)
                print(f"Loaded {valid_fit_count} measurements from {fit_file_path}")
                if not plot_hist:
                    print("If you want to create histograms, you need to activate this function and run again!")
        except Exception as e:
            print(f"Could not access file {analysis_fit_file} correctly. Maybe you are trying to access npz or txt file when you want the other one. Abort!\n{e}")
            quit()

    '''
        Start Loops
    '''
    loops_started = False
    if run_loop:
        use_keyboard_stop = (data_taking_number == None and data_taking_time == None)

        if use_keyboard_stop:
            keyboard_thread = threading.Thread(target=look_for_stop, daemon=True)
            keyboard_thread.start()
        loop_thread = threading.Thread(target=measurement_loop)
        loop_thread.start()
        loops_started = True

    if loops_started: loop_thread.join()

    stop_event.set()
    if run_loop:
        print("\nAll Measurements finished.")
        end_time = time.time()
        elapsed_time = end_time - start_time
        print("Total measurement duration: {:.3f} seconds".format(elapsed_time))
        connected_channels.pop()
        channel_use.pop()
        
        # Load measurements from h5 file to preserve memory during campaign
        if store_raw_data:
            h5_file_path = data_folder + filename + file_suffix + ".h5"
            print(f"Loading measurements from {h5_file_path}...")
            stored_results = load_measurements_from_h5(h5_file_path, fit_file_path, connected_channels, fit_param_names)
            print(f"Loaded {len(stored_results)} measurements from h5 file")
        else:
            # If raw data wasn't stored, load from txt file (fit results only, no raw waveforms)
            print(f"Loading fit results from {fit_file_path}...")
            stored_results = load_fit_results_from_txt(fit_file_path, connected_channels, fit_param_names)
            print(f"Loaded {len(stored_results)} measurements from txt file")
        
    else:
        elapsed_time = 0.0


    if plot_wave and len(stored_results) > 0 and len(stored_results) < 20:
        make_wave_plots(stored_results, chip)

    if plot_hist and len(stored_results) > 0:
        # if run_loop:
        #     stored_results = load_fit_results_from_txt(fit_file_path, connected_channels, fit_param_names)
        make_hist_plots(stored_results, plot_correlation=plot_corr_map, plot_pulse_width=plot_width_comparison, chip=chip)
