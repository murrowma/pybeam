import numpy as np

def subject_rt(subject_number, data, subjects, min_rt, max_rt):

    # determine how many data points are in the file
    N_data = len(data[:,3])
   
    # preallocate arrays for each coherence condition, one array each for the upper and lower threhsolds
    rt0_upper = np.zeros(N_data)
    rt0_lower = np.zeros(N_data)
    rt5_upper = np.zeros(N_data)
    rt5_lower = np.zeros(N_data)
    rt10_upper = np.zeros(N_data)
    rt10_lower = np.zeros(N_data)
    rt40_upper = np.zeros(N_data)
    rt40_lower = np.zeros(N_data)

    # load data for the subject into the arrays
    for ii in range(N_data):
        if (data[ii,0] == int(subject_number)):

            if (data[ii,2] == 1) & (data[ii,3] < max_rt) & (data[ii,3] > min_rt):
                if (data[ii,1] == 0):
                    rt0_upper[ii] = data[ii,3]
                if (data[ii,1] == 5):
                    rt5_upper[ii] = data[ii,3]
                if (data[ii,1] == 10):
                    rt10_upper[ii] = data[ii,3]
                if (data[ii,1] == 40):
                    rt40_upper[ii] = data[ii,3]

            if (data[ii,2] == 0) & (data[ii,3] < max_rt) & (data[ii,3] > min_rt):
                if (data[ii,1] == 0):
                    rt0_lower[ii] = data[ii,3]
                if (data[ii,1] == 5):
                    rt5_lower[ii] = data[ii,3]
                if (data[ii,1] == 10):
                    rt10_lower[ii] = data[ii,3]
                if (data[ii,1] == 40):
                    rt40_lower[ii] = data[ii,3]

    # place rt data into dictionaries of the right for to be used with pybeam
    rt0 = { 'rt_upper' : rt0_upper[rt0_upper != 0] , 'rt_lower' : rt0_lower[rt0_lower != 0] }
    rt5 = { 'rt_upper' : rt5_upper[rt5_upper != 0] , 'rt_lower' : rt5_lower[rt5_lower != 0] }
    rt10 = { 'rt_upper' : rt10_upper[rt10_upper != 0] , 'rt_lower' : rt10_lower[rt10_lower != 0] }
    rt40 = { 'rt_upper' : rt40_upper[rt40_upper != 0] , 'rt_lower' : rt40_lower[rt40_lower != 0] }
    
    # return a dictionary containing these dictionaries
    rt = {'rt0' : rt0, 'rt5' : rt5, 'rt10' : rt10, 'rt40' : rt40}
    return rt
