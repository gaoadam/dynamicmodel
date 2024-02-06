#Import libraries

#Numeric/Computational libraries
import numpy as np

#Visualization libraries
from matplotlib import pyplot as plt


def step(x, dt, f, t, args):
    """
    Description:
        step function to return vector (array) x at t + dt given x and other parameters
    Args:
        x: 1 dimensional array, shape n (numpy array)
            contains n values at time t
        dt: float; scalar value
            timestep quantity
        f: 1 dimensional list, shape n
            contains n functions to be used on x to calculate time derivatives of x
        t: float; time value, passed into functions in f
        args: dict;  arguments to determine f
        
    """
    
    #Validate that the input shapes match
    try:
        if x.shape[0] == len(f):
            #Calculate x at t+dt
            x_next = x.copy()
            for n in range(0, x.shape[0]):
                f_next = f[n]
                x_next[n] = x[n] + f_next(x, t, args)*dt

            return x_next

    except ValueError:
            print("x and f do not have matching lengths")

def x_iterate(x_t0, dt, N, f, args={}):
    """
    Description:
        generate an array containing values of vector x over time, ie a 3d array
    Args:
        x_t0: 1 dimensional array (numpy array), shape n
            contains n values at time t = 0
        dt: float; scalar value
            timestep quantity
        N: int;
            number of time steps to be iterated through
        f: 1 dimensional list, shape n
            contains n functions to be used on x in a step function, to calculate x at each time step t+dt
        args: dict; arguments to determine f
    """
    
    #Initialize tensor at time = 0
    x_current = x_t0.copy()
    
    #Initiate output data structure for x  by adding a dimension
    x_full = np.expand_dims(x_current, axis=0)
    
    #Compute x at next time step N times
    for n in range(0, N):
        #compute current time
        t = dt*n
        #Calculate x at t+dt
        x_next = step(x_current, dt, f, t, args)
        #Append x at t+dt
        x_full = np.concatenate((x_full.copy(), np.expand_dims(x_next, axis=0)))
        #Initialize current tensor for next iteration
        x_current = x_next.copy()
        
    #add time as 1st column
    x_full = np.concatenate((np.expand_dims(np.linspace(0, N*dt, N+1), axis=1), x_full), axis=1)
    
    #return output value
    return x_full

def harmonics(t, args):
    #Get arguments
    a_list, n_list = args['a_list'],args['n_list']
    #Initialize output value
    t_out = 0
    #iterate through harmonic values
    for n in range(0, len(n_list)):
        t_out = t_out + a_list[n]*np.sin(n_list[n]*2*np.pi*t)
    return t_out

def x_driven(x_t0, dt, N, m, b, k, u1, args):
    """
    Description:
        Simulate damped driven oscillator, with driving force u_1(t)
    x_t0: 1 dimensional array, shape 2
            contains 2 values at time t = 0: position and velocity
    dt: scalar value
            timestep quantity
    N: scalar value
        number of time steps to be iterated through
    m: scalar value
        mass constant
    b: scalar value
        friction constant
    k: scalar value
        spring constant
    u: 1 dimensional list, shape n
        contains n functions to be used on x; input functions in a system
    args: dictionary
        arguments for function u
    """
    #Prepare system functions
    #Velocity
    def xdot(x, t, args):
        return x[1]
    #Acceleration
    def xdotdot(x, t, args):
        return (u1(t, args) - b*x[1] - k*x[0])/m
    
    #Iterate through the time steps to calculate the variables using the system functions
    x_full = x_iterate(x_t0=x_t0, dt=dt, N=N, f=[xdot, xdotdot], args=args)
    
    return x_full

def x_rlc(x_t0, dt, N, r, l, c, u1, args):
    """
    Description:
        Simulate RLC circuit, with voltage u_1(t)
    x_t0: 1 dimensional array, shape 2: capacitor voltage and inductor current
            contains 2 values at time t = 0: 
    dt: scalar value
            timestep quantity
    N: scalar value
        number of time steps to be iterated through
    r: scalar value
        resistance constant
    l: scalar value
        inductor constant
    c: scalar value
        capacitance constant
    u: 1 dimensional list, shape n
        contains n functions to be used on x; input functions in a system
    args: dictionary
        arguments for function u
    """
    #Prepare system functions
    #function for v_c ie capacitor voltage
    def x1dot(x, t, args):
        return x[1]/c
    
    #function for i ie inductor current
    def x2dot(x, t, args):
        return -1/l*x[0] - r/l*x[1] + 1/l*u1(t, args)
    
    #Iterate through the time steps to calculate the variables using the system functions
    x_full = x_iterate(x_t0=x_t0, dt=dt, N=N, f=[x1dot, x2dot], args=args)
    
    return x_full

def plot_df(df, title):
    """
    Description:
        Take dataframe of specific specifications, convert to torch tensor, and plot values.
        Time 't' should ALWAYS be the first column of the dataframe.

    Args:
        df: pandas Dataframe, column values are the signal values over time
        columns: list of str; columns to plot
        title: str, plot title
    """
    #convert to numpy array
    x_array = df.to_numpy()
    
    #Get array shape to determine plot quantity
    n = np.shape(x_array)[1]
    
    #initialize subplots
    fig, axs = plt.subplots(n-1, constrained_layout=True)
    fig.suptitle(title)
    
    #plot each variable with respect to t
    if n-1>1:
        for i in range(1,n):
        #plot x(t)
            axs[i-1].plot(x_array[:,0], x_array[:,i])
            axs[i-1].set_xlabel('t')
            axs[i-1].set_ylabel(df.columns[i])
    elif n-1==1:
        axs.plot(x_array[:,0], x_array[:,1])
        axs.set_xlabel('t')
        axs.set_ylabel(df.columns[1])
def plot_x(x_array, columns, title):
    """
    Description:
        Take numpy array, and plot values.
        Time 't' should ALWAYS be the first index and dimension 1 (ie the 2nd dimension).

    Args:
        x_array: numpy array
        columns: list of str; column names to plot
        title: str; plot title
    """
    #Get array shape to determine plot quantity
    n = np.shape(x_array)[1]
    
    #initialize subplots
    fig, axs = plt.subplots(n-1,  constrained_layout=True)
    fig.suptitle(title)
    
    #plot each variable with respect to t
    if n-1>1:
        for i in range(1,n):
        #plot x(t)
            axs[i-1].plot(x_array[:,0], x_array[:,i])
            axs[i-1].set_xlabel('t')
            axs[i-1].set_ylabel(df.columns[i])
    elif n-1==1:
        axs.plot(x_array[:,0], x_array[:,1])
        axs.set_xlabel('t')
        axs.set_ylabel(columns[1])