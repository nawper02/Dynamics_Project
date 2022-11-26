import numpy as np

def rk4(odefun, tspan, y0, h, params):
    """Uses RK4 Method to calculate the solution to an ODE.
        Parameters
        ----------
        odefun : callable
            A callable function to the derivative function defining the system.
        tspan : array_like
            An array of times (or other independent variable) at which to evaluate
            Euler's Method.  The nTimes values in this array determine the size of
            the outputs.
        y0 : array_like
            Array containing the initial conditions to be used in evaluating the
            odefun.  Must be the same size as that expected for the
            second input of odefun.
        h : step size
        params : array_like
            Array containing the parameters to be used in evaluating the
            odefun.  Must be the same size as that expected for the
            third input of odefun.
        Returns
        -------
        t : ndarray
            t[i] is the ith time (or other independent variable) at which
            the RK4 Method was evaluated.
        y : ndarray
            y[i] is the ith dependent variable at which the RK4 Method was
            evaluated.
        Notes
        -----
        I'm very proud of this function as it should be able to handle any number of equations!
        """

    tspan = np.asarray(tspan)

    # Determine the number of items in outputs
    num_times = len(tspan)
    num_equations = len(y0)

    # Initialize outputs
    t = np.zeros(num_times, dtype=float)
    y = np.zeros((num_times, num_equations), dtype=float)

    # Assign first row of outputs
    t[0] = tspan[0]
    y[0] = y0

    # Assign other rows of output
    for n in range(num_times-1):

        # Create empty vectors for k values
        # Each entry of k is the k value for that equation (IE, k1[0] is k1 for first equation)
        k1 = np.zeros(num_equations, dtype=float)
        k2 = np.zeros(num_equations, dtype=float)
        k3 = np.zeros(num_equations, dtype=float)
        k4 = np.zeros(num_equations, dtype=float)

        for j in range(num_equations-1):
            # Compute k1 for all equations
            dydt = odefun(t[n], y[n], params)
            k1[j] = h * dydt[j]

        for j in range(num_equations-1):
            # Compute k2 for all equations
            input_y = y[n].copy()
            for index, element in enumerate(y[n]):
                input_y[index] += 0.5 * k1[index]
            dydt = odefun(t[n] + (h/2), input_y, params)
            k2[j] = h * dydt[j]

        for j in range(num_equations-1):
            # Compute k3 for all equations
            input_y = y[n].copy()
            for index, element in enumerate(y[n]):
                input_y[index] += 0.5 * k2[index]
            dydt = odefun(t[n] + (h/2), input_y, params)
            k3[j] = h * dydt[j]

        for j in range(num_equations-1):
            # Compute k4 for all equations
            input_y = y[n].copy()
            for index, element in enumerate(y[n]):
                input_y[index] += k3[index]
            dydt = odefun(t[n] + h, input_y, params)
            k4[j] = h * dydt[j]

        # Calculate the next state
        t[n + 1] = t[n] + h
        y[n+1, 2] = y[n, 1] # Save last acceleration for func_block slip condition
        for i in range(num_equations-1):
            y[n+1, i] = y[n, i] + (1.0/6.0) * (k1[i] + (2 * k2[i]) + (2 * k3[i]) + k4[i])

        # Set last acceleration for this state to the acceleration from the last state

    return t, y
