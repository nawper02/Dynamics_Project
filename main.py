import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from rk4 import rk4
from Ball import Ball


# TODO: Check analytical solution / derivation for rubber ball


# Define system of ODE's
def func_block(t, y, ball):

    s = y[0]
    v = y[1]
    a = y[2] # last acceleration, used to check slip condition

    theta = s / R
    Fn_tot = ((ball.m / (R-ball.d)) * pow(v, 2) + ball.m * 9.81 * np.sin(theta)) / np.sin(ball.gamma)

    dydt = np.zeros(len(y))
    dydt[0] = v

    # if the ball is slipping, the acceleration is

    if (((2/5) * ball.m * pow(ball.r, 2) * a)/pow(ball.d, 2)) > ball.mu_s * Fn_tot: # TODO: CHECK THIS
        #ball.m * 9.81 * np.cos(theta) > ball.mu_s * Fn_tot: ?
        dydt[1] = 9.81 * np.cos(theta) - ball.mu_k * (((1/(R-ball.d)) * pow(v, 2) + 9.81 * np.sin(theta)) / np.sin(ball.gamma))

    # if the ball is rolling, the acceleration is
    else:
        dydt[1] = (9.81 * pow(ball.d, 2) * np.cos(theta)) / (((2.0/5.0) * pow(ball.r, 2)) + pow(ball.d, 2))

    #dydt[2] = dydt[1] # This is last acceleration

    return dydt


# A function that computes if the simulation passes the top of the loop with non-negative normal force
def passes_tests(H, ball):
    # Compute various lists using H
    y0 = compute_initial_conditions(H, ball, step)

    # Solve ODE
    t_res, y = rk4(func_block, t, y0, step, ball)
    pos = y[:, 0]  # pos [m]
    vel = y[:, 1]  # velocity [m/s]

    ang, fNorm, angDeg = extrapolate_data(pos, vel, ball)

    # Loop over angDeg array
    for index, element in enumerate(angDeg):
        # When an element is found in this range,
        if 260 < element < 280:
            # Check to see if the normal force there is greater than zero
            if fNorm[index] > 0:
                print(f"fNorm at top: {fNorm[index]}")
                print(f"H = {H / 0.0254}")
                return True
    else:
        return False


def compute_initial_conditions(H, ball, dt):
    # Compute various lists using H
    theta = np.radians(angRamp)

    #hLoop = H - R * (1 - np.sin(angInitRad))
    #sRamp = hLoop / np.cos(angInitRad)
    sLoop = R * angInitRad
    vLoop, aLoop = compute_vLoop(ball, theta, H, dt)

    y0 = [sLoop, vLoop, aLoop]
    return y0


def extrapolate_data(pos, vel, ball):

    ang = pos / R  # loop angle [rad]
    angDeg = (ang * 180) / np.pi  # loop angle (90 deg is bottom of loop) [deg]

    fNorm_tot = ((ball.m / (R-ball.d)) * pow(vel, 2) + ball.m * 9.81 * np.sin(ang)) / np.sin(ball.gamma)

    return ang, fNorm_tot, angDeg


def compute_drop_height(ball, h_guess):
    # Loop to increment H and check if the tests are passed
    while True:
        h_guess += 0.0001
        print(h_guess/0.0254)
        tests_passed = passes_tests(h_guess, ball)
        if tests_passed:
            print(f"Drop height found! It is {h_guess/0.0254} inches")
            return h_guess


def compute_vLoop(ball, theta, H, dt):
    """
    :param dt:
    :param ball: ball object
    :param theta: angle of ramp (rad)
    :param H: drop height (m)
    :return: velocity at loop entry
    """

    l = (H - R + (R * np.cos(theta))) / np.sin(theta)

    sg = 0
    vg = 0
    ag = 0

    Fn_tot = ball.m * 9.81 * np.cos(theta) / np.sin(ball.gamma)

    while sg < l:
        # If the ball is slipping
        if (((2/5) * ball.m * pow(ball.r, 2) * ag)/pow(ball.d, 2)) > ball.mu_s * Fn_tot:
            # Old condition: 0.5 * (m * 9.81 * np.sin(theta) - m * ag) > ((mu_s * m * 9.81) / (2 * np.sin(y)))
            ag = 9.81 * (np.sin(theta) - (ball.mu_k * np.cos(theta)/np.sin(ball.gamma)))
        # If the ball is rolling
        else:
            ag = (9.81 * np.sin(theta)) / (1 + ((2 * np.power(ball.r, 2)) / (5 * np.power(ball.d, 2))))
        vg += ag * dt
        sg += vg * dt
    return vg, ag


def compute_rubber_ball(ball):
    H = R * (2 + 0.2 * ball.r + ((pow(ball.d, 2)) / (2 * ball.r)))
    print(f"Height ratio for rubber ball (Should be close to 2.7): {H/R}")
    # TODO: Check. should be close to 2.7R
    return H


def run_sim_and_plot(H, ball, name):
    y0 = compute_initial_conditions(H, ball, step)
    print(y0)

    # Solve ODE
    t_res, y = rk4(func_block, t, y0, step, ball)
    pos = y[:, 0]  # pos [m]
    vel = y[:, 1]  # velocity [m/s]

    # Extrapolate data
    ang, fNorm, angDeg = extrapolate_data(pos, vel, ball)

    # Create First Plot

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.plot(t, angDeg, '-b', label='Angle (deg)')
    ax2.plot(t, vel, '-r', label='Velocity')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Angle (deg)')
    ax2.set_ylabel('Velocity')
    plt.title(f'Simulation of {name} -- angle and velocity vs time')
    fig.legend()
    plt.savefig(f"{name}_ang_vel.pdf")

    # Create Second Plot

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.plot(angDeg, fNorm, '-b', label='Normal Force')
    ax2.plot(angDeg, vel, '-r', label='Velocity')
    plt.title(f'Simulation of {name} -- normal force and velocity vs angle')
    ax1.set_xlabel('Angle (deg)')
    ax1.set_ylabel("Normal Force")
    ax2.set_ylabel('Velocity')
    fig.legend()
    plt.savefig(f"{name}_fnorm_vel.pdf")

    plt.show()


def main():

    # Find min H for rubber ball
    H_rubber = compute_rubber_ball(rubber)

    # Use rubber ball drop height as initial guess
    h_guess = H_rubber

    # Compute drop height for steel ball
    H_stainless_steel = h_guess #compute_drop_height(stainless_steel, h_guess)

    # Compute drop height for plastic ball
    H_plastic = h_guess #compute_drop_height(plastic, h_guess)

    # Simulate and plot all 3 balls
    run_sim_and_plot(H_rubber, rubber, "rubber")
    run_sim_and_plot(H_plastic, plastic, "plastic")
    run_sim_and_plot(H_stainless_steel, stainless_steel, "stainless_steel")
    """
    Height ratio for rubber ball (Should be close to 2.7): 2.0033839300411524
    [0.08866272600131193, 1.8657995712490034, 6.296416183003713]
    [0.08866272600131193, 1.8633017911431131, 6.296416183003712]
    [0.08866272600131193, 1.8594643041348715, 6.296416183003713]
    """


# Entry point
if __name__ == "__main__":

    # Initialize global variables
    R = 5   * 0.0254            # IN TO M
    angRamp = 40                # DEG
    #mu = 0.213
    step = 0.01
    t = np.arange(0, 2, step)    # S
    angInit = 90 - angRamp
    angInitRad = np.radians(angInit)

    # Create ball objects
    stainless_steel = Ball(14e-3/2, 11e-3)
    plastic = Ball(15e-3/2, 3e-3)
    rubber = Ball(15.9e-3/2, 2.8e-3)

    main()
