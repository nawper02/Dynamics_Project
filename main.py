import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from rk4 import rk4
from Ball import Ball


# TODO: Check slipping conditions in loop
# TODO: Check derivation of Fntot in loop
# TODO: Check roll / slip eq's
# TODO: Check mu_s computation -- account for different normal forces?

# TODO: Possibly keep track of omega? Compute work done by friction and check cons nrg?


# Define system of ODE's
def func_block(t, y, ball):

    slipped = False

    s = y[0]
    v = y[1]
    a = y[2] # last acceleration, used to check slip condition

    theta = s / R
    Fn_tot = ((ball.m / (R-ball.d)) * pow(v, 2) + ball.m * 9.81 * np.sin(theta)) / np.sin(ball.gamma)

    dydt = np.zeros(len(y))
    dydt[0] = v

    # If ball is rubber, never slips
    if ball.name == "rubber":
        dydt[1] = (9.81 * pow(ball.d, 2) * np.cos(theta)) / (((2.0/5.0) * pow(ball.r, 2)) + pow(ball.d, 2))

    else: # if the ball is not rubber
        # if the ball is slipping, the acceleration is
        if (((2/5) * ball.m * pow(ball.r, 2) * a)/pow(ball.d, 2)) > ball.mu_s * Fn_tot: # TODO: CHECK THIS
            #ball.m * 9.81 * np.cos(theta) > ball.mu_s * Fn_tot: ?
            dydt[1] = 9.81 * np.cos(theta) - ball.mu_k * (((1/(R-ball.d)) * pow(v, 2) + 9.81 * np.sin(theta)) / np.sin(ball.gamma))
            slipped = True

        # if the ball is rolling, the acceleration is
        else:
            dydt[1] = (9.81 * pow(ball.d, 2) * np.cos(theta)) / (((2.0/5.0) * pow(ball.r, 2)) + pow(ball.d, 2))

    return dydt, slipped


# A function that computes if the simulation passes the top of the loop with non-negative normal force
def passes_tests(H, ball):
    # Compute various lists using H
    y0 = compute_initial_conditions(H, ball, step)

    # Solve ODE
    t_res, y, slipped_arr = rk4(func_block, t, y0, step, ball)
    pos = y[:, 0]  # pos [m]
    vel = y[:, 1]  # velocity [m/s]

    ang, fNorm, angDeg = extrapolate_data(pos, vel, ball)

    # Loop over angDeg array
    for index, element in enumerate(angDeg):
        # When an element is found in this range,
        if 268 < element < 272:
            # Check to see if the normal force there is greater than zero
            #if fNorm[index] > 0:
            if fNorm[index] > 0:
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
    i = 0
    while True:
        i += 1
        h_guess += 0.0001
        if i % 10 == 0: # Print every 10 iterations
            # print name of ball
            print(f"Ball: {ball.name}")
            print(f"\tCurrent guess: {h_guess/0.0254} inches")
        tests_passed = passes_tests(h_guess, ball)
        if tests_passed:
            print(f"Drop height found for {ball.name} ball! It is {h_guess/0.0254} inches")
            return h_guess


def compute_drop_height_bisection(ball, a, b):
    i = 0
    while True:
        i += 1

        c = (a + b) / 2
        tp_c = passes_tests(c, ball)
        if tp_c:
            b = c
        if not tp_c:
            a = c
        print(f"Current guess: {c/0.0254} inches")

        if abs(a - b) < 0.0001:
            print(f"Drop height found for {ball.name} ball! It is {c/0.0254} inches")
            return c


def compute_vLoop(ball, theta, H, dt):
    """
    :param dt:
    :param ball: ball object
    :param theta: angle of ramp (rad)
    :param H: drop height (m)
    :return: velocity at loop entry
    """

    l = (H - (R - (R * np.sin(np.radians(90) - theta)))) / np.sin(theta)

    sg = 0
    vg = 0
    ag = 0

    Fn_tot = ball.m * 9.81 * np.cos(theta) / np.sin(ball.gamma)

    while sg < l:
        # if the ball is rubber, it never slips
        if ball.name == "rubber":
            ag = (9.81 * np.sin(theta) * ball.d * ball.r) / (((2/5) * pow(ball.r, 2)) + pow(ball.d, 2))
        else: # if the ball is not rubber
            # If the ball is slipping
            if ball.mu_s * Fn_tot < (2 * ball.m * ball.r * ag) / (5 * ball.d):
                ag = 9.81 * (np.sin(theta) - (ball.mu_k * np.cos(theta) / np.sin(ball.gamma)))
            # If the ball is rolling
            else:
                ag = (9.81 * np.sin(theta) * ball.d * ball.r) / (((2/5) * pow(ball.r, 2)) + pow(ball.d, 2))
        vg += ag * dt
        sg += vg * dt
    return vg, ag


def compute_rubber_ball(ball):
    H = (2 * R - ball.d) + (pow(ball.r, 2) / (10 * pow(ball.d, 2))) * (R - ball.d)
    print(f"Height ratio for rubber ball (Should be close to 2.7): {H/R}")
    return H


def run_sim_and_plot(H, ball, name):
    y0 = compute_initial_conditions(H, ball, step)

    # Solve ODE
    t_res, y, slipped_arr = rk4(func_block, t, y0, step, ball)
    pos = y[:, 0]  # pos [m]
    vel = y[:, 1]  # velocity [m/s]

    # Extrapolate data
    ang, fNorm, angDeg = extrapolate_data(pos, vel, ball)

    # Create First Plot

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.plot(t, angDeg, '-b', label='Angle (deg)')
    ax2.plot(t, vel, '-r', label='Velocity')
    ax2.plot(t, slipped_arr, '-g', label='Slip')
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
    """
    # Find min H for rubber ball
    H_rubber = compute_rubber_ball(rubber)
    print(f"Min H for rubber ball (analytical): {H_rubber/0.0254} inches")

    # Use rubber ball drop height as initial guess
    h_guess = H_rubber
    """

    # Compute analytical drop height for rubber ball
    H_rubber_analytical = compute_rubber_ball(rubber)
    print(f"Rubber ball analytical drop height: {H_rubber_analytical/0.0254}")
    # print rubber ball.r
    print(f"R: {R}")
    # print rubber ball.d
    print(f"Rubber ball d: {rubber.d}")

    # Compute drop height for rubber ball (Should be 2.7 * R), last time got 14.514
    H_rubber = compute_drop_height_bisection(rubber, 0.1, 20.0)

    # Compute drop height for steel ball
    H_stainless_steel = compute_drop_height_bisection(stainless_steel, 0.1*0.0254, 20.0*0.0254)

    # Compute drop height for plastic ball
    H_plastic = compute_drop_height_bisection(plastic, 0.1*0.0254, 20.0*0.0254)

    # Print all 3 drop heights
    print(f"\nDrop height for rubber ball: {H_rubber/0.0254} inches")
    print(f"Drop height for stainless steel ball: {H_stainless_steel/0.0254} inches")
    print(f"Drop height for plastic ball: {H_plastic/0.0254} inches")

    # Simulate and plot all 3 balls
    run_sim_and_plot(H_rubber, rubber, "rubber")
    run_sim_and_plot(H_plastic, plastic, "plastic")
    run_sim_and_plot(H_stainless_steel, stainless_steel, "stainless_steel")


# Entry point
if __name__ == "__main__":

    # Initialize global variables
    R = 5   * 0.0254            # IN TO M
    angRamp = 50                # DEG
    #mu = 0.213
    step = 0.0005
    t = np.arange(0, 1, step)    # S
    angInit = 90 - angRamp
    angInitRad = np.radians(angInit)

    # Create ball objects
    stainless_steel = Ball(14e-3/2, 11e-3, "stainless")
    plastic = Ball(15e-3/2, 3e-3, "plastic")
    rubber = Ball(15.9e-3/2, 2.8e-3, "rubber")

    main()
