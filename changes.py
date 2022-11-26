import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from rk4 import rk4
from ball import Ball


# TODO: make work for mu = >1
# create a class that describes a rolling ball
import numpy as np


class Ball():
    def __init__(self, radius, mass):
        self.mass = mass
        self.radius = radius
        self.gamma = np.arcsin(np.sqrt(np.power((self.radius+1.5), 2)+49)/(self.radius+1.5))
        self.d = self.radius * np.sin(self.gamma)
        self.mu_s = np.sin(np.radians(12)) * np.sin(self.gamma)
        self.mu_k = self.mu_s / 2.0 # GUESS

# Define system of ODE's
def func_block(t, y):
    """
    :param t: time
    :param y: current state vector, y = [s, v]
    :returns: derivative of state vector, dydt = [v, a]
    """

    s = y[0]
    v = y[1]

    theta = s / R
    Fn = mass * (9.81 * np.sin(theta) + ((v ** 2) / R))  # normal force

    dydt = np.zeros(len(y))
    dydt[0] = y[1]

    # acceleration has two modes, rolling without slipping, and slipping

    # q: when does rolling without slipping occur?
    # a: when the normal force is greater than the friction force

    # if the ball is slipping, the acceleration is
    if (5/7) * mass * 9.81 * np.cos(theta) > mu * Fn: #mass * 9.81 * np.cos(theta) > mu * Fn:
        dydt[1] = (9.81 * np.cos(theta)) - (mu * 9.81 * np.sin(theta) * np.sign(v)) - ((mu * (v ** 2)) / (R * np.sign(v)))

    # if the ball is rolling, the acceleration is
    else:
        # acceleration of a ball inside of a loop if its rolling
        dydt[1] = (5/7) * 9.81 * np.cos(theta)


    return dydt


# A function that computes if the simulation passes the top of the loop with non-negative normal force
def passes_tests(H):
    # Compute various lists using H
    hLoop = H - R * (1 - np.sin(angInitRad))
    sRamp = hLoop / np.cos(angInitRad)
    vLoop = np.sqrt(2 * 9.81 * hLoop * (1 - mu * np.tan(angInitRad)))
    sLoop = R * angInitRad
    y0 = [sLoop, vLoop]

    # Solve ODE
    y = odeint(func_block, y0, t, tfirst=True)
    pos = y[:, 0]    # pos [m]
    vel = y[:, 1]    # velocity [m/s]

    ang = pos/R     # loop angle [rad]
    fNorm = mass * (9.81 * np.sin(ang) + ((vel**2) / R))    # normal force
    angDeg = (ang*180) / np.pi  # loop angle (90 deg is bottom of loop) [deg]

    # Loop over angDeg array
    for index, element in enumerate(angDeg):
        # When an element is found in this range,
        if 260 < element < 280:
            # Check to see if the normal force there is greater than zero
            if fNorm[index] > 0:
                print(f"fNorm at top: {fNorm[index]}")
                print(f"H = {H/0.0254}")
                return True
    else:
        return False


def vLoop(ball, theta, H, dt):
    """
    :param ball: ball object
    :param theta: angle of ramp (rad)
    :param H: drop height (m)
    :return: velocity at loop entry
    """
    d = ball.d
    rb = ball.r
    m = ball.m
    mu_s = ball.mu_s
    mu_k = ball.mu_k
    y = ball.gamma
    l = (H - R + (R * np.cos(theta))) / np.sin(theta)
    sg = 0
    vg = 0
    ag = 0
    while sg < l:
        if 0.5 * (m * 9.81 * np.sin(theta) - m * ag) > ((mu_s * m * 9.81) / (2 * np.sin(y))):
            ag = 9.81 * (np.sin(theta) - (mu_k/np.sin(y)))
        else:
            ag = (9.81 * np.sin(theta)) / (1 + ((2 * m * np.power(rb, 2)) / (5 * np.power(d, 2))))
            # TODO: Check derivation. hw13 19.5. 5/7 g sin theta?
        vg += ag * dt
        sg += vg * dt
    return vg


def compute_rubber_ball(ball):
    H = R * (2 + 0.2 * ball.radius + ((pow(ball.d, 2)) / (2 * ball.radius)))

    # TODO: Check. should be close to 2.7R
    return H


def compute_drop_height(ball, h_guess):
    H = h_guess   # initial guess

    # Loop to increment H and check if the tests are passed
    while True:
        H += 0.0001
        print(H/0.0254)
        tests_passed = passes_tests(H)
        if tests_passed:
            print(f"Drop height found! It is {H/0.0254} inches")
            break


def main():
    H_rubber = compute_rubber_ball(rubber)

    h_guess = H_rubber

    H_rubber = compute_drop_height(rubber, h_guess)

    H_plastic = compute_drop_height(plastic, h_guess)


    # Perform final solving of the system for plots

    # Compute various lists using H
    hLoop = H - R * (1 - np.sin(angInitRad))

    vloop = vLoop(plastic, angInitRad, H, 0.01)

    sLoop = R * angInitRad
    y0 = [sLoop, vLoop]

    # Solve ODE one last time (we're gonna celebrate)
    y = odeint(func_block, y0, t, tfirst=True)
    pos = y[:, 0]    # pos [m]
    vel = y[:, 1]    # velocity [m/s]

    ang = pos/R     # loop angle [rad]
    fNorm = mass * (9.81 * np.sin(ang) + ((vel**2) / R))    # normal force
    angDeg = (ang*180) / np.pi  # loop angle (90 deg is bottom of loop) [deg]

    # Create First Plot

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.plot(t, angDeg, '-b', label='Angle (deg)')
    ax2.plot(t, vel, '-r', label='Velocity')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Angle (deg)')
    ax2.set_ylabel('Velocity')
    plt.title('Simulation of block -- angle and velocity vs time')
    fig.legend()
    plt.savefig("ang_vel.pdf")

    # Create Second Plot

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.plot(angDeg, fNorm, '-b', label='Normal Force')
    ax2.plot(angDeg, vel, '-r', label='Velocity')
    plt.title('Simulation of block -- normal force and velocity vs angle')
    ax1.set_xlabel('Angle (deg)')
    ax1.set_ylabel("Normal Force")
    ax2.set_ylabel('Velocity')
    fig.legend()
    plt.savefig("fnorm_vel.pdf")

    plt.show()


# Entry point
if __name__ == "__main__":
    # Stuff that doesn't need to be re-computed
    R = 5   * 0.0254            # IN TO M
    angRamp = 50                # DEG
    mu = 0.213
    mass = 1
    t = np.arange(0, 2, 0.01)    # S
    angInit = 90 - angRamp
    angInitRad = angInit * np.pi / 180
    rad_ball = 0.5 * 0.0254      # IN TO M

    # Create ball objects
    stainless_steel = Ball(14e-3, 11e-3)
    plastic = Ball(15e-5, 3e-3)
    rubber = Ball(2.8e-3, 15.9e-3)


    main()
