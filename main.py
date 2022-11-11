import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define system of ODE's
def func_block(t, y):
    s = y[0]
    v = y[1]
    theta = s / R
    dydt = np.zeros(len(y))
    dydt[0] = y[1]
    dydt[1] = (9.81 * np.cos(theta)) - (mu * 9.81 * np.sin(theta) * np.sign(v)) - ((mu * (v ** 2)) / (R * np.sign(v)))
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

    ang = pos / R  # loop angle [rad]
    fNorm = mass * ((9.81 * np.sin(ang) + (vel ** 2)) / R)  # normal force
    angDeg = (ang * 180) / np.pi  # loop angle (90 deg is bottom of loop) [deg]

    # Loop over angDeg array
    for index, element in enumerate(angDeg):
        # When an element is found in this range,
        if 269.5 < element < 270.5:
            print("Passes top")
            # Check to see if the normal force there is greater than zero
            if fNorm[index] > -2:
                print(f"fNorm at top: {fNorm[index]}")
                print(f"H = {H/0.0254}")
                print("passes top with non-negative normal force")
                return True
    else:
        return False


def main():
    H = 23.8 * 0.0254   # initial guess

    # Loop to increment H and check if the tests are passed
    #while True:
    #    H += 0.001
        #print(H/0.0254)
    #    tests_passed = passes_tests(H)
    #    if tests_passed:
    #        print(f"Drop height found! It is {H/0.0254} inches")
    #        break

    # Perform final solving of the system

    # Compute various lists using H
    hLoop = H - R * (1 - np.sin(angInitRad))
    vLoop = np.sqrt(2 * 9.81 * hLoop * (1 - mu * np.tan(angInitRad)))
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
    mu = 0.1
    mass = 1
    t = np.arange(0, 2, 0.01)    # S
    angInit = 90 - angRamp
    angInitRad = angInit * np.pi / 180

    main()
