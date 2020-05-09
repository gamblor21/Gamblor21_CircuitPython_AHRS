##
## My own test file, quick and dirty to see if the algorithm worked
## Only calibrates for the gryo and hardiron offsets
## Set to run on the Adafruit LSM9DS1 over I2C cause that is what I have
##

import mahony

import board
import busio
import time
import adafruit_lsm9ds1

MAG_MIN = [-0.5764, 0.0097, -0.5362]
MAG_MAX = [0.4725, 0.9919, 0.4743]


def map_range(x, in_min, in_max, out_min, out_max):
    """
    Maps a number from one range to another.
    :return: Returns value mapped to new range
    :rtype: float
    """
    mapped = (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min
    if out_min <= out_max:
        return max(min(mapped, out_max), out_min)

    return min(max(mapped, out_max), out_min)


m = mahony.mahony()
m.twoKp = 2 * 25
m.twoKi = 2 * 0
print(m.twoKp)
print(m.twoKi)

sensor = adafruit_lsm9ds1.LSM9DS1_I2C(board.I2C())
# first values always seem off so discard them and wait a second
mx, my, mz = sensor.magnetic
gx, gy, gz = sensor.gyro
ax, ay, az = sensor.acceleration
time.sleep(0.5)

count = 0
lastPrint = time.monotonic()
timestamp = time.monotonic_ns()

while True:
    if (time.monotonic_ns() - timestamp) > 6500000:
        mx, my, mz = sensor.magnetic
        # adjust for magnetic calibration - hardiron only
        mx = map_range(mx, MAG_MIN[0], MAG_MAX[0], -1, 1)
        my = map_range(my, MAG_MIN[1], MAG_MAX[1], -1, 1)
        mz = map_range(mz, MAG_MIN[2], MAG_MAX[2], -1, 1)

        gx, gy, gz = sensor.gyro
        # adjust for my gyro calibration values
        gx -= 1.1250
        gy -= 3.8732
        gz += 1.2834

        ax, ay, az = sensor.acceleration

        m.update(gx, gy, gz, ax, ay, az, mx, my, mz)

        count += 1
        timestamp = time.monotonic_ns()

    if time.monotonic() > lastPrint + 0.1:
        m.computeAngles()
        print(
            "Orientation: ",
            m.yaw * 57.29578 + 180,
            ", ",
            m.pitch * 57.29578,
            ", ",
            m.roll * 57.29578,
        )
        print("Quaternion: ", m.q0, ", ", m.q1, ", ", m.q2, ", ", m.q3)
        # print("Count: ", count)
        lastPrint = time.monotonic()
        count = 0
