# The MIT License (MIT)
#
# Copyright (c) 2020 Mark Komus
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""
`ahrs`
================================================================================

AHRS library for CircuitPython
Mahony Algorithm

Madgwick's implementation of Mayhony's AHRS algorithm.
See: http:##www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/

From the x-io website "Open-source resources available on this website are
provided under the GNU General Public Licence unless an alternative licence
is provided in source."

Original Information
Date			Author			Notes
29/09/2011	SOH Madgwick    Initial release
02/10/2011	SOH Madgwick	Optimised for reduced CPU load
Algorithm paper:
http:##ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=4608934&url=http%3A%2F%2Fieeexplore.ieee.org%2Fstamp%2Fstamp.jsp%3Ftp%3D%26arnumber%3D4608934

This version based upon AdaFruit AHRS https://github.com/adafruit/Adafruit_AHRS

* Author(s): Mark Komus

Implementation Notes
--------------------

**Hardware:**

    `Any 9DOF sensor`_

**Software and Dependencies:**

* Adafruit CircuitPython firmware for the supported boards:
  https://github.com/adafruit/circuitpython/releases

"""

# imports
import math

__version__ = "0.0.0-auto.0"
__repo__ = "https://github.com/gamblor21/CircuitPython_AHRS.git"

#============================================================================================
# Functions

#-------------------------------------------------------------------------------------------
# AHRS algorithm update

class mahony(object):
    def __init__(self):
        self.twoKp = 2.0 * 0.5  # 2 * proportional gain (Kp)
        self.twoKi = 2.0 * 0.0  # 2 * integral gain (Ki)
        self.q0 = 1.0
        self.q1 = 0.0
        self.q2 = 0.0
        self.q3 = 0.0
        self.integralFBx = 0.0
        self.integralFBy = 0.0
        self.integralFBz = 0.0
        self.anglesComputed = 0
        self.sample_freq = 100.0
        self.invSampleFreq = 1.0 / self.sample_freq

        self.roll = 0.0
        self.yaw = 0.0
        self.pitch = 0.0

    def invSqrt(self, x):
        return x ** -0.5

    def update(self, gx, gy, gz, ax, ay, az, mx, my, mz):
        recipNorm = 0
        q0q0 = q0q1 = q0q2 = q0q3 = q1q1 = q1q2 = q1q3 = q2q2 = q2q3 = q3q3 = 0
        hx = hy = bx = bz = 0
        halfvx = halfvy = halfvz = halfwx = halfwy = halfwz = 0
        halfex = halfey = halfez = 0
        qa = qb = qc = 0

        # Use IMU algorithm if magnetometer measurement invalid
        # (avoids NaN in magnetometer normalisation)
        if (mx == 0.0) and (my == 0.0) and (mz == 0.0):
            updateIMU(gx, gy, gz, ax, ay, az)
            return

        # Convert gyroscope degrees/sec to radians/sec
        gx *= 0.0174533
        gy *= 0.0174533
        gz *= 0.0174533

        # Compute feedback only if accelerometer measurement valid
        # (avoids NaN in accelerometer normalisation)
        if not ((ax == 0.0) and (ay == 0.0) and (az == 0.0)):
            # Normalise accelerometer measurement
            recipNorm = self.invSqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm

            # Normalise magnetometer measurement
            recipNorm = self.invSqrt(mx * mx + my * my + mz * mz)
            mx *= recipNorm
            my *= recipNorm
            mz *= recipNorm

            # Auxiliary variables to avoid repeated arithmetic
            q0q0 = self.q0 * self.q0
            q0q1  = self.q0 * self.q1
            q0q2  = self.q0 * self.q2
            q0q3  = self.q0 * self.q3
            q1q1  = self.q1  * self.q1
            q1q2  = self.q1  * self.q2
            q1q3  = self.q1  * self.q3
            q2q2  = self.q2  * self.q2
            q2q3  = self.q2  * self.q3
            q3q3  = self.q3  * self.q3

            # Reference direction of Earth's magnetic field
            hx = 2.0 * (mx * (0.5 - q2q2 - q3q3) + my * (q1q2 - q0q3) + mz * (q1q3 + q0q2))
            hy = 2.0 * (mx * (q1q2 + q0q3) + my * (0.5 - q1q1 - q3q3) + mz * (q2q3 - q0q1))
            bx = math.sqrt(hx * hx + hy * hy)
            bz = 2.0 * (mx * (q1q3 - q0q2) + my * (q2q3 + q0q1) + mz * (0.5 - q1q1 - q2q2))

            # Estimated direction of gravity and magnetic field
            halfvx = q1q3 - q0q2
            halfvy = q0q1 + q2q3
            halfvz = q0q0 - 0.5 + q3q3
            halfwx = bx * (0.5 - q2q2 - q3q3) + bz * (q1q3 - q0q2)
            halfwy = bx * (q1q2 - q0q3) + bz * (q0q1 + q2q3)
            halfwz = bx * (q0q2 + q1q3) + bz * (0.5 - q1q1 - q2q2)

            # Error is sum of cross product between estimated direction
            # and measured direction of field vectors
            halfex = (ay * halfvz - az * halfvy) + (my * halfwz - mz * halfwy)
            halfey = (az * halfvx - ax * halfvz) + (mz * halfwx - mx * halfwz)
            halfez = (ax * halfvy - ay * halfvx) + (mx * halfwy - my * halfwx)

            # Compute and apply integral feedback if enabled
            if self.twoKi > 0.0:
                # integral error scaled by Ki
                self.integralFBx += self.twoKi * halfex * self.invSampleFreq
                self.integralFBy += self.twoKi * halfey * self.invSampleFreq
                self.integralFBz += self.twoKi * halfez * self.invSampleFreq
                gx += self.integralFBx  # apply integral feedback
                gy += self.integralFBy
                gz += self.integralFBz
            else:
                self.integralFBx = 0.0  # prevent integral windup
                self.integralFBy = 0.0
                self.integralFBz = 0.0

            # Apply proportional feedback
            gx += self.twoKp * halfex
            gy += self.twoKp * halfey
            gz += self.twoKp * halfez

        # Integrate rate of change of quaternion
        gx *= (0.5 * self.invSampleFreq)  # pre-multiply common factors
        gy *= (0.5 * self.invSampleFreq)
        gz *= (0.5 * self.invSampleFreq)
        qa = self.q0
        qb = self.q1
        qc = self.q2
        self.q0 += (-qb * gx - qc * gy - self.q3 * gz)
        self.q1 += (qa * gx + qc * gz - self.q3 * gy)
        self.q2 += (qa * gy - qb * gz + self.q3 * gx)
        self.q3 += (qa * gz + qb * gy - qc * gx)

        # Normalise quaternion
        recipNorm = self.invSqrt(self.q0 * self.q0 + self.q1  * self.q1  + self.q2  * self.q2  + self.q3  * self.q3 )
        self.q0 *= recipNorm
        self.q1  *= recipNorm
        self.q2  *= recipNorm
        self.q3  *= recipNorm
        self.anglesComputed = 0

#-------------------------------------------------------------------------------------------
# IMU algorithm update

    def updateIMU(self, gx, gy, gz, ax, ay, az):
        recipNorm = 0
        halfvx = halfvy = halfvz =0
        halfex = halfey = halfez =0
        qa = qb = qc =0

        # Convert gyroscope degrees/sec to radians/sec
        gx *= 0.0174533
        gy *= 0.0174533
        gz *= 0.0174533

        # Compute feedback only if accelerometer measurement valid
        # (avoids NaN in accelerometer normalisation)
        if not ((ax == 0.0) and (ay == 0.0) and (az == 0.0)):
            # Normalise accelerometer measurement
            recipNorm = self.invSqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm

            # Estimated direction of gravity
            halfvx = self.q1  * self.q3  - self.q0 * self.q2
            halfvy = self.q0 * self.q1  + self.q2  * self.q3
            halfvz = self.q0 * self.q0 - 0.5 + self.q3  * self.q3

            # Error is sum of cross product between estimated
            # and measured direction of gravity
            halfex = (ay * halfvz - az * halfvy)
            halfey = (az * halfvx - ax * halfvz)
            halfez = (ax * halfvy - ay * halfvx)

            # Compute and apply integral feedback if enabled
            if self.twoKi > 0.0:
                # integral error scaled by Ki
                self.integralFBx += self.twoKi * halfex * self.invSampleFreq
                self.integralFBy += self.twoKi * halfey * self.invSampleFreq
                self.integralFBz += self.twoKi * halfez * self.invSampleFreq
                gx += integralFBx  # apply integral feedback
                gy += integralFBy
                gz += integralFBz
            else:
                self.integralFBx = 0.0  # prevent integral windup
                self.integralFBy = 0.0
                self.integralFBz = 0.0

            # Apply proportional feedback
            gx += self.twoKp * halfex
            gy += self.twoKp * halfey
            gz += self.twoKp * halfez

        # Integrate rate of change of quaternion
        gx *= (0.5 * self.invSampleFreq)  # pre-multiply common factors
        gy *= (0.5 * self.invSampleFreq)
        gz *= (0.5 * self.invSampleFreq)
        qa = self.q0
        qb = self.q1
        qc = self.q2
        self.q0 += (-qb * gx - qc * gy - self.q3  * gz)
        self.q1  += (qa * gx + qc * gz - self.q3  * gy)
        self.q2  += (qa * gy - qb * gz + self.q3  * gx)
        self.q3  += (qa * gz + qb * gy - qc * gx)

        # Normalise quaternion
        recipNorm = self.invSqrt(self.q0 * self.q0 + self.q1  * self.q1  + self.q2  * self.q2  + self.q3  * self.q3 )
        self.q0 *= recipNorm
        self.q1  *= recipNorm
        self.q2  *= recipNorm
        self.q3  *= recipNorm
        self.anglesComputed = 0

    def computeAngles(self):
        self.roll = math.atan2(self.q0 * self.q1  + self.q2  * self.q3 , 0.5 - self.q1  * self.q1  - self.q2  * self.q2 )
        self.pitch = math.asin(-2.0 * (self.q1  * self.q3 - self.q0 * self.q2))
        self.yaw = math.atan2(self.q1 * self.q2 + self.q0 * self.q3, 0.5 - self.q2 * self.q2 - self.q3 * self.q3)
        self.anglesComputed = 1
