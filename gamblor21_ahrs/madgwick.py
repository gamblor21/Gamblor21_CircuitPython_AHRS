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
`madgwick`
================================================================================

AHRS library for CircuitPython
Madgwick Algorithm

Implementation of Madgwick's IMU and AHRS algorithms.
See: http:##www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/

From the x-io website "Open-source resources available on this website are
provided under the GNU General Public Licence unless an alternative licence
is provided in source."

Date			Author          Notes
29/09/2011	SOH Madgwick    Initial release
02/10/2011	SOH Madgwick	Optimised for reduced CPU load
19/02/2012	SOH Madgwick	Magnetometer measurement is normalised

This version based upon AdaFruit AHRS https://github.com/adafruit/Adafruit_AHRS

* Author(s): Mark Komus

Implementation Notes
--------------------

**Hardware:**

    Any 9DOF sensor

**Software and Dependencies:**

* Adafruit CircuitPython firmware for the supported boards:
  https://github.com/adafruit/circuitpython/releases

"""

# imports
import math


class Madgwick:
    """AHRS Madgwick algorithm.
    """

    def __init__(self, beta=0.1, sample_freq=100):
        self._beta = beta
        self.q0 = 1.0
        self.q1 = 0.0
        self.q2 = 0.0
        self.q3 = 0.0

        self._sample_freq = sample_freq
        self.invSampleFreq = 1.0 / self.sample_freq

        self._roll = 0.0
        self._yaw = 0.0
        self._pitch = 0.0
        self._anglesComputed = False

    def _inv_sqrt(self, x):
        return x ** -0.5

    def update(self, gx, gy, gz, ax, ay, az, mx, my, mz):
        """Call this function sample_freq times a second with values from your sensor
        The units of the accelerometer and magnetometer do not matter for this alogirthm
        The gryoscope must be in degrees/sec

        :param gx, gy, gz: Gyroscope values in degrees/sec
        :param ax, ay, az: Accelerometer values
        :param mx, my, mz: Magnetometer values
        """
        recipNorm = 0.0
        s0 = s1 = s2 = s3 = 0.0
        qDot1 = qDot2 = qDot3 = qDot4 = 0.0
        hx = hy = 0.0
        _2q0mx = (
            _2q0my
        ) = (
            _2q0mz
        ) = (
            _2q1mx
        ) = (
            _2bx
        ) = (
            _2bz
        ) = (
            _4bx
        ) = (
            _4bz
        ) = (
            _2q0
        ) = (
            _2q1
        ) = (
            _2q2
        ) = (
            _2q3
        ) = (
            _2q0q2
        ) = (
            _2q2q3
        ) = q0q0 = q0q1 = q0q2 = q0q3 = q1q1 = q1q2 = q1q3 = q2q2 = q2q3 = q3q3 = 0.0

        # Use IMU algorithm if magnetometer measurement invalid (avoids NaN in
        # magnetometer normalisation)
        if (mx == 0.0) and (my == 0.0) and (mz == 0.0):
            self.updateIMU(gx, gy, gz, ax, ay, az)
            return

        ## Convert gyroscope degrees/sec to radians/sec
        gx *= 0.0174533
        gy *= 0.0174533
        gz *= 0.0174533

        ## Rate of change of quaternion from gyroscope
        qDot1 = 0.5 * (-self.q1 * gx - self.q2 * gy - self.q3 * gz)
        qDot2 = 0.5 * (self.q0 * gx + self.q2 * gz - self.q3 * gy)
        qDot3 = 0.5 * (self.q0 * gy - self.q1 * gz + self.q3 * gx)
        qDot4 = 0.5 * (self.q0 * gz + self.q1 * gy - self.q2 * gx)

        ## Compute feedback only if accelerometer measurement valid (avoids NaN in
        ## accelerometer normalisation)
        if not ((ax == 0.0) and (ay == 0.0) and (az == 0.0)):
            ## Normalise accelerometer measurement
            recipNorm = self._inv_sqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm

            ## Normalise magnetometer measurement
            recipNorm = self._inv_sqrt(mx * mx + my * my + mz * mz)
            mx *= recipNorm
            my *= recipNorm
            mz *= recipNorm

            ## Auxiliary variables to avoid repeated arithmetic
            _2q0mx = 2.0 * self.q0 * mx
            _2q0my = 2.0 * self.q0 * my
            _2q0mz = 2.0 * self.q0 * mz
            _2q1mx = 2.0 * self.q1 * mx
            _2q0 = 2.0 * self.q0
            _2q1 = 2.0 * self.q1
            _2q2 = 2.0 * self.q2
            _2q3 = 2.0 * self.q3
            _2q0q2 = 2.0 * self.q0 * self.q2
            _2q2q3 = 2.0 * self.q2 * self.q3
            q0q0 = self.q0 * self.q0
            q0q1 = self.q0 * self.q1
            q0q2 = self.q0 * self.q2
            q0q3 = self.q0 * self.q3
            q1q1 = self.q1 * self.q1
            q1q2 = self.q1 * self.q2
            q1q3 = self.q1 * self.q3
            q2q2 = self.q2 * self.q2
            q2q3 = self.q2 * self.q3
            q3q3 = self.q3 * self.q3

            ## Reference direction of Earth's magnetic field
            hx = (
                mx * q0q0
                - _2q0my * self.q3
                + _2q0mz * self.q2
                + mx * q1q1
                + _2q1 * my * self.q2
                + _2q1 * mz * self.q3
                - mx * q2q2
                - mx * q3q3
            )
            hy = (
                _2q0mx * self.q3
                + my * q0q0
                - _2q0mz * self.q1
                + _2q1mx * self.q2
                - my * q1q1
                + my * q2q2
                + _2q2 * mz * self.q3
                - my * q3q3
            )
            _2bx = math.sqrt(hx * hx + hy * hy)
            _2bz = (
                -_2q0mx * self.q2
                + _2q0my * self.q1
                + mz * q0q0
                + _2q1mx * self.q3
                - mz * q1q1
                + _2q2 * my * self.q3
                - mz * q2q2
                + mz * q3q3
            )
            _4bx = 2.0 * _2bx
            _4bz = 2.0 * _2bz

            ## Gradient decent algorithm corrective step
            s0 = (
                -_2q2 * (2.0 * q1q3 - _2q0q2 - ax)
                + _2q1 * (2.0 * q0q1 + _2q2q3 - ay)
                - _2bz
                * self.q2
                * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx)
                + (-_2bx * self.q3 + _2bz * self.q1)
                * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my)
                + _2bx
                * self.q2
                * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz)
            )
            s1 = (
                _2q3 * (2.0 * q1q3 - _2q0q2 - ax)
                + _2q0 * (2.0 * q0q1 + _2q2q3 - ay)
                - 4.0 * self.q1 * (1 - 2.0 * q1q1 - 2.0 * q2q2 - az)
                + _2bz
                * self.q3
                * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx)
                + (_2bx * self.q2 + _2bz * self.q0)
                * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my)
                + (_2bx * self.q3 - _4bz * self.q1)
                * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz)
            )
            s2 = (
                -_2q0 * (2.0 * q1q3 - _2q0q2 - ax)
                + _2q3 * (2.0 * q0q1 + _2q2q3 - ay)
                - 4.0 * self.q2 * (1 - 2.0 * q1q1 - 2.0 * q2q2 - az)
                + (-_4bx * self.q2 - _2bz * self.q0)
                * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx)
                + (_2bx * self.q1 + _2bz * self.q3)
                * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my)
                + (_2bx * self.q0 - _4bz * self.q2)
                * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz)
            )
            s3 = (
                _2q1 * (2.0 * q1q3 - _2q0q2 - ax)
                + _2q2 * (2.0 * q0q1 + _2q2q3 - ay)
                + (-_4bx * self.q3 + _2bz * self.q1)
                * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx)
                + (-_2bx * self.q0 + _2bz * self.q2)
                * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my)
                + _2bx
                * self.q1
                * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz)
            )
            recipNorm = self._inv_sqrt(
                s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3
            )  ## normalise step magnitude
            s0 *= recipNorm
            s1 *= recipNorm
            s2 *= recipNorm
            s3 *= recipNorm

            ## Apply feedback step
            qDot1 -= self._beta * s0
            qDot2 -= self._beta * s1
            qDot3 -= self._beta * s2
            qDot4 -= self._beta * s3

        ## Integrate rate of change of quaternion to yield quaternion
        self.q0 += qDot1 * self.invSampleFreq
        self.q1 += qDot2 * self.invSampleFreq
        self.q2 += qDot3 * self.invSampleFreq
        self.q3 += qDot4 * self.invSampleFreq

        ## Normalise quaternion
        recipNorm = self._inv_sqrt(
            self.q0 * self.q0
            + self.q1 * self.q1
            + self.q2 * self.q2
            + self.q3 * self.q3
        )
        self.q0 *= recipNorm
        self.q1 *= recipNorm
        self.q2 *= recipNorm
        self.q3 *= recipNorm

        self._anglesComputed = False

    ##-------------------------------------------------------------------------------------------
    ## IMU algorithm update

    def updateIMU(self, gx, gy, gz, ax, ay, az):
        """
        Called is was have no mag reading (internal use)
        """
        recipNorm = 0.0
        s0 = s1 = s2 = s3 = 0.0
        qDot1 = qDot2 = qDot3 = qDot4 = 0.0
        _2q0 = (
            _2q1
        ) = (
            _2q2
        ) = _2q3 = _4q0 = _4q1 = _4q2 = _8q1 = _8q2 = q0q0 = q1q1 = q2q2 = q3q3 = 0.0

        ## Convert gyroscope degrees/sec to radians/sec
        gx *= 0.0174533
        gy *= 0.0174533
        gz *= 0.0174533

        ## Rate of change of quaternion from gyroscope
        qDot1 = 0.5 * (-self.q1 * gx - self.q2 * gy - self.q3 * gz)
        qDot2 = 0.5 * (self.q0 * gx + self.q2 * gz - self.q3 * gy)
        qDot3 = 0.5 * (self.q0 * gy - self.q1 * gz + self.q3 * gx)
        qDot4 = 0.5 * (self.q0 * gz + self.q1 * gy - self.q2 * gx)

        ## Compute feedback only if accelerometer measurement valid (avoids NaN in
        ## accelerometer normalisation)
        if not ((ax == 0.0) and (ay == 0.0) and (az == 0.0)):
            ## Normalise accelerometer measurement
            recipNorm = self._inv_sqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm

            ## Auxiliary variables to avoid repeated arithmetic
            _2q0 = 2.0 * self.q0
            _2q1 = 2.0 * self.q1
            _2q2 = 2.0 * self.q2
            _2q3 = 2.0 * self.q3
            _4q0 = 4.0 * self.q0
            _4q1 = 4.0 * self.q1
            _4q2 = 4.0 * self.q2
            _8q1 = 8.0 * self.q1
            _8q2 = 8.0 * self.q2
            q0q0 = self.q0 * self.q0
            q1q1 = self.q1 * self.q1
            q2q2 = self.q2 * self.q2
            q3q3 = self.q3 * self.q3

            ## Gradient decent algorithm corrective step
            s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay
            s1 = (
                _4q1 * q3q3
                - _2q3 * ax
                + 4.0 * q0q0 * self.q1
                - _2q0 * ay
                - _4q1
                + _8q1 * q1q1
                + _8q1 * q2q2
                + _4q1 * az
            )
            s2 = (
                4.0 * q0q0 * self.q2
                + _2q0 * ax
                + _4q2 * q3q3
                - _2q3 * ay
                - _4q2
                + _8q2 * q1q1
                + _8q2 * q2q2
                + _4q2 * az
            )
            s3 = 4.0 * q1q1 * self.q3 - _2q1 * ax + 4.0 * q2q2 * self.q3 - _2q2 * ay
            recipNorm = self._inv_sqrt(
                s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3
            )  ## normalise step magnitude
            s0 *= recipNorm
            s1 *= recipNorm
            s2 *= recipNorm
            s3 *= recipNorm

            ## Apply feedback step
            qDot1 -= self._beta * s0
            qDot2 -= self._beta * s1
            qDot3 -= self._beta * s2
            qDot4 -= self._beta * s3

        ## Integrate rate of change of quaternion to yield quaternion
        self.q0 += qDot1 * self.invSampleFreq
        self.q1 += qDot2 * self.invSampleFreq
        self.q2 += qDot3 * self.invSampleFreq
        self.q3 += qDot4 * self.invSampleFreq

        ## Normalise quaternion
        recipNorm = self._inv_sqrt(
            self.q0 * self.q0
            + self.q1 * self.q1
            + self.q2 * self.q2
            + self.q3 * self.q3
        )
        self.q0 *= recipNorm
        self.q1 *= recipNorm
        self.q2 *= recipNorm
        self.q3 *= recipNorm

        self._anglesComputed = False

    def compute_angles(self):
        """
        Compute all the angles if there have been new samples (internal use)
        """
        self._roll = math.atan2(
            self.q0 * self.q1 + self.q2 * self.q3,
            0.5 - self.q1 * self.q1 - self.q2 * self.q2,
        )
        self._pitch = math.asin(-2.0 * (self.q1 * self.q3 - self.q0 * self.q2))
        self._yaw = math.atan2(
            self.q1 * self.q2 + self.q0 * self.q3,
            0.5 - self.q2 * self.q2 - self.q3 * self.q3,
        )
        self._anglesComputed = True

    @property
    def yaw(self):
        """
        Current yaw (z-axis) value in radians/sec. (read-only)
        """
        if not self._anglesComputed:
            self.compute_angles()
        return self._yaw

    @property
    def pitch(self):
        """
        Current pitch (y-axis) value in radians/sec. (read-only)
        """
        if not self._anglesComputed:
            self.compute_angles()
        return self._pitch

    @property
    def roll(self):
        """
        Current roll (x-axis) value in radians/sec. (read-only)
        """
        if not self._anglesComputed:
            self.compute_angles()
        return self._roll

    @property
    def sample_freq(self):
        """The current sample frequency value in Hertz."""
        return self._sample_freq

    @sample_freq.setter
    def sample_freq(self, value):
        self._sample_freq = value
        self.invSampleFreq = 1.0 / self.sample_freq

    @property
    def beta(self):
        """The current beta value (Proportional gain)."""
        return self._beta

    @beta.setter
    def beta(self, value):
        self._beta = value
