Introduction
============

.. image:: https://readthedocs.org/projects/circuitpython-ahrs/badge/?version=latest
    :target: https://circuitpython-ahrs.readthedocs.io/
    :alt: Documentation Status

.. image:: https://img.shields.io/discord/327254708534116352.svg
    :target: https://discord.gg/nBQh6qu
    :alt: Discord

.. image:: https://github.com/gamblor21/CircuitPython_AHRS/workflows/Build%20CI/badge.svg
    :target: https://github.com/gamblor21/CircuitPython_AHRS/actions
    :alt: Build Status

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: Code Style: Black

AHRS library for CircuitPython

This library contains right now one alogrithm for AHRS - Attitude and Heading Reference System.
It is used to combine multiple sensor values to give a heading, pitch and roll value such as used
by aircraft.


Dependencies
=============
This driver depends on:

* `Adafruit CircuitPython <https://github.com/adafruit/circuitpython>`_

Please ensure all dependencies are available on the CircuitPython filesystem.
This is easily achieved by downloading
`the Adafruit library and driver bundle <https://circuitpython.org/libraries>`_.

Installing from PyPI
=====================
.. note:: This library is not available on PyPI yet. Install documentation is included
   as a standard element. Stay tuned for PyPI availability!

.. todo:: Remove the above note if PyPI version is/will be available at time of release.
   If the library is not planned for PyPI, remove the entire 'Installing from PyPI' section.

On supported GNU/Linux systems like the Raspberry Pi, you can install the driver locally `from
PyPI <https://pypi.org/project/adafruit-circuitpython-ahrs/>`_. To install for current user:

.. code-block:: shell

    pip3 install adafruit-circuitpython-ahrs

To install system-wide (this may be required in some cases):

.. code-block:: shell

    sudo pip3 install adafruit-circuitpython-ahrs

To install in a virtual environment in your current project:

.. code-block:: shell

    mkdir project-name && cd project-name
    python3 -m venv .env
    source .env/bin/activate
    pip3 install adafruit-circuitpython-ahrs

Usage Example
=============

.. todo:: Add a quick, simple example. It and other examples should live in the examples folder and be included in docs/examples.rst.

Contributing
============

Contributions are welcome! Please read our `Code of Conduct
<https://github.com/gamblor21/CircuitPython_AHRS/blob/master/CODE_OF_CONDUCT.md>`_
before contributing to help this project stay welcoming.

Documentation
=============


