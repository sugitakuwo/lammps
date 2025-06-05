.. index:: fix nve/sphere/sllod

fix nve/sphere/sllod command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID nve/sphere/sllod

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/sphere/sllod = style name of this fix command
* same optional keywords as :doc:`fix nve/sphere <fix_nve_sphere>` may be appended

Description
"""""""""""

Perform constant NVE integration for finite-size spherical particles
using the SLLOD equations of motion. The fix must be used together
with :doc:`fix deform <fix_deform>` employing ``remap v`` to impose the
shear flow. Particle positions, velocities, and angular velocities are
integrated without thermostatting.

Related commands
""""""""""""""""

:doc:`fix nve/sphere <fix_nve_sphere>`, :doc:`fix nvt/sllod <fix_nvt_sllod>`

Default
"""""""

none
