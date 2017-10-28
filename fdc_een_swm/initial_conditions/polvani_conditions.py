#!/usr/bin/env/ python3
"""
This module contains the necessary objects to create and initially-balanced
turbulent flow field and to analyse the statistics of that field.

Classes
-------
- `TurbulentICs` : Creates the appropriate flow field

Miscellaneous Functions
-----------------------
- `build_initial_condition` : Creates the turbulent, initially-balanced condition

| Author: Adam G. Peddle
| Contact: ap553@exeter.ac.uk
| Version: 1.0
"""

import logging
import sys
import getopt
import random
import numpy as np
from spectral_toolbox import *

class TurbulentICs:
    """
    This class encapsulates the necessary methods to create initially-balanced turbulent
    initial conditions, following Polvani et al 1994.

    **Parameters**

    - `control` : a dictionary containing the necessary parameters to set the attributes
    - `toolbox` : a toolbox of methods for working with spectral methods

    **Attributes**

    - `N` : The domain size, such that the domain consists of N X N equally-spaced gridpoints
    - `L` : The dimensional sidelength of the square domain
    - `k0` : The peak location of the initial kinetic energy spectrum
    - `m` : The power used in creating the spectrum (see energy_spectrum)
    - `outname` : The stem to be used in naming any output files created
    - `psi` : An array containing the streamfunction in spectral space
    - `chi` : An array containing the velocity potential in spectral space
    - `h` : The perturbation height of the fluid
    - `vx`, `vx` : The x- and y-velocities of the fluid flow

    - `st` : The spectral toolbox used to handle the standard spectral techniques

    **Methods**

    - `energy_spectrum` : compute desired kinetic energy in shell
    - `construct_streamfunction` : create the streamfn, psi
    - `create_height` : create the initial height
    - `jacobian_t` : compute the jacobian_sub_t
    - `update_psi_t` : compute the linearised dpsi/dt
    - `create_velocity_potential` : create chi by iteration
    - `compute_velocities` : compute velocities by Helmholtz
    - `check_potential` : check of accurate chi iteration
    - `check_height` : check of accuracy in height field

    """
    def __init__(self, control, toolbox):
        self.N = int(control['Nx'])
        self.L = control['Lx']
        self.k0 = 14.0
        self.m = 25.0
        self.outname = control['outname']

        self.psi = np.zeros((control['Nx'], control['Nx']), dtype = complex) # Streamfunction
        self.chi = np.zeros((control['Nx'], control['Nx']), dtype = complex) # Velocity potential

        self.h = np.zeros((control['Nx'], control['Nx']), dtype = complex)
        self.vx = np.zeros((control['Nx'], control['Nx']), dtype = complex)
        self.vy = np.zeros((control['Nx'], control['Nx']), dtype = complex)

        self.st = toolbox

    def energy_spectrum(self, k1, k2):
        """
        This method implements the shape of the initial energy spectrum:

        .. math:: E(k) \\approx \\frac{k^{m/2}}{(k + k_{0})^{m}}

        **Parameters**

        - `k1`, `k2` : The wavenumbers in the x- and y-directions

        **Returns**

        - `E` : the energy in the desired shell

        """
        # Compute the modulus of the wavenumbers
        k_mod = np.sqrt(k1**2 + k2**2)

        # Normalise the energy to yield initial RMS values of 1/sqrt{2}
        norm_val = 1.045456659e21*2.0*np.pi*self.k0*19.45
        return norm_val*k_mod**(self.m/2)/(k_mod + self.k0)**self.m

    def construct_streamfunction(self):
        """
        Construct the initially-balanced streamfunction.

        Sets the attribute psi with the streamfunction computed according to

        .. math:: 1/2 k^{2}|\\psi_{k}|^{2} = E_{K}(k)

        Proceeds by creating a vorticity field with initially random phases
        in Fourier space and moduli of 1. This field is then integrated to
        find the corresponding streamfunction which is 'tuned' to have the
        desired kinetic energy spectrum. The spectrum is tuned such that
        the desired spectrum occurs in shells.

        """

        # Initialise the random phases, e^{i*theta}
        N = self.N
        phases = np.random.uniform(0.0, 2*np.pi, (N-1,N-1))
        vorticity_hat = np.zeros((N,N),dtype=complex)
        vorticity_buffer = np.exp(1j*phases)

        Nyquist_1 = np.exp(1j*np.random.uniform(0.0,2*np.pi, (N//2-1)))
        Nyquist_2 = np.exp(1j*np.random.uniform(0.0,2*np.pi, (N//2-1)))

        # Ensure that the initial vorticity field has a Hermitian
        # spectrum in Fourier space for a real-valued solution.
        mid = (N-1)//2
        mid2 = 2*mid
        for k1 in range(N-1):
            for k2 in range(k1 + 1, N-1):
                vorticity_buffer[k1, k2] = np.conj(vorticity_buffer[mid2 - k1,mid2 - k2])

        for k in range(mid):
            vorticity_buffer[k,k] = np.conj(vorticity_buffer[mid2 - k, mid2 - k])

        vorticity_hat[1:,1:] = vorticity_buffer

        # Set up the Nyqusist frequencies to yield a real-valued,
        # C_{inf} vorticity field in realspace.
        vorticity_hat[0,1:N//2] = Nyquist_1
        vorticity_hat[0,N//2+1:N] = np.conj(Nyquist_1[::-1])
        vorticity_hat[1:N//2,0] = Nyquist_2
        vorticity_hat[N//2+1:N,0] = np.conj(Nyquist_2[::-1])

        vorticity_hat[N//2, N//2] = 0.0
        vorticity_hat[0,0] = 0.0
        vorticity_hat[0,N//2] = 0.0
        vorticity_hat[N//2,0] = 0.0

        vorticity = self.st.inverse_fft(vorticity_hat)

        vorticity/=np.amax(np.real(vorticity))
        vorticity_hat = self.st.forward_fft(-vorticity)

        # Find the associated streamfunction via inverse Helmholtz
        self.psi = self.st.solve_inverse_laplacian(vorticity_hat, 0.0)

        # Calculate the total kinetic energy in each wavenumber shell
        E_k = dict()
        for k1 in range(-N//2, N//2):
            for k2 in range(-N//2, N//2):
                k = np.sqrt(k1**2 + k2**2)
                k_ctr = int(round(k))

                psi_value = self.psi[k1 + N//2, k2 + N//2]
                local_energy = np.real(0.5*((2*np.pi)**2)*k**2*psi_value*np.conj(psi_value))
                try:
                    E_k[k_ctr] += local_energy
                except KeyError:
                    E_k[k_ctr] = local_energy

        for k1 in range(-N//2, N//2):
            for k2 in range(-N//2, N//2):
                k = np.sqrt(k1**2 + k2**2)
                k_ctr = int(round(k))
                # Compute the desired kinetic energy level for the
                # current shell
                E_target = self.energy_spectrum(k1, k2)

                # `Tune' the kinetic energy spectrum, shell by shell
                if k_ctr > 0 and E_k[k_ctr] > 1e-12:
                    self.psi[k1 + N//2, k2 + N//2] *= np.sqrt(E_target/E_k[k_ctr])
                else:
                    self.psi[k1 + N//2, k2 + N//2] = 0.0

    def create_height(self, control):
        """
        Compute the initially-balanced nondimensional perturbation height field

        The initially-balanced perturbation height is computed by inverting the
        LHS of the following equation:

        .. math:: \\nabla^{2}h = \\nabla^{2}\\psi + 2RJ(\\psi_{x},\\psi_{y})

        **Parameters**

        - `control` : dictionary containing necessary parameters

        **Returns**


        - `h` : nondimensional perturbation height (set in object)

        """
        # Compute Right-Hand Side
        RHS = self.st.laplacian(self.psi) + \
              2.0*control['Rossby']*self.st.jacobian1(self.psi)

        self.h = self.st.solve_inverse_laplacian(RHS, 0)

    def jacobian_t(self, psi, psi_t):
        """
        Compute the 'Jacobian-sub-t'.

        Computes the Jacobian-sub-t used by Polvani et al (1994):

        .. math:: J(\\psi_{x},\\psi_{y})_{t} = (\\psi_{t})_{xx}\\psi_{yy} +
                                               \\psi_{xx}(\\psi_{t})_{yy} -
                                               2\\psi_{xy}(\\psi_{t})_{xy}

        **Parameters**

        - `psi` : the streamfunction at the current iteration
        - `psi_t` : the value of the linearised time derivative of the streamfunction

        **Returns**

        - `J` : the Fourier representation of the Jacobian-sub-t

        """

        # Compute the spatial derivatives of psi_t
        psi_t_xx_hat = self.st.calc_derivative(psi_t, 'x', 'x')
        psi_t_yy_hat = self.st.calc_derivative(psi_t, 'y', 'y')
        psi_t_xy_hat = self.st.calc_derivative(psi_t, 'x', 'y')

        # Compute the spatial derivatives of psi
        psi_xx_hat = self.st.calc_derivative(psi, 'x', 'x')
        psi_yy_hat = self.st.calc_derivative(psi, 'y', 'y')
        psi_xy_hat = self.st.calc_derivative(psi, 'x', 'y')

        # Transform to real-space for computation of nonlinear
        # multiplication
        psi_t_xx = self.st.inverse_fft(self.st.dealias_pad(psi_t_xx_hat))
        psi_t_yy = self.st.inverse_fft(self.st.dealias_pad(psi_t_yy_hat))
        psi_t_xy = self.st.inverse_fft(self.st.dealias_pad(psi_t_xy_hat))

        psi_xx = self.st.inverse_fft(self.st.dealias_pad(psi_xx_hat))
        psi_yy = self.st.inverse_fft(self.st.dealias_pad(psi_yy_hat))
        psi_xy = self.st.inverse_fft(self.st.dealias_pad(psi_xy_hat))

        # Computation of realspace representation of Jacobian
        J_canonical = (psi_t_xx*psi_yy) + (psi_t_yy*psi_xx) - 2.0*(psi_t_xy*psi_xy)
        return self.st.dealias_unpad(self.st.forward_fft(J_canonical))

    def update_psi_t(self, control):
        """
        Updates the locally linearised estimate of psi_t

        Updates the value of psi_t for use by the Jacobian-sub-t by
        inverting the following equation:

        .. math:: \\nabla^{2}\\psi_{t} = J(\\psi, \\nabla^{2}\\psi) + R^{-1}\\nabla^{2}\\chi + \\nabla \\cdot (\\nabla^{2}\\psi, \\nabla\\chi)

        **Parameters**

        - `control` : dictionary containing required parameters (Rossby)

        **Returns**

        - `psi_t` : updated value of psi_t

        """
        # Initialise psi_t
        psi_t = np.zeros((self.N, self.N), dtype = complex)

        # Computation of first two terms on right-hand side
        RHS = self.st.jacobian(self.psi, self.st.laplacian(self.psi)) + \
              control['Rossby']**-1*self.st.laplacian(self.chi)

        # Computation of last term on RHS (divergence term)
        term1 = self.st.multiply_nonlinear(self.st.laplacian(self.psi), self.st.calc_derivative(self.chi, 'x'))
        term2 = self.st.multiply_nonlinear(self.st.laplacian(self.psi), self.st.calc_derivative(self.chi, 'y'))
        RHS += self.st.calc_derivative(term1, 'x') + self.st.calc_derivative(term2, 'y')

        # Computation of psi_t by inverting LHS
        psi_t = self.st.solve_inverse_laplacian(RHS, 0)
        return psi_t

    def create_velocity_potential(self, control):
        """
        Create the initially-balanced turbulent velocity potential.

        Employs an iterative procedure to find the velocity potential, chi,
        which satisfies the following equation:

        .. math:: R^{-1}(1 - B\\nabla^{2})\\nabla^{2}\\chi = -J(\\psi, \\nabla^{2}\\chi) + \
                                                             \\nabla^{2}J(\\psi, h) + \
                                                             2RJ(\\psi_{x},\\psi_{y})_{t} - \
                                                             \\nabla \\cdot(\\nabla^{2}\\psi, \\nabla \\chi) + \
                                                             \\nabla^{2}[\\nabla \\cdot(h \\nabla \\chi)]

        The procedure iterates until the L_inf error between the current and previous values of chi are
        less than the tolerance specified in the control dictionary. The initial value of chi for the
        procedure is taken to be 0 everywhere.

        **Parameters**

        - `control` : the control dictionary, as before

        **Returns**

        - `chi` : the velocity potential (set in object, not explicitly returned)

        **See Also**

        update_psi_t, jacobian_t

        """
        # Initialise old value of chi for use in convergence checking
        chi_old = np.zeros((control['Nx'], control['Nx']), dtype = complex)

        L_inf = 1.0e10
        # Skip convergence check on first iteration to avoid
        # division-by-zero
        first_iteration_flag = True
        while L_inf > control['tolerance']:

            # Compute terms on the right hand side
            RHS = -self.st.jacobian(self.psi, self.st.laplacian(self.chi))
            RHS += self.st.laplacian(self.st.jacobian(self.psi, self.h))
            psi_t = self.update_psi_t(control)
            RHS += 2.0*control['Rossby']*self.jacobian_t(self.psi, psi_t)

            # Compute divergence terms on RHS as well
            term1 = self.st.multiply_nonlinear(self.st.laplacian(self.psi), self.st.calc_derivative(self.chi, 'x'))
            term2 = self.st.multiply_nonlinear(self.st.laplacian(self.psi), self.st.calc_derivative(self.chi, 'y'))
            RHS -= self.st.calc_derivative(term1, 'x') + self.st.calc_derivative(term2, 'y')

            term1 = self.st.multiply_nonlinear(self.h, self.st.calc_derivative(self.chi, 'x'))
            term2 = self.st.multiply_nonlinear(self.h, self.st.calc_derivative(self.chi, 'y'))
            buff_mat = self.st.calc_derivative(term1, 'x') + self.st.calc_derivative(term2, 'y')
            RHS += self.st.laplacian(buff_mat)

            # Compute the laplacian and hyperlaplacian on the LHS
            LHS = (self.st.laplace_op - control['Burger']*self.st.laplace_op*self.st.laplace_op)/control['Rossby']

            # Constant Fourier mode is a free choice. Set the corresponding LHS value to 1
            # to prevent divide-by-zero error...
            LHS[self.N//2, self.N//2] = 1.0
            self.chi = RHS/LHS

            # ...and set the constant mode of the potential to zero.
            self.chi[self.N//2, self.N//2] = 0.0

            chi_canonical = self.st.inverse_fft(self.chi)

            if not first_iteration_flag: # Skip convergence checking on first iteration
                # Compute L_inf between current and previous iteration:
                L_inf_mat = (chi_canonical - chi_old)/(chi_old)
                L_inf = np.max(np.absolute(L_inf_mat))
                print('L_inf: ', L_inf)
            else:
                first_iteration_flag = False
                print('First iteration complete')

            # Set old value of chi for convergence checking
            chi_old = chi_canonical.copy()

    def compute_velocities(self, control):
        """
        Compute the flow velocities from the streamfunction and potential.

        This method implements the velocity computations from the
        potential and streamfunction through the Helmholtz equation:

        .. math:: u = k \\times \\nabla \\psi  + \\epsilon \\nabla \\chi

        where epsilon is computed as:

        .. math:: \\epsilon = R\\max(1,R)/\\max(1,B)

        which leads to a uniform definition of balance when R or F are much
        smaller than unity.

        **Parameters**

        - `control` : the control dictionary

        **Returns**

        - `vx`, `vy` : the x- and y-velocities, set in the object

        """

        # Compute epsilon for ratio between divergent and rotational
        # components
        epsilon = control['Rossby']*max(1,control['Rossby'])/max(1, control['Burger'])

        # Compute x- and y- derivatives of the streamfunction
        psi_x = self.st.calc_derivative(self.psi, 'x')
        psi_y = self.st.calc_derivative(self.psi, 'y')

        # Compute x- and y- derivatives of the potential
        chi_x = self.st.calc_derivative(self.chi, 'x')
        chi_y = self.st.calc_derivative(self.chi, 'y')

        # Compute velocities by Helmholtz equation
        self.vx = -psi_y + epsilon*chi_x
        self.vy = psi_x + epsilon*chi_y


def build_initial_condition(rossby, froude, nx, outfilename, csvFlag):
    """
    Main method to compute initially-balanced turbulent flow field.

    Creates initially-balanced randomised flow field plus height field for
    use in rotating shallow water equations. Follows the procedure implemented
    in Polvani et al (1994).

    Relies on the methods contained in the TurbulentICs class as well as those in
    the spectral_toolbox.

    **Parameters**

    - `rossby` : the desired Rossby number, defined as :math: `R = U/(f_{0}L)`
    - `froude` : the Froude number, defined as :math: `F = U/\\sqrt{gH}`
    - `nx` : the number of points to use, such that the domain is n x n
    - `outfilename` : the stem to be used for creating outfiles
    - `csvFlag` : optional flag to create MATLAB-syle output

    **Returns**

    - `outfilename.bin` : data file for use by RSWE code
    - `outfilename_u.txt` : CSV file with v1
    - `outfilename_v.txt` : CSV file with v2
    - `outfilename_h.txt` : CSV file with h
    - `outfilename_data.txt` : CSV file with other data

    **See Also**

    TurbulentICs, SpectralToolbox

    """
    control = dict()
    np.set_printoptions(precision=2,linewidth=200)

    # Set control parameters:
    control['Rossby'] = rossby
    control['Froude'] = froude
    control['Nx'] = nx
    control['Lx'] = 2.0*np.pi*14.0
    control['outname'] = outfilename

    # Compute the Burger number
    control['Burger'] = (control['Rossby']/control['Froude'])**2

    # Hardcode tolerance for iterative solver
    control['tolerance'] = 1e-6
    logging.info('Control Structure Created Successfully')

    # Create spectral toolbox object for use in TurbulentICs
    st = SpectralToolbox(control['Nx'], control['Lx'])

    # Create turbulence handler
    turb_cond = TurbulentICs(control, st)
    logging.info('Turbulence Handler Created Successfully')

    # Create streamfunction field, store in turb_cond
    turb_cond.construct_streamfunction()
    logging.info('Streamfunction Created Successfully')

    # Create height field, store in turb_cond
    turb_cond.create_height(control)
    logging.info('Height Computed Successfully')

    # Create velocity potential field, store in turb_cond
    logging.info('Beginning Computation of Velocity Potential')
    turb_cond.create_velocity_potential(control)
    logging.info('Potential Created Successfully')

    # Compute velocities
    turb_cond.compute_velocities(control)

    # Save ICs to File
    if csvFlag:
        h = np.real(turb_cond.st.inverse_fft(turb_cond.h))
        h = 1.0 + control['Rossby']*h/control['Burger']

        v1 = np.real(turb_cond.st.inverse_fft(turb_cond.vx))
        v2 = np.real(turb_cond.st.inverse_fft(turb_cond.vy))

        np.savetxt(str(outfilename) + "_u.txt", v1, delimiter = "\t")
        np.savetxt(str(outfilename) + "_v.txt", v2, delimiter = "\t")
        np.savetxt(str(outfilename) + "_h.txt", h, delimiter = "\t")

        with open(str(outfilename) + '_data.txt', 'w') as f:
            f.write("Rossby {}\n".format(control['Rossby']))
            f.write("Froude {}\n".format(control['Froude']))
            f.write("Burger {}\n".format(control['Burger']))
            f.write("Nx {}\n".format(control['Nx']))

        logging.info("Initial Condition saved to successfully")

    else:
         with open(str(outfilename) + '.bin', 'wb') as f:
            np.save(f, control['Nx'])

            # Compute total nondimensional height
            height = 1.0 + (control['Rossby']/control['Burger'])*st.inverse_fft(turb_cond.h)
            np.save(f, st.forward_fft(height))
            np.save(f, turb_cond.vx)
            np.save(f, turb_cond.vy)

            logging.info("Initial Condition saved to {}".format(str(outfilename) + '.bin'))


if __name__ == "__main__":

    # Set up logging
    logging.basicConfig(level = logging.INFO)
    csvFlag = False

    # Get optional arguments to pass to initial condition builder
    opts, args = getopt.gnu_getopt(sys.argv[1:], '', ['csv'])
    for o, a in opts:
        if o in ("--csv"):
            csvFlag = True
        else:
            logging.warning('Unknown Option')

    # Call through to create intial condition
    build_initial_condition(float(args[0]), float(args[1]), int(args[2]), args[3], csvFlag)

